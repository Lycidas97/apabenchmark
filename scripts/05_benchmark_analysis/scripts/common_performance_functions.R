bedtools_path <- Sys.getenv("APABENCHMARK_BEDTOOLS_BIN", unset = "bedtools")

resolve_bedtools_bin <- function(bedtools_path) {
  if (dir.exists(bedtools_path)) {
    return(file.path(bedtools_path, "bedtools"))
  }
  if (file.exists(bedtools_path)) {
    return(bedtools_path)
  }

  resolved <- Sys.which(bedtools_path)
  if (nzchar(resolved)) {
    return(resolved)
  }
  stop(sprintf("bedtools binary not found: %s", bedtools_path))
}

safe_divide <- function(numerator, denominator) {
  if (is.na(numerator) || is.na(denominator) || denominator <= 0) {
    return(0)
  }
  return(numerator / denominator)
}

safe_f1 <- function(precision, recall) {
  return(safe_divide(2 * precision * recall, precision + recall))
}

read_table_if_exists <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.table())
  }
  return(fread(path))
}

read_count_matrix_dt <- function(path, first_col_name = "cell_id") {
  if (!file.exists(path)) {
    stop(sprintf("Matrix file not found: %s", path))
  }
  header_line <- readLines(path, n = 1, warn = FALSE)
  if (length(header_line) == 0) {
    stop(sprintf("Empty matrix file: %s", path))
  }

  header_fields <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  dt <- fread(path, header = FALSE, skip = 1, data.table = TRUE)

  if (ncol(dt) == length(header_fields) + 1) {
    header_fields <- c(first_col_name, header_fields)
  } else if (ncol(dt) != length(header_fields)) {
    stop(sprintf("Header/data column mismatch for %s: header=%d data=%d", path, length(header_fields), ncol(dt)))
  }

  setnames(dt, header_fields)
  return(dt)
}

rename_matched_columns <- function(match_dt, gt_col_num, res_col_num) {
  column_names <- c(
    "gt_chr", "gt_start", "gt_end", "gt_name", "gt_score", "gt_strand",
    "res_chr", "res_start", "res_end", "res_name", "res_score", "res_strand"
  )

  if (gt_col_num > 6) {
    for (i in 6:(gt_col_num - 1)) {
      column_names <- c(column_names[1:i], paste0("gt_", i + 1), column_names[(i + 1):length(column_names)])
    }
  }

  if (res_col_num > 6) {
    for (i in 6:(res_col_num - 1)) {
      column_names <- c(column_names, paste0("res_", i + 1))
    }
  }
  setnames(match_dt, old = names(match_dt), new = column_names)
}

calculate_index <- function(dt, index_type) {
  pos <- sapply(colnames(dt), function(col) {
    parts <- unlist(strsplit(col, split = "\\|"))
    if (length(parts) > 1) {
      pos_parts <- unlist(strsplit(parts[2], split = ":"))
      as.integer(pos_parts[2])
    } else {
      0
    }
  })

  strand <- sapply(colnames(dt), function(col) {
    parts <- unlist(strsplit(col, split = "\\|"))
    if (length(parts) > 1) {
      pos_parts <- unlist(strsplit(parts[2], split = ":"))
      pos_parts[length(pos_parts)]
    } else {
      "+"
    }
  })

  if (all(strand == "+")) {
    pos <- as.numeric(abs(pos - min(pos)))
  } else if (all(strand == "-")) {
    pos <- as.numeric(abs(pos - max(pos)))
  } else {
    stop("Mixed strands within the same data.table")
  }
  if (index_type == "PDUI") {
    index_weight <- as.numeric(pos == max(pos))
  } else if (index_type == "PPUI") {
    index_weight <- as.numeric(pos == min(pos))
  } else if (index_type == "RWUI") {
    index_weight <- as.numeric(order(pos)) - 1
    index_weight <- index_weight / max(index_weight)
  } else if (index_type == "DWUI") {
    index_weight <- pos / sum(pos)
  } else {
    stop("Index type not found, please check")
  }

  dt_matrix <- matrix(as.numeric(as.matrix(dt)), nrow = nrow(dt), ncol = ncol(dt))
  index_weight_matrix <- matrix(index_weight, nrow = nrow(dt), ncol = ncol(dt), byrow = TRUE)
  weighted_sum <- rowSums(dt_matrix * index_weight_matrix)
  row_sums <- rowSums(dt_matrix)
  index <- weighted_sum / row_sums
  names(index) <- rownames(dt)
  return(index)
}

wilcox_test_by_cell_type <- function(index_vec, cell_type_map) {
  names(index_vec) <- cell_type_map[names(index_vec)]
  index_vec <- index_vec[!is.na(index_vec)]
  index_vec <- split(index_vec, names(index_vec))
  if (length(index_vec) > 1) {
    return(wilcox.test(as.numeric(index_vec[[1]]), as.numeric(index_vec[[2]]))$p.value)
  }
  return(NA)
}

build_pseudobulk_by_split <- function(cell_feature_matrix, split_num, cell_type_label) {
  if (!is.numeric(split_num) || length(split_num) != 1 || split_num < 1 || split_num != as.integer(split_num)) {
    stop(sprintf("Invalid dexseq_split_num: %s", as.character(split_num)))
  }
  split_num <- as.integer(split_num)

  n_cells <- nrow(cell_feature_matrix)
  if (n_cells < split_num) {
    stop(sprintf("dexseq_split_num=%d requires at least %d cells for cell type %s, but got %d", split_num, split_num, cell_type_label, n_cells))
  }

  group_ids <- rep(seq_len(split_num), length.out = n_cells)
  idx_list <- split(seq_len(n_cells), group_ids)
  agg <- vapply(
    idx_list,
    function(idx) colSums(cell_feature_matrix[idx, , drop = FALSE]),
    FUN.VALUE = numeric(ncol(cell_feature_matrix))
  )
  agg <- as.matrix(agg)
  rownames(agg) <- colnames(cell_feature_matrix)
  colnames(agg) <- paste0(cell_type_label, "_", seq_len(split_num))
  return(agg)
}

run_dexseq_pipeline <- function(dxd, fit_type = NULL) {
  dxd <- DEXSeq::estimateSizeFactors(dxd, locfunc = genefilter::shorth)
  disp_args <- list(object = dxd)
  if (!is.null(fit_type)) {
    disp_args$fitType <- fit_type
  }
  dxd <- do.call(DEXSeq::estimateDispersions, disp_args)
  dxd <- DEXSeq::testForDEU(dxd)
  dxd <- DEXSeq::estimateExonFoldChanges(dxd)
  dxr <- as.data.table(DEXSeq::DEXSeqResults(dxd))
  return(dxr)
}

standardize_dexseq_result <- function(dxr) {
  if (!is.data.table(dxr)) {
    dxr <- as.data.table(dxr)
  }
  if (!("groupID" %in% names(dxr))) {
    dxr[, groupID := NA_character_]
  }
  if (!("padj" %in% names(dxr))) {
    dxr[, padj := NA_real_]
  }

  if (!("log2fold_T2_T1" %in% names(dxr))) {
    log2fc_cols <- grep("^log2fold_", names(dxr), value = TRUE)
    if (length(log2fc_cols) >= 1) {
      if (length(log2fc_cols) > 1) {
        warning(
          sprintf("Multiple DEXSeq log2fold columns detected (%s); using %s", paste(log2fc_cols, collapse = ", "), log2fc_cols[1]),
          call. = FALSE
        )
      }
      setnames(dxr, log2fc_cols[1], "log2fold_T2_T1")
    } else {
      dxr[, log2fold_T2_T1 := NA_real_]
    }
  }

  return(dxr[, .(groupID, padj, log2fold_T2_T1)])
}

dexseq_test_by_cell_type <- function(input_mtx, cell_type_map, split_num = 6) {
  mtx_for_dex <- mapply(function(mtx, te, cell_type_map) {
    te <- gsub(":", "_", te)
    colnames(mtx) <- paste0(te, ":", seq_along(colnames(mtx)))
    rownames(mtx) <- cell_type_map[rownames(mtx)]
    mtx <- matrix(as.numeric(as.matrix(mtx)), nrow = nrow(mtx), ncol = ncol(mtx), dimnames = list(rownames(mtx), colnames(mtx)))
  },
  input_mtx,
  names(input_mtx),
  MoreArgs = list(cell_type_map = cell_type_map),
  SIMPLIFY = FALSE) %>% do.call(cbind, .)

  cell_types <- unique(rownames(mtx_for_dex))
  if (length(cell_types) != 2) {
    stop(sprintf("DEXSeq requires exactly 2 cell types, got %d", length(cell_types)))
  }

  ct_dex_1 <- build_pseudobulk_by_split(
    cell_feature_matrix = mtx_for_dex[rownames(mtx_for_dex) == cell_types[1], , drop = FALSE],
    split_num = split_num,
    cell_type_label = cell_types[1]
  )
  ct_dex_2 <- build_pseudobulk_by_split(
    cell_feature_matrix = mtx_for_dex[rownames(mtx_for_dex) == cell_types[2], , drop = FALSE],
    split_num = split_num,
    cell_type_label = cell_types[2]
  )

  dex_mtx <- cbind(ct_dex_1, ct_dex_2)
  sample_table <- data.frame(
    row.names = colnames(dex_mtx),
    condition = factor(c(rep(cell_types[1], split_num), rep(cell_types[2], split_num)), levels = cell_types)
  )

  feature_id <- rownames(dex_mtx)
  group_id <- rownames(dex_mtx) %>% strsplit(., ":") %>% sapply(., function(x) {
    x[1]
  })

  dxd <- DEXSeqDataSet(
    countData = dex_mtx,
    sampleData = sample_table,
    design = ~ sample + exon + condition:exon,
    featureID = feature_id,
    groupID = group_id,
    featureRanges = NULL,
    transcripts = NULL,
    alternativeCountData = NULL
  )
  attempts <- list(
    list(label = "default", fit_type = NULL),
    list(label = "mean", fit_type = "mean")
  )
  attempt_errors <- character(0)

  for (attempt in attempts) {
    dxr <- tryCatch({
      run_dexseq_pipeline(dxd, fit_type = attempt$fit_type)
    }, error = function(e) {
      attempt_errors <<- c(attempt_errors, sprintf("%s: %s", attempt$label, conditionMessage(e)))
      NULL
    })

    if (!is.null(dxr)) {
      dxr <- standardize_dexseq_result(dxr)
      if (attempt$label != "default") {
        warning(sprintf("DEXSeq fallback succeeded with fitType='%s'", attempt$fit_type), call. = FALSE)
      }
      return(dxr)
    }
  }

  warning(
    sprintf("DEXSeq failed after fallback attempts; returning empty result. %s", paste(attempt_errors, collapse = " | ")),
    call. = FALSE
  )
  dxr_empty <- data.table(
    groupID = character(),
    padj = numeric(),
    log2fold_T2_T1 = numeric()
  )
  return(dxr_empty)
}

fisher_test_by_cell_type <- function(mtx, cell_type_map, simulate_p_value = TRUE, fisher_B = 2000) {
  if (ncol(mtx) <= 1) {
    return(NA_real_)
  }
  row_name <- as.vector(cell_type_map[rownames(mtx)])
  keep_rows <- !is.na(row_name)
  if (!any(keep_rows)) {
    return(NA_real_)
  }

  mtx <- mtx[keep_rows, , drop = FALSE]
  row_name <- row_name[keep_rows]
  mtx <- matrix(
    as.numeric(as.matrix(mtx)),
    nrow = nrow(mtx),
    ncol = ncol(mtx),
    dimnames = list(row_name, colnames(mtx))
  )
  ct_mtx <- rowsum(mtx, group = rownames(mtx))
  ct_mtx <- ct_mtx[rowSums(ct_mtx) > 0, , drop = FALSE]
  ct_mtx <- ct_mtx[, colSums(ct_mtx) > 0, drop = FALSE]
  if (nrow(ct_mtx) < 2 || ncol(ct_mtx) < 2) {
    return(NA_real_)
  }

  fisher_p <- tryCatch(
    fisher.test(ct_mtx, simulate.p.value = simulate_p_value, B = fisher_B)$p.value,
    error = function(e) NA_real_
  )
  return(fisher_p)
}

summarize_fisher_input_validity <- function(grouped_mtx, cell_type_map, source_label, match_type_label) {
  total_te <- as.integer(length(grouped_mtx))
  tested_te_ncol_gt1 <- 0L
  fisher_valid <- 0L
  fisher_skipped_ncol_le1 <- 0L
  fisher_skipped_missing_cell_type <- 0L
  fisher_skipped_invalid_contingency <- 0L

  if (length(grouped_mtx) > 0) {
    for (mtx in grouped_mtx) {
      if (is.null(mtx) || ncol(mtx) <= 1) {
        fisher_skipped_ncol_le1 <- fisher_skipped_ncol_le1 + 1L
        next
      }
      tested_te_ncol_gt1 <- tested_te_ncol_gt1 + 1L

      row_name <- as.vector(cell_type_map[rownames(mtx)])
      keep_rows <- !is.na(row_name)
      if (!any(keep_rows)) {
        fisher_skipped_missing_cell_type <- fisher_skipped_missing_cell_type + 1L
        next
      }

      mtx <- mtx[keep_rows, , drop = FALSE]
      row_name <- row_name[keep_rows]
      mtx <- matrix(
        as.numeric(as.matrix(mtx)),
        nrow = nrow(mtx),
        ncol = ncol(mtx),
        dimnames = list(row_name, colnames(mtx))
      )
      ct_mtx <- rowsum(mtx, group = rownames(mtx))
      ct_mtx <- ct_mtx[rowSums(ct_mtx) > 0, , drop = FALSE]
      ct_mtx <- ct_mtx[, colSums(ct_mtx) > 0, drop = FALSE]
      if (nrow(ct_mtx) < 2 || ncol(ct_mtx) < 2) {
        fisher_skipped_invalid_contingency <- fisher_skipped_invalid_contingency + 1L
      } else {
        fisher_valid <- fisher_valid + 1L
      }
    }
  }

  return(data.table(
    source = source_label,
    match_type = match_type_label,
    total_te = total_te,
    tested_te_ncol_gt1 = tested_te_ncol_gt1,
    fisher_valid = fisher_valid,
    fisher_skipped_ncol_le1 = fisher_skipped_ncol_le1,
    fisher_skipped_missing_cell_type = fisher_skipped_missing_cell_type,
    fisher_skipped_invalid_contingency = fisher_skipped_invalid_contingency
  ))
}

calculate_diff <- function(mtx, cell_type_map, index_type) {
  row_name <- cell_type_map[rownames(mtx)]
  mtx <- apply(mtx, 2, as.numeric)
  rownames(mtx) <- row_name
  ct_mtx <- rowsum(mtx, group = rownames(mtx))

  pos <- sapply(colnames(ct_mtx), function(col) {
    parts <- unlist(strsplit(col, split = "\\|"))
    if (length(parts) > 1) {
      pos_parts <- unlist(strsplit(parts[2], split = ":"))
      as.integer(pos_parts[2])
    } else {
      0
    }
  })

  strand <- sapply(colnames(ct_mtx), function(col) {
    parts <- unlist(strsplit(col, split = "\\|"))
    if (length(parts) > 1) {
      pos_parts <- unlist(strsplit(parts[2], split = ":"))
      pos_parts[length(pos_parts)]
    } else {
      "+"
    }
  })
  if (all(strand == "+")) {
    pos <- as.numeric(abs(pos - min(pos)))
  } else if (all(strand == "-")) {
    pos <- as.numeric(abs(pos - max(pos)))
  } else {
    stop("Mixed strands within the same data.table")
  }
  if (index_type == "PDUI") {
    index_weight <- as.numeric(pos == max(pos))
  } else if (index_type == "PPUI") {
    index_weight <- as.numeric(pos == min(pos))
  } else if (index_type == "RWUI") {
    index_weight <- as.numeric(order(pos)) - 1
    index_weight <- index_weight / max(index_weight)
  } else if (index_type == "DWUI") {
    index_weight <- pos / sum(pos)
  } else if (index_type == "MPRO") {
    ct_prop <- sweep(ct_mtx, 1, rowSums(ct_mtx), "/")
    mpro <- max(abs(ct_prop[1, ] - ct_prop[2, ]))
    return(mpro)
  } else {
    stop("Index type not found, please check")
  }
  index_weight_matrix <- matrix(index_weight, nrow = nrow(ct_mtx), ncol = ncol(ct_mtx), byrow = TRUE)
  weighted_sum <- rowSums(ct_mtx * index_weight_matrix)
  row_sums <- rowSums(ct_mtx)
  index <- weighted_sum / row_sums
  return(abs(index[1] - index[2]))
}

calculate_apa_detect_performance <- function(gt_filter_group, pd_filter_group, filter_type_1, filter_type_2) {
  gt <- intersect(gt_filter_group[[filter_type_1]], gt_filter_group[[filter_type_2]])
  pd <- intersect(pd_filter_group[[filter_type_1]], pd_filter_group[[filter_type_2]])
  overlap <- length(intersect(gt, pd))
  precision <- safe_divide(overlap, length(pd))
  recall <- safe_divide(overlap, length(gt))
  f1 <- safe_f1(precision, recall)
  result <- data.table(
    filter_type_1 = filter_type_1,
    filter_type_2 = filter_type_2,
    precision = precision,
    recall = recall,
    f1 = f1,
    gt_count = length(gt),
    pd_count = length(pd)
  )
  return(result)
}

compute_apa_statistics <- function(
    grouped_mtx,
    cell_type_map,
    index_types,
    diff_types,
    dexseq_split_num,
    adjustment_method,
    fisher_simulate_p_value,
    fisher_B,
    fisher_seed = NULL,
    catch_fisher_errors = FALSE
) {
  p_list <- list()
  suppressWarnings({
    for (index in index_types) {
      index_values <- lapply(grouped_mtx, calculate_index, index_type = index)
      p_list[[paste0("wilcox_", index)]] <- lapply(index_values, wilcox_test_by_cell_type, cell_type_map = cell_type_map) %>% p.adjust(method = adjustment_method)
    }
  })

  if (fisher_simulate_p_value && !is.null(fisher_seed)) {
    set.seed(fisher_seed)
  }

  p_list[["fisher"]] <- lapply(grouped_mtx, function(mtx, cell_type_map) {
    if (catch_fisher_errors) {
      return(tryCatch({
        fisher_test_by_cell_type(
          mtx = mtx,
          cell_type_map = cell_type_map,
          simulate_p_value = fisher_simulate_p_value,
          fisher_B = fisher_B
        )
      }, error = function(e) {
        return(NA)
      }))
    }
    fisher_test_by_cell_type(
      mtx = mtx,
      cell_type_map = cell_type_map,
      simulate_p_value = fisher_simulate_p_value,
      fisher_B = fisher_B
    )
  }, cell_type_map = cell_type_map) %>% p.adjust(method = adjustment_method)

  dxr <- dexseq_test_by_cell_type(grouped_mtx, cell_type_map, split_num = dexseq_split_num)

  diff_list <- list()
  for (diff_type in diff_types) {
    diff_list[[diff_type]] <- lapply(grouped_mtx, calculate_diff, cell_type_map = cell_type_map, index_type = diff_type)
  }

  return(list(
    p_list = p_list,
    dxr = dxr,
    diff_list = diff_list
  ))
}

build_filter_group <- function(
    p_list,
    dxr,
    diff_list,
    p_thresholds,
    log2fc_thresholds,
    index_diff_thresholds
) {
  filter_group <- list()

  for (p_threshold in p_thresholds) {
    for (stats in names(p_list)) {
      p_vec <- p_list[[stats]]
      filter_group[[paste0(stats, "_", p_threshold)]] <- p_vec[!is.na(p_vec) & p_vec < p_threshold] %>% names
    }
    filter_group[[paste0("dexseq", "_", p_threshold)]] <- dxr[!is.na(padj) & padj < p_threshold, unique(groupID)] %>% as.vector %>% gsub(., pattern = "_", replacement = ":")
  }

  for (log2fc_threshold in log2fc_thresholds) {
    filter_group[[paste0("dexseq_log2fc", "_", log2fc_threshold)]] <- dxr[!is.na(log2fold_T2_T1) & abs(log2fold_T2_T1) > log2fc_threshold, unique(groupID)] %>% as.vector %>% gsub(., pattern = "_", replacement = ":")
  }

  for (index_diff_threshold in index_diff_thresholds) {
    for (diff_type in names(diff_list)) {
      diff_vec <- diff_list[[diff_type]]
      filter_group[[paste0(diff_type, "_", index_diff_threshold)]] <- diff_vec[!is.na(diff_vec) & diff_vec > index_diff_threshold] %>% names
    }
  }

  return(filter_group)
}

precompute_gt_apa_statistics <- function(
    grouped_gt_mtx,
    cell_type_map,
    index_types = c("PDUI", "PPUI", "RWUI", "DWUI"),
    diff_types = c("PDUI", "PPUI", "RWUI", "DWUI", "MPRO"),
    p_thresholds = c(0.05, 0.01),
    log2fc_thresholds = c(0.25, 0.5, 0.75, 1, 1.25),
    index_diff_thresholds = c(0.1, 0.2, 0.3, 0.4, 0.5),
    dexseq_split_num = 6,
    adjustment_method = "fdr",
    fisher_simulate_p_value = TRUE,
    fisher_B = 2000,
    fisher_seed = 1234
) {
  gt_stats <- compute_apa_statistics(
    grouped_mtx = grouped_gt_mtx,
    cell_type_map = cell_type_map,
    index_types = index_types,
    diff_types = diff_types,
    dexseq_split_num = dexseq_split_num,
    adjustment_method = adjustment_method,
    fisher_simulate_p_value = fisher_simulate_p_value,
    fisher_B = fisher_B,
    fisher_seed = fisher_seed,
    catch_fisher_errors = FALSE
  )

  gt_filter_group <- build_filter_group(
    p_list = gt_stats$p_list,
    dxr = gt_stats$dxr,
    diff_list = gt_stats$diff_list,
    p_thresholds = p_thresholds,
    log2fc_thresholds = log2fc_thresholds,
    index_diff_thresholds = index_diff_thresholds
  )

  return(list(
    stats = gt_stats,
    filter_group = gt_filter_group
  ))
}

perform_full_apa_analysis <- function(
    grouped_gt_mtx,
    grouped_pd_mtx,
    cell_type_map,
    sample,
    tool,
    protocol,
    index_types = c("PDUI", "PPUI", "RWUI", "DWUI"),
    diff_types = c("PDUI", "PPUI", "RWUI", "DWUI", "MPRO"),
    p_thresholds = c(0.05, 0.01),
    log2fc_thresholds = c(0.25, 0.5, 0.75, 1, 1.25),
    index_diff_thresholds = c(0.1, 0.2, 0.3, 0.4, 0.5),
    dexseq_split_num = 6,
    adjustment_method = "fdr",
    fisher_simulate_p_value = TRUE,
    fisher_B = 2000,
    fisher_seed = 1234,
    precomputed_gt = NULL
) {
  apa_detect_performance <- data.table()

  stopifnot(
    is.list(grouped_gt_mtx),
    is.list(grouped_pd_mtx),
    !is.null(names(grouped_gt_mtx)),
    !is.null(names(grouped_pd_mtx))
  )

  pd_stats <- compute_apa_statistics(
    grouped_mtx = grouped_pd_mtx,
    cell_type_map = cell_type_map,
    index_types = index_types,
    diff_types = diff_types,
    dexseq_split_num = dexseq_split_num,
    adjustment_method = adjustment_method,
    fisher_simulate_p_value = fisher_simulate_p_value,
    fisher_B = fisher_B,
    fisher_seed = fisher_seed,
    catch_fisher_errors = TRUE
  )

  pd_filter_group <- build_filter_group(
    p_list = pd_stats$p_list,
    dxr = pd_stats$dxr,
    diff_list = pd_stats$diff_list,
    p_thresholds = p_thresholds,
    log2fc_thresholds = log2fc_thresholds,
    index_diff_thresholds = index_diff_thresholds
  )

  if (!is.null(precomputed_gt)) {
    gt_filter_group <- precomputed_gt$filter_group
  } else {
    gt_precomputed <- precompute_gt_apa_statistics(
      grouped_gt_mtx = grouped_gt_mtx,
      cell_type_map = cell_type_map,
      index_types = index_types,
      diff_types = diff_types,
      p_thresholds = p_thresholds,
      log2fc_thresholds = log2fc_thresholds,
      index_diff_thresholds = index_diff_thresholds,
      dexseq_split_num = dexseq_split_num,
      adjustment_method = adjustment_method,
      fisher_simulate_p_value = fisher_simulate_p_value,
      fisher_B = fisher_B,
      fisher_seed = fisher_seed
    )
    gt_filter_group <- gt_precomputed$filter_group
  }

  p_stats <- c(paste0("wilcox_", index_types), "fisher", "dexseq")
  p_stats <- expand.grid(p_stats, p_thresholds, stringsAsFactors = FALSE) %>% apply(1, paste0, collapse = "_")
  diff_stats <- diff_types
  diff_stats <- expand.grid(diff_stats, index_diff_thresholds, stringsAsFactors = FALSE) %>% apply(1, paste0, collapse = "_")
  fc_stats <- c("dexseq_log2fc")
  fc_stats <- expand.grid(fc_stats, as.character(log2fc_thresholds), stringsAsFactors = FALSE) %>% apply(1, paste0, collapse = "_")
  diff_stats <- c(diff_stats, fc_stats)

  for (filter_type_1 in p_stats) {
    for (filter_type_2 in diff_stats) {
      result <- calculate_apa_detect_performance(gt_filter_group, pd_filter_group, filter_type_1, filter_type_2)
      apa_detect_performance <- rbind(apa_detect_performance, result)
    }
  }
  apa_detect_performance[, sample := sample]
  apa_detect_performance[, tool := tool]
  apa_detect_performance[, protocol := protocol]

  return(apa_detect_performance)
}

global_analysis_config <- function() {
  return(list(
    window_sizes = c(50, 100, 150, 200),
    index_types = c("PDUI", "PPUI", "RWUI", "DWUI"),
    diff_types = c("PDUI", "PPUI", "RWUI", "DWUI", "MPRO"),
    p_thresholds = c(0.05, 0.01),
    log2fc_thresholds = c(0.25, 0.5, 0.75, 1, 1.25),
    index_diff_thresholds = c(0.1, 0.2, 0.3, 0.4, 0.5),
    dexseq_split_num = 6,
    adjustment_method = "fdr",
    fisher_simulate_p_value = TRUE,
    fisher_B = 2000,
    fisher_seed = 1234
  ))
}
