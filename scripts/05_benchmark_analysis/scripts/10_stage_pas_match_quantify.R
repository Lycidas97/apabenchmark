library(optparse)
library(Seurat)
library(DEXSeq)
library(data.table)
library(magrittr)
library(dplyr)

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

script_arg <- commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]
script_path <- sub("--file=", "", script_arg[1])
script_dir <- dirname(normalizePath(script_path))
source(file.path(script_dir, "common_performance_functions.R"))

option_list <- list(
  make_option(c("--tool"), dest = "tool", type = "character"),
  make_option(c("--sample"), dest = "sample", type = "character"),
  make_option(c("--gt_pas"), dest = "gt_pas", type = "character"),
  make_option(c("--pd_pas"), dest = "pd_pas", type = "character"),
  make_option(c("--gt_mtx"), dest = "gt_mtx", type = "character"),
  make_option(c("--pd_mtx"), dest = "pd_mtx", type = "character"),
  make_option(c("--output_dir"), dest = "output_dir", type = "character"),
  make_option(c("--dexseq_split_num"), dest = "dexseq_split_num", type = "integer", default = NA_integer_)
)

opt <- parse_args(OptionParser(option_list = option_list))
required_args <- c("tool", "sample", "gt_pas", "pd_pas", "gt_mtx", "pd_mtx", "output_dir", "dexseq_split_num")
for (arg in required_args) {
  if (is.null(opt[[arg]]) || opt[[arg]] == "" || (arg == "dexseq_split_num" && is.na(opt[[arg]]))) {
    stop(sprintf("Missing required argument: --%s", arg))
  }
}

tool <- opt$tool
sample <- opt$sample
gt_pas_path <- opt$gt_pas
pd_pas_path <- opt$pd_pas
gt_mtx_path <- opt$gt_mtx
pd_mtx_path <- opt$pd_mtx
output_dir <- opt$output_dir

for (input_path in c(gt_pas_path, pd_pas_path, gt_mtx_path, pd_mtx_path)) {
  if (!file.exists(input_path)) {
    stop(sprintf("Input file not found: %s", input_path))
  }
}

bedtools_bin <- resolve_bedtools_bin(bedtools_path)

protocol <- unlist(strsplit(sample, "_"))[1]
cfg <- global_analysis_config()
if (opt$dexseq_split_num < 1) {
  stop("--dexseq_split_num must be >= 1")
}
cfg$dexseq_split_num <- as.integer(opt$dexseq_split_num)

quality_issues <- data.table(
  sample = character(),
  tool = character(),
  protocol = character(),
  stage = character(),
  issue_code = character(),
  severity = character(),
  detail = character(),
  metric_name = character(),
  metric_value = character()
)
append_quality_issue <- function(stage, issue_code, severity = "warning", detail = "", metric_name = NA_character_, metric_value = NA_character_) {
  if (is.null(metric_name) || length(metric_name) == 0) {
    metric_name <- NA_character_
  } else {
    metric_name <- as.character(metric_name)
  }
  if (is.null(metric_value) || length(metric_value) == 0) {
    metric_value <- NA_character_
  } else {
    metric_value <- as.character(metric_value)
  }
  quality_issues <<- rbind(
    quality_issues,
    data.table(
      sample = sample,
      tool = tool,
      protocol = protocol,
      stage = stage,
      issue_code = issue_code,
      severity = severity,
      detail = as.character(detail),
      metric_name = metric_name,
      metric_value = metric_value
    ),
    use.names = TRUE,
    fill = TRUE
  )
}

gt_mtx <- read_count_matrix_dt(gt_mtx_path, first_col_name = "barcode") %>% as.matrix
rownames(gt_mtx) <- gt_mtx[, 1]
gt_mtx <- gt_mtx[, -1]

gt_pas_dt <- fread(gt_pas_path)
setnames(gt_pas_dt, c("exon_id", "pas_type"), c("te_name", "te_type"))
gt_pas_dt[, pas_name := paste0(te_name, "|", chr, ":", start, ":", end, ":", strand)]
gt_pas_dt[, pas_num := .N, by = te_name]
gt_col_num <- ncol(gt_pas_dt)

gt_te_dt <- gt_pas_dt[, c("te_name", "pas_num"), with = FALSE]
gt_te_dt <- gt_te_dt[!duplicated(te_name)]
gt_te_dt[, chr := sapply(strsplit(te_name, ":"), `[`, 3), by = .(te_name)]
gt_te_dt[, start := as.integer(sapply(strsplit(te_name, ":"), `[`, 4)), by = .(te_name)]
gt_te_dt[, end := as.integer(sapply(strsplit(te_name, ":"), `[`, 5)), by = .(te_name)]
gt_te_dt[, strand := sapply(strsplit(te_name, ":"), `[`, 6), by = .(te_name)]
gt_te_dt <- gt_te_dt[order(chr, start, end, strand)]
gt_te_dt <- gt_te_dt[, c("chr", "start", "end", "te_name", "pas_num", "strand"), with = FALSE]

if (!("cell_type" %in% colnames(gt_mtx))) {
  stop(sprintf("cell_type column not found in GT matrix: %s", gt_mtx_path))
}

gt_pas_cols <- colnames(gt_mtx)[grepl("\\|", colnames(gt_mtx))]
if (length(gt_pas_cols) == 0) {
  stop(sprintf("No GT PAS columns detected in matrix: %s", gt_mtx_path))
}

gt_te_pas_map <- split(gt_pas_cols, sapply(gt_pas_cols, function(col) {
  unlist(strsplit(col, split = "\\|"))[1]
}))
cell_type_map <- as.list(gt_mtx[, "cell_type"])
grouped_gt_mtx <- lapply(gt_te_pas_map, function(cols) {
  gt_mtx[, cols, drop = FALSE]
})

pd_mtx <- matrix(
  numeric(0),
  nrow = nrow(gt_mtx),
  ncol = 0,
  dimnames = list(rownames(gt_mtx), character(0))
)

pd_mtx_raw <- NULL
pd_mtx_error <- NULL
pd_file_size <- file.info(pd_mtx_path)$size
if (is.na(pd_file_size) || pd_file_size == 0) {
  append_quality_issue(
    stage = "stage10",
    issue_code = "pd_matrix_file_empty",
    severity = "warning",
    detail = "PD matrix file is empty; continuing with empty PD matrix."
  )
  warning(
    sprintf("PD matrix file is empty for %s/%s; continuing with empty PD matrix.", tool, sample),
    call. = FALSE
  )
} else {
  pd_mtx_raw <- tryCatch(
    read_count_matrix_dt(pd_mtx_path, first_col_name = "barcode") %>% as.matrix,
    error = function(e) {
      pd_mtx_error <<- conditionMessage(e)
      NULL
    }
  )
}

if (is.null(pd_mtx_raw)) {
  if (!is.null(pd_mtx_error)) {
    append_quality_issue(
      stage = "stage10",
      issue_code = "pd_matrix_parse_failed",
      severity = "warning",
      detail = sprintf("Failed to parse PD matrix: %s", pd_mtx_error)
    )
    warning(
      sprintf(
        "Failed to parse PD matrix for %s/%s (%s); continuing with empty PD matrix.",
        tool, sample, pd_mtx_error
      ),
      call. = FALSE
    )
  }
} else if (ncol(pd_mtx_raw) < 1) {
  append_quality_issue(
    stage = "stage10",
    issue_code = "pd_matrix_no_columns",
    severity = "warning",
    detail = "PD matrix has no columns; continuing with empty PD matrix."
  )
  warning(
    sprintf("PD matrix has no columns for %s/%s; continuing with empty PD matrix.", tool, sample),
    call. = FALSE
  )
} else {
  rownames(pd_mtx_raw) <- as.character(pd_mtx_raw[, 1])
  if (ncol(pd_mtx_raw) <= 1) {
    append_quality_issue(
      stage = "stage10",
      issue_code = "pd_matrix_no_pas_columns",
      severity = "warning",
      detail = "PD matrix has no PAS columns; continuing with empty PD matrix."
    )
    warning(
      sprintf("PD matrix has no PAS columns for %s/%s; continuing with empty PD matrix.", tool, sample),
      call. = FALSE
    )
  } else {
    pd_feature_mtx <- pd_mtx_raw[, -1, drop = FALSE]
    pd_feature_mtx <- matrix(
      as.numeric(pd_feature_mtx),
      nrow = nrow(pd_feature_mtx),
      ncol = ncol(pd_feature_mtx),
      dimnames = dimnames(pd_feature_mtx)
    )

    matched_idx <- match(rownames(gt_mtx), rownames(pd_feature_mtx))
    pd_mtx <- matrix(
      0,
      nrow = nrow(gt_mtx),
      ncol = ncol(pd_feature_mtx),
      dimnames = list(rownames(gt_mtx), colnames(pd_feature_mtx))
    )
    valid_idx <- !is.na(matched_idx)
    if (any(valid_idx)) {
      pd_mtx[valid_idx, ] <- pd_feature_mtx[matched_idx[valid_idx], , drop = FALSE]
    } else {
      append_quality_issue(
        stage = "stage10",
        issue_code = "pd_matrix_no_barcode_overlap",
        severity = "warning",
        detail = "No overlapping barcodes between GT and PD matrices; continuing with empty PD matrix."
      )
      warning(
        sprintf("No overlapping barcodes between GT and PD matrices for %s/%s; continuing with empty PD matrix.", tool, sample),
        call. = FALSE
      )
    }

    pd_col_before_nonzero_filter <- ncol(pd_mtx)
    pd_mtx <- pd_mtx[, colSums(pd_mtx != 0) > 0, drop = FALSE]
    pd_col_after_nonzero_filter <- ncol(pd_mtx)
    dropped_zero_cols <- pd_col_before_nonzero_filter - pd_col_after_nonzero_filter
    if (dropped_zero_cols > 0) {
      append_quality_issue(
        stage = "stage10",
        issue_code = "pd_matrix_zero_columns_dropped",
        severity = "info",
        detail = sprintf("Dropped %d all-zero PD PAS columns after GT barcode alignment.", dropped_zero_cols),
        metric_name = "dropped_zero_pd_columns",
        metric_value = dropped_zero_cols
      )
    }
    if (pd_col_after_nonzero_filter == 0) {
      append_quality_issue(
        stage = "stage10",
        issue_code = "pd_matrix_all_columns_zero_after_alignment",
        severity = "warning",
        detail = "All PD PAS columns are zero after GT barcode alignment."
      )
    }
  }
}

pd_pas_raw <- read_table_if_exists(pd_pas_path)
if (nrow(pd_pas_raw) == 0 || ncol(pd_pas_raw) == 0) {
  append_quality_issue(
    stage = "stage10",
    issue_code = "pd_pas_bed_empty",
    severity = "warning",
    detail = "PD PAS BED is empty; downstream matching will be empty."
  )
  pd_pas_dt <- data.table(
    chr = character(),
    start = integer(),
    end = integer(),
    name = character(),
    pos = character(),
    strand = character()
  )
} else {
  if (ncol(pd_pas_raw) < 6) {
    stop(sprintf("PD PAS BED requires at least 6 columns, got %d: %s", ncol(pd_pas_raw), pd_pas_path))
  }
  pd_pas_dt <- pd_pas_raw[, 1:6, with = FALSE]
  setnames(pd_pas_dt, names(pd_pas_dt), c("chr", "start", "end", "name", "pos", "strand"))
}

pd_pas_dt <- pd_pas_dt[order(chr, start, end)]
if (nrow(pd_pas_dt) > 0) {
  pd_pas_dt[, rank := {
    strand_first <- strand[1]
    if (!is.na(strand_first) && strand_first == "-") {
      (seq_along(strand)) - 1
    } else {
      (.N - seq_along(strand))
    }
  }, by = name]
} else {
  pd_pas_dt[, rank := integer()]
}
pd_pas_dt <- pd_pas_dt[rank == 0]
pd_pas_dt[, pas_name := paste0("|", chr, ":", start, ":", end, ":", strand)]
pd_col_num <- ncol(pd_pas_dt)
pd_pas_rows_before_matrix_filter <- nrow(pd_pas_dt)
pd_pas_dt <- pd_pas_dt[name %in% colnames(pd_mtx), ]
pd_pas_rows_after_matrix_filter <- nrow(pd_pas_dt)
if (pd_pas_rows_before_matrix_filter > pd_pas_rows_after_matrix_filter) {
  append_quality_issue(
    stage = "stage10",
    issue_code = "pd_pas_bed_rows_not_in_matrix",
    severity = "info",
    detail = sprintf(
      "Dropped %d PD PAS BED rows not present in PD matrix columns.",
      pd_pas_rows_before_matrix_filter - pd_pas_rows_after_matrix_filter
    ),
    metric_name = "dropped_pd_pas_bed_rows",
    metric_value = pd_pas_rows_before_matrix_filter - pd_pas_rows_after_matrix_filter
  )
}
if (ncol(pd_mtx) > 0) {
  if (nrow(pd_pas_dt) == 0) {
    append_quality_issue(
      stage = "stage10",
      issue_code = "pd_matrix_no_pas_coordinate_mapping",
      severity = "warning",
      detail = "No PAS coordinate mapping for PD matrix columns; dropping all PD columns."
    )
    warning(
      sprintf("No PAS coordinate mapping for PD matrix columns in %s/%s; dropping all PD columns.", tool, sample),
      call. = FALSE
    )
    pd_mtx <- pd_mtx[, 0, drop = FALSE]
  } else {
    matched_pas_idx <- match(colnames(pd_mtx), pd_pas_dt$name)
    keep_cols <- !is.na(matched_pas_idx)
    if (!all(keep_cols)) {
      append_quality_issue(
        stage = "stage10",
        issue_code = "pd_matrix_columns_without_pas_mapping",
        severity = "warning",
        detail = sprintf("Dropping %d PD columns without PAS mapping.", sum(!keep_cols)),
        metric_name = "dropped_unmapped_pd_columns",
        metric_value = sum(!keep_cols)
      )
      warning(
        sprintf("Dropping %d PD columns without PAS mapping for %s/%s.", sum(!keep_cols), tool, sample),
        call. = FALSE
      )
      pd_mtx <- pd_mtx[, keep_cols, drop = FALSE]
      matched_pas_idx <- matched_pas_idx[keep_cols]
    }
    if (ncol(pd_mtx) > 0) {
      colnames(pd_mtx) <- pd_pas_dt$pas_name[matched_pas_idx]
    }
  }
}

match_performance <- data.table()
pas_quantify_performance <- data.table()
apa_inputs <- list()

for (window_size in cfg$window_sizes) {
  gt_pas_bed_tmp <- tempfile(fileext = ".bed")
  predict_pas_tmp <- tempfile(fileext = ".bed")
  fwrite(gt_pas_dt, gt_pas_bed_tmp, sep = "\t", quote = FALSE, col.names = FALSE)
  fwrite(pd_pas_dt, predict_pas_tmp, sep = "\t", quote = FALSE, col.names = FALSE)
  output_file <- tempfile()
  window_stderr <- tempfile()
  window_status <- system2(
    bedtools_bin,
    args = c("window", "-sw", "-sm", "-a", gt_pas_bed_tmp, "-b", predict_pas_tmp, "-w", as.character(window_size / 2)),
    stdout = output_file,
    stderr = window_stderr
  )
  if (!identical(window_status, 0L)) {
    err_msg <- if (file.exists(window_stderr)) paste(readLines(window_stderr, warn = FALSE), collapse = "\n") else ""
    stop(sprintf("bedtools window failed for window_size=%d: %s", window_size, err_msg))
  }
  unlink(c(window_stderr))

  match_dt <- read_table_if_exists(output_file)
  unlink(c(gt_pas_bed_tmp, predict_pas_tmp, output_file))
  mt <- paste0("pas_", window_size / 2)
  pd_mt <- paste0("pd_", window_size / 2)
  gt_mt <- paste0("gt_", window_size / 2)

  if (nrow(match_dt) == 0) {
    append_quality_issue(
      stage = "stage10",
      issue_code = "pas_window_match_empty",
      severity = "warning",
      detail = sprintf("No PAS window matches at %s.", mt),
      metric_name = "match_type",
      metric_value = mt
    )
    match_performance <- rbind(match_performance, data.table(sample = sample, tool = tool, protocol = protocol, match_type = mt, tp = 0, fp = 0, fn = nrow(gt_pas_dt), precision = 0, recall = 0, f1 = 0))
    pas_quantify_performance <- rbind(pas_quantify_performance, data.table(sample = sample, tool = tool, protocol = protocol, match_type = mt, cor_pas = 0, rmse_pas = 0, mae_pas = 0, mape_pas = 0, rmse_pas_ct = 0, mae_pas_ct = 0, mape_pas_ct = 0))
    apa_inputs[[pd_mt]] <- list()
    apa_inputs[[gt_mt]] <- list()
    next
  }

  rename_matched_columns(match_dt, gt_col_num, pd_col_num)
  match_dt[, distance := abs(gt_start - res_start)]
  match_dt[, best_match := distance == min(distance), by = res_name]
  match_dt <- match_dt[best_match == TRUE]
  match_dt[, best_match := distance == min(distance), by = gt_name]

  pas_merged_dt <- merge(gt_pas_dt, match_dt,
    by.x = c("chr", "start", "end", "strand"),
    by.y = c("gt_chr", "gt_start", "gt_end", "gt_strand"),
    all.x = TRUE
  )
  gt_cols <- names(pas_merged_dt)[grepl("^gt_", names(pas_merged_dt))]
  pas_merged_dt[, (gt_cols) := NULL]
  pas_merged_dt[, all_match := !is.na(best_match)]
  pas_merged_dt[is.na(best_match), best_match := FALSE]
  setnames(pas_merged_dt, paste0("res_", pd_col_num), "pd_pas_name")

  pas_match_tp <- pas_merged_dt[all_match == TRUE, .N]
  pas_match_fp <- pd_pas_dt[, .N] - pas_match_tp
  pas_match_fn <- gt_pas_dt[, .N] - pas_match_tp
  pas_match_precision <- safe_divide(pas_match_tp, pas_match_tp + pas_match_fp)
  pas_match_recall <- safe_divide(pas_match_tp, pas_match_tp + pas_match_fn)
  pas_match_f1 <- safe_f1(pas_match_precision, pas_match_recall)
  match_performance <- rbind(match_performance, data.table(sample = sample, tool = tool, protocol = protocol, match_type = mt, tp = pas_match_tp, fp = pas_match_fp, fn = pas_match_fn, precision = pas_match_precision, recall = pas_match_recall, f1 = pas_match_f1))

  current_pas_map <- pas_merged_dt[all_match == TRUE, .(te_name, gt_pas_name = pas_name, pd_pas_name)]
  current_pas_map[, pd_pas_count := .N, by = te_name]
  current_pas_map[, gt_pas_count := uniqueN(gt_pas_name), by = te_name]

  multipd_pd_te_pas_map <- split(current_pas_map[(pd_pas_count > 1), pd_pas_name], current_pas_map[(pd_pas_count > 1), te_name])
  grouped_pd_mtx_mapped_by_pd <- lapply(multipd_pd_te_pas_map, function(cols) pd_mtx[, cols, drop = FALSE])
  if (length(grouped_pd_mtx_mapped_by_pd) > 0) {
    keep_pd <- vapply(grouped_pd_mtx_mapped_by_pd, function(x) !is.null(x) && ncol(x) > 0, logical(1))
    grouped_pd_mtx_mapped_by_pd <- grouped_pd_mtx_mapped_by_pd[keep_pd]
  }
  apa_inputs[[pd_mt]] <- grouped_pd_mtx_mapped_by_pd

  gt_pas_mapping <- current_pas_map[gt_pas_count > 1, .(pd_pas_list = list(unique(pd_pas_name))), by = .(te_name, gt_pas_name)]
  grouped_pd_mtx_mapped_by_gt <- lapply(split(gt_pas_mapping, by = "te_name"), function(te_data) {
    all_pd_cols <- unique(unlist(te_data$pd_pas_list))
    valid_cols <- intersect(all_pd_cols, colnames(pd_mtx))
    sub_mtx <- as.matrix(pd_mtx[, valid_cols, drop = FALSE])
    storage.mode(sub_mtx) <- "numeric"
    result_mtx <- matrix(nrow = nrow(sub_mtx), ncol = nrow(te_data), dimnames = list(rownames(sub_mtx), te_data$gt_pas_name))
    for (i in seq_len(nrow(te_data))) {
      target_cols <- intersect(te_data$pd_pas_list[[i]], valid_cols)
      if (length(target_cols) > 0) {
        result_mtx[, i] <- rowSums(sub_mtx[, target_cols, drop = FALSE], na.rm = TRUE)
      } else {
        result_mtx[, i] <- NA_real_
      }
    }
    result_mtx[, colSums(is.na(result_mtx)) < nrow(result_mtx), drop = FALSE]
  })
  if (length(grouped_pd_mtx_mapped_by_gt) > 0) {
    keep_gt <- vapply(grouped_pd_mtx_mapped_by_gt, function(x) !is.null(x) && ncol(x) > 0, logical(1))
    grouped_pd_mtx_mapped_by_gt <- grouped_pd_mtx_mapped_by_gt[keep_gt]
  }
  apa_inputs[[gt_mt]] <- grouped_pd_mtx_mapped_by_gt

  mapped_pd_mtx <- pd_mtx %>% apply(., 2, as.numeric)
  colnames(mapped_pd_mtx) <- pas_merged_dt[match(colnames(pd_mtx), pas_merged_dt$pd_pas_name), pas_name]
  rownames(mapped_pd_mtx) <- rownames(pd_mtx)
  mapped_pd_mtx <- mapped_pd_mtx[, !is.na(colnames(mapped_pd_mtx)), drop = FALSE]
  mapped_pd_mtx <- mapped_pd_mtx %>% t %>% rowsum(., group = colnames(mapped_pd_mtx)) %>% t

  mapped_gt_mtx <- gt_mtx[, -1, drop = FALSE] %>% apply(., 2, as.numeric)
  colnames(mapped_gt_mtx) <- colnames(gt_mtx[, -1, drop = FALSE])
  rownames(mapped_gt_mtx) <- rownames(gt_mtx)
  shared_pas <- intersect(colnames(mapped_pd_mtx), colnames(mapped_gt_mtx))
  if (length(shared_pas) == 0) {
    append_quality_issue(
      stage = "stage10",
      issue_code = "pas_quantify_no_shared_pas",
      severity = "warning",
      detail = sprintf("No shared PAS between mapped GT and PD matrices at %s.", mt),
      metric_name = "match_type",
      metric_value = mt
    )
    pas_quantify_performance <- rbind(pas_quantify_performance, data.table(sample = sample, tool = tool, protocol = protocol, match_type = mt, cor_pas = 0, rmse_pas = 0, mae_pas = 0, mape_pas = 0, rmse_pas_ct = 0, mae_pas_ct = 0, mape_pas_ct = 0))
    next
  }
  mapped_pd_mtx <- mapped_pd_mtx[, shared_pas, drop = FALSE]
  mapped_gt_mtx <- mapped_gt_mtx[rownames(mapped_pd_mtx), shared_pas, drop = FALSE]

  mapped_ct_gt_mtx <- mapped_gt_mtx
  mapped_ct_pd_mtx <- mapped_pd_mtx
  rownames(mapped_ct_gt_mtx) <- cell_type_map[rownames(mapped_gt_mtx)]
  rownames(mapped_ct_pd_mtx) <- cell_type_map[rownames(mapped_pd_mtx)]
  mapped_ct_gt_mtx <- rowsum(mapped_ct_gt_mtx, group = rownames(mapped_ct_gt_mtx))
  mapped_ct_pd_mtx <- rowsum(mapped_ct_pd_mtx, group = rownames(mapped_ct_pd_mtx))

  cor_pas <- mapply(cor, as.data.frame(mapped_pd_mtx), as.data.frame(mapped_gt_mtx), MoreArgs = list(method = "pearson")) %>% mean
  rmse_pas <- sqrt(mean((mapped_gt_mtx - mapped_pd_mtx)^2))
  mae_pas <- mean(abs(mapped_gt_mtx - mapped_pd_mtx))
  mape_pas <- mean(ifelse(mapped_gt_mtx == 0, 0, abs(mapped_gt_mtx - mapped_pd_mtx) / mapped_gt_mtx))
  rmse_pas_ct <- sqrt(mean((mapped_ct_gt_mtx - mapped_ct_pd_mtx)^2))
  mae_pas_ct <- mean(abs(mapped_ct_gt_mtx - mapped_ct_pd_mtx))
  mape_pas_ct <- mean(ifelse(mapped_ct_gt_mtx == 0, 0, abs(mapped_ct_gt_mtx - mapped_ct_pd_mtx) / mapped_ct_gt_mtx))
  pas_quantify_performance <- rbind(pas_quantify_performance, data.table(sample = sample, tool = tool, protocol = protocol, match_type = mt, cor_pas = cor_pas, rmse_pas = rmse_pas, mae_pas = mae_pas, mape_pas = mape_pas, rmse_pas_ct = rmse_pas_ct, mae_pas_ct = mae_pas_ct, mape_pas_ct = mape_pas_ct))
}

stage_dir <- file.path(output_dir, tool)
if (!dir.exists(stage_dir)) {
  dir.create(stage_dir, recursive = TRUE)
}
prefix <- file.path(stage_dir, sample)
bundle_path <- paste0(prefix, "_stage10_bundle.rds")

stage10_bundle <- list(
  sample = sample,
  tool = tool,
  protocol = protocol,
  output_dir = output_dir,
  bedtools_bin = bedtools_bin,
  cfg = cfg,
  gt_te_dt = gt_te_dt,
  pd_pas_dt = pd_pas_dt,
  pd_mtx = pd_mtx,
  pd_col_num = pd_col_num,
  grouped_gt_mtx = grouped_gt_mtx,
  cell_type_map = cell_type_map,
  match_performance = match_performance,
  pas_quantify_performance = pas_quantify_performance,
  apa_inputs = apa_inputs,
  quality_issues = quality_issues
)

saveRDS(stage10_bundle, bundle_path)
fwrite(match_performance, paste0(prefix, "_stage10_match_performance.tsv"), sep = "\t", quote = FALSE)
fwrite(pas_quantify_performance, paste0(prefix, "_stage10_pas_quantify_performance.tsv"), sep = "\t", quote = FALSE)
fwrite(quality_issues, paste0(prefix, "_stage10_quality_issues.tsv"), sep = "\t", quote = FALSE)

cat(sprintf("stage10 bundle saved: %s\n", bundle_path))
