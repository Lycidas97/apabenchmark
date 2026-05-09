#' @title Analyze APA Detection Performance
#' @description This script performs a comprehensive analysis of APA (Alternative Polyadenylation) detection performance by comparing ground truth (GT) and predicted (PD) PAS (Polyadenylation Site) data. It includes matching PAS positions, quantifying PAS usage, and evaluating differential APA events across cell types.
#' @details The script reads GT and PD PAS data, matches PAS positions using bedtools, and calculates various performance metrics such as precision, recall, and F1 score. It also evaluates the correlation and error between GT and PD PAS quantification. Additionally, it performs statistical tests (Wilcoxon, Fisher, DEXSeq) to detect differential APA events and calculates performance metrics for APA detection.
#' @param gt_pas_path Path to the ground truth PAS file in BED format.
#' @param pd_pas_path Path to the predicted PAS file in BED format.
#' @param gt_mtx_path Path to the ground truth PAS quantification matrix.
#' @param pd_mtx_path Path to the predicted PAS quantification matrix.
#' @param output_dir Directory to save the output performance metrics.
#' @param tool Name of the tool used for PAS prediction.
#' @param sample Name of the sample being analyzed.
#' @param protocol Protocol used for data generation (e.g., SpatialTranscriptomics).
#' @return Writes performance metrics to output files in the specified directory.
#' @examples
#' # Run the script with specified paths and parameters
#' # Rscript calculate_benchmark_performance.R --tool scapatrap --sample SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2 --gt_pas /path/to/apabenchmark_final/data/sim_pas/mm10_sim_pas_gn5000_rep2.bed --pd_pas /path/to/apabenchmark_final/data/sim_bam_result/scapatrap/SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2/pas.bed --gt_mtx /path/to/apabenchmark_final/data/sim_bam/SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2.bam.expr.tsv --pd_mtx /path/to/apabenchmark_final/data/sim_bam_result/scapatrap/SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2/pas_counts.tsv --output_dir /path/to/apabenchmark_final/data/performance

library(optparse)
library(Seurat)
library(DEXSeq)
library(data.table)
library(magrittr)
library(dplyr)

bedtools_bin <- Sys.getenv("APABENCHMARK_BEDTOOLS_BIN", "bedtools")
if (dir.exists(bedtools_bin)) {
  bedtools_bin <- file.path(bedtools_bin, "bedtools")
}

# Define the functions
rename_matched_columns <- function(match_dt, gt_col_num, res_col_num) {
  column_names <- c(
    'gt_chr', 'gt_start', 'gt_end', 'gt_name', 'gt_score', 'gt_strand',
    'res_chr', 'res_start', 'res_end', 'res_name', 'res_score', 'res_strand'
  )
  
  if (gt_col_num > 6) {
    for (i in 6:(gt_col_num - 1)) {
      column_names <- c(column_names[1:i], paste0('gt_', i+1), column_names[(i+1):length(column_names)])
    }
  }
  
  if (res_col_num > 6) {
    for (i in 6:(res_col_num - 1)) {
      column_names <- c(column_names, paste0('res_', i+1))
    }
  }
  setnames(match_dt, old = names(match_dt), new = column_names)
}

calculate_index <- function(dt, index_type) {
  col_ids <- as.character(colnames(dt))
  pos <- sapply(col_ids, function(col) {
    parts <- unlist(strsplit(col, split = "\\|"))
    if (length(parts) > 1) {
      pos_parts <- unlist(strsplit(parts[2], split = ":"))
      as.integer(pos_parts[2])
    } else {
      0
    }
  })

  strand <- sapply(col_ids, function(col) {
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
  if (index_type == "PDUI"){
    index_weight <- as.numeric(pos == max(pos))
  } else if (index_type == "PPUI") {
    index_weight <- as.numeric(pos == min(pos))
  } else if (index_type == "RWUI") {
    index_weight <- as.numeric(order(pos)) - 1
    index_weight <- index_weight / max(index_weight)
  } else if (index_type == "DWUI"){
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

wilcox_test_by_cell_type <- function(index_vec, cell_type_map){
    names(index_vec) <- cell_type_map[names(index_vec)]
    index_vec <- index_vec[!is.na(index_vec)]
    index_vec <- split(index_vec,names(index_vec))
    if (length(index_vec) > 1) {
        return(wilcox.test(as.numeric(index_vec[[1]], index_vec[[2]]))$p.value)
    } else {
        return(NA)
    }
}

dexseq_test_by_cell_type <- function(input_mtx, cell_type_map, split_num={
    env_split <- suppressWarnings(as.integer(Sys.getenv("APABENCHMARK_DEXSEQ_SPLIT_NUM", "")))
    if (is.na(env_split) || env_split < 1) 10 else env_split
}){
    mtx_for_dex <- mapply(function(mtx, te, cell_type_map){
            te <- gsub(":", "_", te)
            colnames(mtx) <- paste0(te, ":", seq_along(colnames(mtx)))
            rownames(mtx) <- cell_type_map[rownames(mtx)]
            mtx <- matrix(as.numeric(as.matrix(mtx)), nrow = nrow(mtx), ncol = ncol(mtx), dimnames = list(rownames(mtx), colnames(mtx)))
        },
        input_mtx,
        names(input_mtx),
        MoreArgs=list(cell_type_map = cell_type_map),
        SIMPLIFY = FALSE
    ) %>% do.call(cbind, .)

    cell_types <- unique(rownames(mtx_for_dex))

    ct_dex_1 <- mtx_for_dex[rownames(mtx_for_dex) == cell_types[1],] %>% split.data.frame(., rep(1:split_num, each=nrow(.) %/% split_num)) %>% lapply(., colSums) %>% do.call(cbind, .)
    ct_dex_2 <- mtx_for_dex[rownames(mtx_for_dex) == cell_types[2],] %>% split.data.frame(., rep(1:split_num, each=nrow(.) %/% split_num)) %>% lapply(., colSums) %>% do.call(cbind, .)

    colnames(ct_dex_1) <- paste0(cell_types[1], "_", 1:split_num)
    colnames(ct_dex_2) <- paste0(cell_types[2], "_", 1:split_num)

    dex_mtx <- cbind(ct_dex_1, ct_dex_2)
    sample_table <- data.frame(row.names = colnames(dex_mtx), condition=c(rep(c(cell_types[1]), split_num), rep(c(cell_types[2]), split_num)) )

    feature_id <- rownames(dex_mtx)
    group_id <- rownames(dex_mtx) %>% strsplit(., ":") %>% sapply(., function(x){x[1]})

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
    dxd = DEXSeq::estimateSizeFactors(dxd, locfunc = genefilter::shorth)
    dxd = DEXSeq::estimateDispersions(dxd)
    dxd = DEXSeq::testForDEU(dxd)
    dxd = DEXSeq::estimateExonFoldChanges(dxd)
    dxr = DEXSeq::DEXSeqResults(dxd)
    return(as.data.table(dxr))
}

fisher_test_by_cell_type <- function(mtx, cell_type_map){
    if (ncol(mtx) <= 1) {
        return(NA)
    }
    row_name <- as.vector(cell_type_map[rownames(mtx)])
    keep <- !is.na(row_name) & nzchar(row_name)
    if (sum(keep) < 2) {
        return(NA)
    }

    mtx <- apply(mtx[keep, , drop = FALSE], 2, as.numeric)
    if (is.null(dim(mtx))) {
        mtx <- matrix(mtx, ncol = 1)
    }
    rownames(mtx) <- row_name[keep]
    ct_mtx <- rowsum(mtx, group = rownames(mtx))

    # fisher.test requires at least two groups with non-zero row sums.
    ct_mtx <- ct_mtx[rowSums(ct_mtx) > 0, , drop = FALSE]
    if (nrow(ct_mtx) < 2 || ncol(ct_mtx) < 2) {
        return(NA)
    }

    fisher_p <- tryCatch(
        fisher.test(ct_mtx, simulate.p.value = TRUE)$p.value,
        error = function(e) NA
    )
    return(fisher_p)
}

calculate_diff <- function(mtx, cell_type_map, index_type){

    row_name <- cell_type_map[rownames(mtx)]
    mtx <-  apply(mtx, 2, as.numeric)
    rownames(mtx) <- row_name
    ct_mtx <- rowsum(mtx, group = rownames(mtx))
    col_ids <- as.character(colnames(ct_mtx))
    # calculate index
    pos <- sapply(col_ids, function(col) {
    parts <- unlist(strsplit(col, split = "\\|"))
    if (length(parts) > 1) {
    pos_parts <- unlist(strsplit(parts[2], split = ":"))
    as.integer(pos_parts[2])
    } else {
    0
    }
    })

    strand <- sapply(col_ids, function(col) {
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
    if (index_type == "PDUI"){
        index_weight <- as.numeric(pos == max(pos))
    } else if (index_type == "PPUI") {
        index_weight <- as.numeric(pos == min(pos))
    } else if (index_type == "RWUI") {
        index_weight <- as.numeric(order(pos)) - 1
        index_weight <- index_weight / max(index_weight)
    } else if (index_type == "DWUI"){
        index_weight <- pos / sum(pos)
    } else if (index_type == "MPRO"){
        ct_prop <- sweep(ct_mtx, 1, rowSums(ct_mtx), `/`)
        mpro <- max(abs(ct_prop[1,] - ct_prop[2,]))
        return(mpro)
    } else
    {
        stop("Index type not found, please check")
    }
  index_weight_matrix <- matrix(index_weight, nrow = nrow(ct_mtx), ncol = ncol(ct_mtx), byrow = TRUE)
  weighted_sum <- rowSums(ct_mtx * index_weight_matrix)
  row_sums <- rowSums(ct_mtx)
  index <- weighted_sum / row_sums
  return(abs(index[1] - index[2]))
}

calculate_apa_detect_performance <- function(gt_filter_group, pd_filter_group, filter_type_1, filter_type_2){
    gt <- intersect(gt_filter_group[[filter_type_1]], gt_filter_group[[filter_type_2]])
    pd <- intersect(pd_filter_group[[filter_type_1]], pd_filter_group[[filter_type_2]])
    precision <- length(intersect(gt, pd)) / length(pd)
    recall <- length(intersect(gt, pd)) / length(gt)
    f1 <- 2 * precision * recall / (precision + recall)
    result <- data.table(filter_type_1 = filter_type_1, 
                         filter_type_2 = filter_type_2,
                         precision = precision, 
                         recall = recall, 
                         f1 = f1,
                         gt_count = length(gt),
                         pd_count = length(pd))
}

perform_full_apa_analysis <- function(
    grouped_gt_mtx,
    grouped_pd_mtx,
    cell_type_map,
    sample,
    tool,
    protocol,
    # 新增参数化配置
    index_types = c("PDUI", "PPUI", "RWUI", "DWUI"),
    diff_types = c("PDUI", "PPUI", "RWUI", "DWUI", "MPRO"),
    p_thresholds = c(0.05, 0.01),
    log2fc_thresholds = c(0.25, 0.5, 0.75, 1, 1.25),
    index_diff_thresholds = c(0.1, 0.2, 0.3, 0.4, 0.5),
    adjustment_method = "fdr"
) {
    # 初始化结果存储
    apa_detect_performance <- data.table()
    
    # 添加参数校验
    stopifnot(
        is.list(grouped_gt_mtx),
        is.list(grouped_pd_mtx),
        !is.null(names(grouped_gt_mtx)),
        !is.null(names(grouped_pd_mtx))
    )
    
# statistical test
    index_types <- c("PDUI", "PPUI", "RWUI", "DWUI")
    gt_p_list <- list()
    pd_p_list <- list()
    suppressWarnings({
    for (index in index_types){
        gt_index <- lapply(grouped_gt_mtx, calculate_index, index_type = index)
        pd_index <- lapply(grouped_pd_mtx, calculate_index, index_type = index)
        gt_p_list[[paste0("wilcox_",index)]] <- lapply(gt_index, wilcox_test_by_cell_type, cell_type_map = cell_type_map) %>% p.adjust(method="fdr")
        pd_p_list[[paste0("wilcox_",index)]] <- lapply(pd_index, wilcox_test_by_cell_type, cell_type_map = cell_type_map) %>% p.adjust(method="fdr")
    }
    })

    gt_p_list[["fisher"]] <- lapply(grouped_gt_mtx, fisher_test_by_cell_type, cell_type_map = cell_type_map) %>% p.adjust(method="fdr")
    # pd_p_list[["fisher"]] <- lapply(grouped_pd_mtx, fisher_test_by_cell_type, cell_type_map = cell_type_map) %>% p.adjust(method="fdr")
    pd_p_list[["fisher"]] <- lapply(grouped_pd_mtx,
    function(mtx, cell_type_map){
        tryCatch({fisher_test_by_cell_type(mtx = mtx, cell_type_map = cell_type_map)},
        error = function(e) {
        print(e)
        return(NA)
    })}, cell_type_map = cell_type_map) %>% p.adjust(method="fdr")

    pd_dxr <- dexseq_test_by_cell_type(grouped_pd_mtx, cell_type_map)
    gt_dxr <- dexseq_test_by_cell_type(grouped_gt_mtx, cell_type_map)
    gt_log2fc_col <- grep("^log2fold_", names(gt_dxr), value = TRUE)[1]
    pd_log2fc_col <- grep("^log2fold_", names(pd_dxr), value = TRUE)[1]

    # calculate index difference
    diff_types <- c("PDUI", "PPUI", "RWUI", "DWUI", "MPRO")
    gt_diff_list <- list()
    pd_diff_list <- list()
    for (diff_type in diff_types){
        gt_diff_list[[diff_type]] <- lapply(grouped_gt_mtx, calculate_diff, cell_type_map = cell_type_map, index_type = diff_type)
        pd_diff_list[[diff_type]] <- lapply(grouped_pd_mtx, calculate_diff, cell_type_map = cell_type_map, index_type = diff_type)
    }

    # grep filter group
    p_thresholds <- c(0.05, 0.01)
    gt_filter_group <- list()
    pd_filter_group <- list()
    for (p_threshold in p_thresholds){
        for (stats in names(gt_p_list)){
            gt_filter_group[[paste0(stats, "_", p_threshold)]] <- gt_p_list[[stats]] %>% .[. < p_threshold] %>% names
            pd_filter_group[[paste0(stats, "_", p_threshold)]] <- pd_p_list[[stats]] %>% .[. < p_threshold] %>% names
        }
        gt_filter_group[[paste0("dexseq", "_", p_threshold)]] <- gt_dxr[padj < p_threshold, unique(groupID)] %>% as.vector %>% gsub(.,pattern = "_",replacement = ":")
        pd_filter_group[[paste0("dexseq", "_", p_threshold)]] <- pd_dxr[padj < p_threshold, unique(groupID)] %>% as.vector %>% gsub(.,pattern = "_",replacement = ":")
    }

    log2fc_thresholds <- c(0.25, 0.5, 0.75, 1, 1.25)
    for (log2fc_threshold in log2fc_thresholds){
        if (!is.na(gt_log2fc_col) && !is.na(pd_log2fc_col)) {
            gt_filter_group[[paste0("dexseq_log2fc", "_", log2fc_threshold)]] <- gt_dxr[abs(get(gt_log2fc_col)) > log2fc_threshold, unique(groupID)] %>% as.vector %>% gsub(.,pattern = "_",replacement = ":")
            pd_filter_group[[paste0("dexseq_log2fc", "_", log2fc_threshold)]] <- pd_dxr[abs(get(pd_log2fc_col)) > log2fc_threshold, unique(groupID)] %>% as.vector %>% gsub(.,pattern = "_",replacement = ":")
        } else {
            gt_filter_group[[paste0("dexseq_log2fc", "_", log2fc_threshold)]] <- character(0)
            pd_filter_group[[paste0("dexseq_log2fc", "_", log2fc_threshold)]] <- character(0)
        }
    }

    index_diff_thresholds <- c(0.1, 0.2, 0.3, 0.4, 0.5)
    for (index_diff_threshold in index_diff_thresholds){
        for (diff_type in names(gt_diff_list)){
            gt_filter_group[[paste0(diff_type, "_", index_diff_threshold)]] <- gt_diff_list[[diff_type]] %>% .[. > index_diff_threshold] %>% names
            pd_filter_group[[paste0(diff_type, "_", index_diff_threshold)]] <- pd_diff_list[[diff_type]] %>% .[. > index_diff_threshold] %>% names
        }
    }

    p_stats <- c("wilcox_PDUI", "wilcox_PPUI", "wilcox_RWUI", "wilcox_DWUI", "fisher", "dexseq")
    p_stats <- expand.grid(p_stats, p_thresholds, stringsAsFactors = FALSE) %>% apply(1, paste0, collapse = "_")
    diff_stats <- c("PDUI", "PPUI", "RWUI", "DWUI", "MPRO")
    diff_stats <- expand.grid(diff_stats, index_diff_thresholds, stringsAsFactors = FALSE) %>% apply(1, paste0, collapse = "_")
    fc_stats <- c("dexseq_log2fc")
    fc_stats <- expand.grid(fc_stats, as.character(log2fc_thresholds), stringsAsFactors = FALSE) %>% apply(1, paste0, collapse = "_")
    diff_stats <- c(diff_stats, fc_stats)

    apa_detect_performance <- data.table()
    for (filter_type_1 in p_stats){
        for (filter_type_2 in diff_stats){
            result <- calculate_apa_detect_performance(gt_filter_group, pd_filter_group, filter_type_1, filter_type_2)
            apa_detect_performance  <- rbind(apa_detect_performance , result)
        }
    }
    apa_detect_performance[, sample := sample]
    apa_detect_performance[, tool := tool]
    apa_detect_performance[, protocol := protocol]
    
    # 优化后的结果返回
    return(apa_detect_performance)
}



# Define the options
option_list <- list(
    make_option(c("--tool"), dest="tool", type="character"),
    make_option(c("--sample"), dest="sample", type="character"),
    make_option(c("--gt_pas"), dest="gt_pas", type="character"),
    make_option(c("--pd_pas"), dest="pd_pas", type="character"),
    make_option(c("--gt_mtx"), dest="gt_mtx", type="character"),
    make_option(c("--pd_mtx"), dest="pd_mtx", type="character"),
    make_option(c("--output_dir"), dest="output_dir", type="character")
)
# run sample
# Rscript calculate_benchmark_performance.R --tool scapatrap --sample SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2 --gt_pas /path/to/apabenchmark_final/data/sim_pas/mm10_sim_pas_gn5000_rep2.bed --pd_pas /path/to/apabenchmark_final/data/sim_bam_result/scapatrap/SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2/pas.bed --gt_mtx /path/to/apabenchmark_final/data/sim_bam/SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2.bam.expr.tsv --pd_mtx /path/to/apabenchmark_final/data/sim_bam_result/scapatrap/SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2/pas_counts.tsv --output_dir /path/to/apabenchmark_final/data/performance

# sample <- "SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2"
# tool <- "scapatrap"
# gt_pas_path <- "/path/to/apabenchmark_final/data/sim_pas/mm10_sim_pas_gn5000_rep2.bed"
# gt_mtx_path <- "/path/to/apabenchmark_final/data/sim_bam/SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2.bam.expr.tsv"
# pd_pas_path <- "/path/to/apabenchmark_final/data/sim_bam_result/scapatrap/SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2/pas.bed"
# pd_mtx_path <- "/path/to/apabenchmark_final/data/sim_bam_result/scapatrap/SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2/pas_counts.tsv"
# output_dir <- "/path/to/apabenchmark_final/data/performance"

opt <- parse_args(OptionParser(option_list=option_list))


tool <- opt$tool
sample <- opt$sample
gt_pas_path <- opt$gt_pas
pd_pas_path <- opt$pd_pas
gt_mtx_path <- opt$gt_mtx
pd_mtx_path <- opt$pd_mtx
output_dir <- opt$output_dir

protocol <- unlist(strsplit(sample, "_"))[1]

# gt 预处理
# read pas counts
gt_mtx <- fread(gt_mtx_path) %>% as.matrix
rownames(gt_mtx) <- gt_mtx[,1]
gt_mtx <- gt_mtx[,-1, drop = FALSE]
colnames(gt_mtx) <- as.character(colnames(gt_mtx))

# prepare gt pas data
gt_pas_dt <- fread(gt_pas_path)
setnames(gt_pas_dt, c("exon_id", "pas_type"), c("te_name", "te_type"))
gt_pas_dt[, pas_name := paste0(te_name, "|", chr, ":", start, ":", end, ":", strand)]
gt_pas_dt[, pas_num := .N, by = te_name]
gt_col_num <- ncol(gt_pas_dt)

# prepare gt te data
gt_te_dt <- gt_pas_dt[, c("te_name", "pas_num"), with = FALSE]
gt_te_dt <- gt_te_dt[!duplicated(te_name)]
gt_te_dt[, chr := sapply(strsplit(te_name, ":"), `[`, 3), by = .(te_name)]
gt_te_dt[, start := as.integer(sapply(strsplit(te_name, ":"), `[`, 4)), by = .(te_name)]
gt_te_dt[, end := as.integer(sapply(strsplit(te_name, ":"), `[`, 5)), by = .(te_name)]
gt_te_dt[, strand := sapply(strsplit(te_name, ":"), `[`, 6), by = .(te_name)]
gt_te_dt <- gt_te_dt[order(chr, start, end, strand)]
gt_te_dt <- gt_te_dt[, c("chr", "start", "end", "te_name", "pas_num", "strand"), with = FALSE]

# prepare pas-te maping
gt_pas_cols <- as.character(colnames(gt_mtx[,-1, drop = FALSE]))
gt_te_pas_map <- split(gt_pas_cols, sapply(gt_pas_cols, function(col) {
  parts <- unlist(strsplit(as.character(col), split = "\\|"))
  if (length(parts) >= 1) parts[1] else ""
}))

# prepare cell type mapping
cell_type_map <- as.list(gt_mtx[,"cell_type"])

# group gt mtx by te
grouped_gt_mtx <- lapply(gt_te_pas_map, function(cols) {
  gt_mtx[, cols, drop=FALSE]
})


# pd 数据预处理
pd_mtx <- fread(pd_mtx_path) %>% as.matrix
rownames(pd_mtx) <- pd_mtx[,1]
pd_mtx <- pd_mtx[rownames(gt_mtx),-1, drop = FALSE]
colnames(pd_mtx) <- as.character(colnames(pd_mtx))
pd_mtx <- pd_mtx[, colSums(pd_mtx != 0) > 0, drop = FALSE]

# prepare pd pas data
pd_pas_dt <- fread(pd_pas_path)
setnames(
    pd_pas_dt,
    c("V1", "V2", "V3", "V4", "V5","V6"),
    c("chr", "start", "end", "name", "pos", "strand")
    )
pd_pas_dt <- pd_pas_dt[order(chr, start, end)]

# remove duplicated pas and keep the most downstream one
pd_pas_dt[, rank := {
  strand_first <- strand[1]
  if (strand_first == "-") {
    (seq_along(strand)) - 1 
  } else {
    (.N - seq_along(strand))
  }
}, by = name]
pd_pas_dt <- pd_pas_dt[rank == 0]
pd_pas_dt[, pas_name := paste0("|", chr, ":", start, ":", end, ":", strand)]
pd_col_num <- ncol(pd_pas_dt)
pd_pas_dt <- pd_pas_dt[name %in% colnames(pd_mtx),]
colnames(pd_mtx) <- pd_pas_dt[match(colnames(pd_mtx), pd_pas_dt$name), pas_name]

# 定义窗口大小列表
window_sizes <- c(50, 100, 150, 200)

# 初始化结果表
match_performance <- data.table()
apa_detect_performance <- data.table()
pas_quantify_performance <- data.table()
# ----------------------------
# 第一部分：循环处理不同窗口大小的PAS匹配
# ----------------------------
# 定义窗口大小列表
window_sizes <- c(50, 100, 150, 200)

# 初始化结果表
match_performance <- data.table()
apa_detect_performance <- data.table()
pas_quantify_performance <- data.table()
# ----------------------------
# 第一部分：循环处理不同窗口大小的PAS匹配
# ----------------------------
for (window_size in window_sizes) {
    # 生成临时文件（每次循环使用新的临时文件）
    gt_pas_bed_tmp <- tempfile(fileext = ".bed")
    predict_pas_tmp <- tempfile(fileext = ".bed")
    fwrite(gt_pas_dt, gt_pas_bed_tmp, sep = "\t", quote = FALSE, col.names = FALSE)
    fwrite(pd_pas_dt, predict_pas_tmp, sep = "\t", quote = FALSE, col.names = FALSE)
    output_file <- tempfile()
    
    # 执行bedtools命令（注意窗口参数为总大小的一半）
    bedtools_command <- sprintf(
        "%s window -sw -sm -a %s -b %s -w %d",
        bedtools_bin,
        gt_pas_bed_tmp,
        predict_pas_tmp,
        window_size / 2  # 例如：200的窗口实际左右各扩展100
    )
    system(paste(bedtools_command, ">", output_file))
    
    # 读取结果并清理临时文件
    match_dt <- fread(output_file)
    unlink(c(gt_pas_bed_tmp, predict_pas_tmp, output_file))
    if (nrow(match_dt) == 0) {
        match_performance <- rbind(
            match_performance,
            data.table(
                sample = sample,
                tool = tool,
                protocol = protocol,
                match_type = paste0("pas_", window_size / 2),  # 关键修改：动态命名
                tp = 0,
                fp = 0,
                fn = nrow(gt_pas_dt),
                precision = 0,
                recall = 0,
                f1 = 0
            ) 
        )
        current_pas_quantify_performance <- data.table(
          sample = sample,
          tool = tool,
          protocol = protocol,
          match_type = paste0("pas_", window_size / 2),  # 关键修改：动态命名
          cor_pas = 0,
          rmse_pas = 0,
          mae_pas = 0,
          mape_pas = 0,
          rmse_pas_ct = 0,
          mae_pas_ct = 0,
          mape_pas_ct = 0
        )
        pas_quantify_performance <- rbind(pas_quantify_performance, current_pas_quantify_performance)
        next  # 跳过当前循环，进入下一个窗口大小的处理
    }
    # 后续处理（与原代码一致）
    rename_matched_columns(match_dt, gt_col_num, pd_col_num)
    match_dt[, distance := abs(gt_start - res_start)]
    match_dt[, best_match := {
        min_distance <- min(distance)
        distance == min_distance
    }, by = res_name]
    match_dt <- match_dt[best_match == TRUE]
    match_dt[, best_match := {
        min_distance <- min(distance)
        distance == min_distance
    }, by = gt_9]
    
    # 合并数据
    pas_merged_dt <- merge(
        gt_pas_dt,
        match_dt,
        by.x = c("chr", "start", "end", "strand"),
        by.y = c("gt_chr", "gt_start", "gt_end", "gt_strand"),
        all.x = TRUE
    )
    gt_cols <- names(pas_merged_dt)[grepl("^gt_", names(pas_merged_dt))]
    pas_merged_dt[, (gt_cols) := NULL]
    pas_merged_dt[, all_match := !is.na(best_match)]
    pas_merged_dt[is.na(best_match), best_match := FALSE]
    setnames(pas_merged_dt, paste0("res_", pd_col_num), "pd_pas_name")
    
    # 计算当前窗口的PAS性能指标
    pas_match_tp <- pas_merged_dt[all_match == TRUE, .N]
    pas_match_fp <- pd_pas_dt[, .N] - pas_match_tp
    pas_match_fn <- gt_pas_dt[, .N] - pas_match_tp
    pas_match_precision <- pas_match_tp / (pas_match_tp + pas_match_fp)
    pas_match_recall <- pas_match_tp / (pas_match_tp + pas_match_fn)
    pas_match_f1 <- 2 * pas_match_precision * pas_match_recall / (pas_match_precision + pas_match_recall)
    
    # 将结果添加到性能表，使用 pas_50/pas_100 等命名
    match_performance <- rbind(
        match_performance,
        data.table(
            sample = sample,
            tool = tool,
            protocol = protocol,
            match_type = paste0("pas_", window_size / 2),  # 关键修改：动态命名
            tp = pas_match_tp,
            fp = pas_match_fp,
            fn = pas_match_fn,
            precision = pas_match_precision,
            recall = pas_match_recall,
            f1 = pas_match_f1
        )
    )

    # 处理te map
    current_pas_map <- pas_merged_dt[all_match == TRUE, .(
        te_name, 
        gt_pas_name=pas_name,
        pd_pas_name
    )]

    current_pas_map[, pd_pas_count := .N, by = te_name]
    current_pas_map[, gt_pas_count := uniqueN(gt_pas_name), by = te_name]

    multipd_pd_te_pas_map <- split(current_pas_map[(pd_pas_count > 1), pd_pas_name], current_pas_map[(pd_pas_count > 1), te_name])

    grouped_pd_mtx_mapped_by_pd <- lapply(multipd_pd_te_pas_map, function(cols) {
    pd_mtx[, cols, drop=FALSE]
    })

    if (length(grouped_pd_mtx_mapped_by_pd) == 0) { 
        apa_detect_performance_mapped_by_pd <- data.table(
            sample = sample,
            tool = tool,
            protocol = protocol,
            match_type = paste0("pd_", window_size / 2)
        )
    } else {
        grouped_pd_mtx_mapped_by_pd <- grouped_pd_mtx_mapped_by_pd[sapply(grouped_pd_mtx_mapped_by_pd, function(x) !is.null(x) && ncol(x) > 0)]
        if (length(grouped_pd_mtx_mapped_by_pd) == 0) {
            apa_detect_performance_mapped_by_pd <- data.table(
            sample = sample,
            tool = tool,
            protocol = protocol,
            match_type = paste0("pd_", window_size / 2)
            )
        } else {
            apa_detect_performance_mapped_by_pd <- perform_full_apa_analysis(grouped_gt_mtx = grouped_gt_mtx, grouped_pd_mtx = grouped_pd_mtx_mapped_by_pd, cell_type_map = cell_type_map, sample = sample, tool = tool, protocol = protocol)
            apa_detect_performance_mapped_by_pd$match_type <- paste0("pd_", window_size / 2)
        }
    }

    gt_pas_mapping <- current_pas_map[
        gt_pas_count > 1,
        .(pd_pas_list = list(unique(pd_pas_name))),
        by = .(te_name, gt_pas_name)
    ]
    grouped_pd_mtx_mapped_by_gt <- lapply(split(gt_pas_mapping, by = "te_name"), function(te_data) {
        # 预先生成所有有效列名
        all_pd_cols <- unique(unlist(te_data$pd_pas_list))
        valid_cols <- intersect(all_pd_cols, colnames(pd_mtx))
        
        # 提取子矩阵并转换为数值矩阵
        sub_mtx <- as.matrix(pd_mtx[, valid_cols, drop = FALSE])
        storage.mode(sub_mtx) <- "numeric"  # 强制转换为数值类型
        # 创建结果矩阵
        result_mtx <- matrix(
            nrow = nrow(sub_mtx),
            ncol = nrow(te_data),
            dimnames = list(rownames(sub_mtx), te_data$gt_pas_name)
        )
        # 向量化聚合操作
        for (i in seq_len(nrow(te_data))) {
            target_cols <- intersect(te_data$pd_pas_list[[i]], valid_cols)
            if (length(target_cols) > 0) {
            result_mtx[, i] <- rowSums(sub_mtx[, target_cols, drop = FALSE], na.rm = TRUE)
            } else {
            result_mtx[, i] <- NA_real_
            }
        }
        
        # 移除全NA的列
        result_mtx[, colSums(is.na(result_mtx)) < nrow(result_mtx), drop = FALSE]
    })


    if (length(grouped_pd_mtx_mapped_by_gt) == 0) {
        apa_detect_performance_mapped_by_gt <- data.table(
            sample = sample,
            tool = tool,
            protocol = protocol,
            match_type = paste0("gt_", window_size / 2)
        )
    } else {
        grouped_pd_mtx_mapped_by_gt <- grouped_pd_mtx_mapped_by_gt[sapply(grouped_pd_mtx_mapped_by_gt, function(x) !is.null(x) && ncol(x) > 0)]
        if (length(grouped_pd_mtx_mapped_by_gt) == 0) {
            apa_detect_performance_mapped_by_gt <- data.table(
            sample = sample,
            tool = tool,
            protocol = protocol,
            match_type = paste0("gt_", window_size / 2)
        )
        } else {
            apa_detect_performance_mapped_by_gt <- perform_full_apa_analysis(grouped_gt_mtx = grouped_gt_mtx, grouped_pd_mtx = grouped_pd_mtx_mapped_by_gt, cell_type_map = cell_type_map, sample = sample, tool = tool, protocol = protocol)
            apa_detect_performance_mapped_by_gt$match_type <- paste0("gt_", window_size / 2)
        }

    }


    apa_detect_performance <- rbind(apa_detect_performance, apa_detect_performance_mapped_by_pd)
    apa_detect_performance <- rbind(apa_detect_performance, apa_detect_performance_mapped_by_gt)


    # calculate quantification performance
    mapped_pd_mtx <- pd_mtx %>% apply(., 2,as.numeric)
    colnames(mapped_pd_mtx) <- pas_merged_dt[match(colnames(pd_mtx), pas_merged_dt$pd_pas_name), pas_name]
    rownames(mapped_pd_mtx) <- rownames(pd_mtx)
    mapped_pd_mtx <- mapped_pd_mtx %>% t %>% rowsum(., group=colnames(mapped_pd_mtx)) %>% t
    mapped_pd_mtx <- mapped_pd_mtx[, !is.na(colnames(mapped_pd_mtx)), drop=FALSE]

    mapped_gt_mtx <- gt_mtx %>% apply(., 2, as.numeric)
    colnames(mapped_gt_mtx) <- colnames(gt_mtx)
    rownames(mapped_gt_mtx) <- rownames(gt_mtx)
    mapped_gt_mtx <- mapped_gt_mtx[,-1, drop = FALSE]

    shared_cells <- intersect(rownames(mapped_gt_mtx), rownames(mapped_pd_mtx))
    shared_pas <- intersect(colnames(mapped_gt_mtx), colnames(mapped_pd_mtx))

    if (length(shared_cells) == 0 || length(shared_pas) == 0) {
      cor_pas <- 0
      rmse_pas <- 0
      mae_pas <- 0
      mape_pas <- 0
      rmse_pas_ct <- 0
      mae_pas_ct <- 0
      mape_pas_ct <- 0
    } else {
      mapped_gt_mtx <- mapped_gt_mtx[shared_cells, shared_pas, drop = FALSE]
      mapped_pd_mtx <- mapped_pd_mtx[shared_cells, shared_pas, drop = FALSE]

      mapped_ct_gt_mtx <- mapped_gt_mtx
      mapped_ct_pd_mtx <- mapped_pd_mtx
      rownames(mapped_ct_gt_mtx) <- cell_type_map[rownames(mapped_gt_mtx)]
      rownames(mapped_ct_pd_mtx) <- cell_type_map[rownames(mapped_pd_mtx)]

      keep_gt <- !is.na(rownames(mapped_ct_gt_mtx)) & nzchar(rownames(mapped_ct_gt_mtx))
      keep_pd <- !is.na(rownames(mapped_ct_pd_mtx)) & nzchar(rownames(mapped_ct_pd_mtx))
      mapped_ct_gt_mtx <- mapped_ct_gt_mtx[keep_gt, , drop = FALSE]
      mapped_ct_pd_mtx <- mapped_ct_pd_mtx[keep_pd, , drop = FALSE]

      if (nrow(mapped_ct_gt_mtx) == 0 || nrow(mapped_ct_pd_mtx) == 0) {
        rmse_pas_ct <- 0
        mae_pas_ct <- 0
        mape_pas_ct <- 0
      } else {
        mapped_ct_gt_mtx <- rowsum(mapped_ct_gt_mtx, group = rownames(mapped_ct_gt_mtx))
        mapped_ct_pd_mtx <- rowsum(mapped_ct_pd_mtx, group = rownames(mapped_ct_pd_mtx))
        ct_rows <- intersect(rownames(mapped_ct_gt_mtx), rownames(mapped_ct_pd_mtx))
        ct_cols <- intersect(colnames(mapped_ct_gt_mtx), colnames(mapped_ct_pd_mtx))
        if (length(ct_rows) == 0 || length(ct_cols) == 0) {
          rmse_pas_ct <- 0
          mae_pas_ct <- 0
          mape_pas_ct <- 0
        } else {
          mapped_ct_gt_mtx <- mapped_ct_gt_mtx[ct_rows, ct_cols, drop = FALSE]
          mapped_ct_pd_mtx <- mapped_ct_pd_mtx[ct_rows, ct_cols, drop = FALSE]
          rmse_pas_ct <- sqrt(mean((mapped_ct_gt_mtx - mapped_ct_pd_mtx)^2))
          mae_pas_ct <- mean(abs(mapped_ct_gt_mtx - mapped_ct_pd_mtx))
          mape_pas_ct <- mean(ifelse(mapped_ct_gt_mtx == 0, 0, abs(mapped_ct_gt_mtx - mapped_ct_pd_mtx) / mapped_ct_gt_mtx))
        }
      }

      cor_pas <- mapply(cor, as.data.frame(mapped_pd_mtx), as.data.frame(mapped_gt_mtx), MoreArgs = list(method="pearson")) %>% mean
      rmse_pas <- sqrt(mean((mapped_gt_mtx - mapped_pd_mtx)^2))
      mae_pas <- mean(abs(mapped_gt_mtx - mapped_pd_mtx))
      mape_pas <- mean(ifelse(mapped_gt_mtx == 0, 0, abs(mapped_gt_mtx - mapped_pd_mtx) / mapped_gt_mtx))
      mape_pas <- mean(mape_pas)
    }

    current_pas_quantify_performance <- data.table(
      sample = sample,
      tool = tool,
      protocol = protocol,
      match_type = paste0("pas_", window_size / 2),  # 关键修改：动态命名
      cor_pas = cor_pas,
      rmse_pas = rmse_pas,
      mae_pas = mae_pas,
      mape_pas = mape_pas,
      rmse_pas_ct = rmse_pas_ct,
      mae_pas_ct = mae_pas_ct,
      mape_pas_ct = mape_pas_ct
    )
    pas_quantify_performance <- rbind(pas_quantify_performance, current_pas_quantify_performance)
}




# ----------------------------
# 第二部分：处理TE（循环外独立执行）
# ----------------------------
gt_te_bed_tmp <- tempfile(fileext = ".bed")
predict_pas_tmp <- tempfile(fileext = ".bed")
fwrite(gt_te_dt, gt_te_bed_tmp, sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(pd_pas_dt, predict_pas_tmp, sep = "\t", quote = FALSE, col.names = FALSE)
output_file <- tempfile()
bedtools_command <- sprintf(
    "%s intersect -s -a %s -b %s -wa -wb",
    bedtools_bin,
    gt_te_bed_tmp,
    predict_pas_tmp
)
system(paste(bedtools_command, ">", output_file))
te_pas_match <- fread(output_file)
unlink(c(gt_te_bed_tmp, predict_pas_tmp, output_file))
rename_matched_columns(te_pas_match, 6, pd_col_num)
te_pas_match[, pd_pas_count := .N, by = gt_name]
setnames(te_pas_match, paste0("res_", pd_col_num), "pd_pas_name")

# 计算TE性能指标
te_match_tp <- te_pas_match[(gt_score > 1) & (pd_pas_count > 1), length(unique(gt_name))]
te_match_fp <- te_pas_match[(gt_score == 1) & (pd_pas_count > 1), length(unique(gt_name))]
te_match_fn <- gt_te_dt[pas_num > 1, length(unique(te_name))] - te_match_tp
te_match_precision <- te_match_tp / (te_match_tp + te_match_fp)
te_match_recall <- te_match_tp / (te_match_tp + te_match_fn)
te_match_f1 <- 2 * te_match_precision * te_match_recall / (te_match_precision + te_match_recall)

# 将TE结果添加到性能表
match_performance <- rbind(
    match_performance,
    data.table(
        sample = sample,
        tool = tool,
        protocol = protocol,
        match_type = "te",  # TE独立命名
        tp = te_match_tp,
        fp = te_match_fp,
        fn = te_match_fn,
        precision = te_match_precision,
        recall = te_match_recall,
        f1 = te_match_f1
    )
)

pd_te_pas_map <- split(te_pas_match[(!is.na(res_name)) & (pd_pas_count > 1), pd_pas_name], te_pas_match[(!is.na(res_name)) & (pd_pas_count > 1), gt_name])

grouped_pd_mtx_mapped_by_te <- lapply(pd_te_pas_map, function(cols) {
  pd_mtx[, cols, drop=FALSE]
})

apa_detect_performance_mapped_by_te <- perform_full_apa_analysis(grouped_gt_mtx = grouped_gt_mtx, grouped_pd_mtx = grouped_pd_mtx_mapped_by_te, cell_type_map = cell_type_map, sample = sample, tool = tool, protocol = protocol)
apa_detect_performance_mapped_by_te$match_type <- "te"

apa_detect_performance <- rbind(
    apa_detect_performance,
    apa_detect_performance_mapped_by_te
)


# 定义所有窗口类型：pas_xx/pd_xx/gt_xx
match_types <- c("pd", "gt")
window_sizes_half <- window_sizes / 2  # e.g., 25,50,75,100

# 定义统计测试相关参数
p_stats <- c("wilcox_PDUI", "wilcox_PPUI", "wilcox_RWUI", "wilcox_DWUI", "fisher", "dexseq")
p_thresholds <- c(0.05, 0.01)
diff_stats <- c("PDUI", "PPUI", "RWUI", "DWUI", "MPRO")
index_diff_thresholds <- c(0.1, 0.2, 0.3, 0.4, 0.5)
log2fc_thresholds <- c(0.25, 0.5, 0.75, 1, 1.25)

# 生成完整filter_type组合
full_filter_types <- expand.grid(
  filter_type_1 = paste0(rep(p_stats, each = length(p_thresholds)), "_", rep(p_thresholds, length(p_stats))),
  filter_type_2 = c(
    paste0(rep(diff_stats, length(index_diff_thresholds)), "_", rep(index_diff_thresholds, each = length(diff_stats))),
    paste0("dexseq_log2fc_", log2fc_thresholds)
  ),
  stringsAsFactors = FALSE
)

# 生成原有窗口类型的组合
window_based_meta <- expand.grid(
    sample = sample,
    tool = tool,
    protocol = protocol,
    match_type = paste0(rep(c("pd", "gt"), each = length(window_sizes_half)), "_", window_sizes_half),
    stringsAsFactors = FALSE
)
# 添加单独的 "te" 类型（不关联窗口）
additional_te_block <- data.frame(
    sample = sample,
    tool = tool,
    protocol = protocol,
    match_type = "te",
    stringsAsFactors = FALSE
)
# 合并所有可能的 match_type 组合
meta_combos <- rbind(
    window_based_meta,
    additional_te_block
)


apa_complete_meta <- merge(
    meta_combos,  # 含 sample/tool/protocol/match_type
    full_filter_types,  # 含 filter_type_1 和 filter_type_2
    by = character(),  # 交叉连接所有行
    allow.cartesian = TRUE
)
apa_detect_performance_full <- merge(
    apa_complete_meta,
    apa_detect_performance,
    by = c("sample", "tool", "protocol", "match_type", "filter_type_1", "filter_type_2"),
    all.x = TRUE
)


# fill na with zero
apa_detect_performance_full[is.na(apa_detect_performance_full)] <- 0

# save results

tool_output_dir <- file.path(output_dir, tool)
if (!dir.exists(tool_output_dir)) {
  dir.create(tool_output_dir, recursive = TRUE, showWarnings = FALSE)
}

output_path_path <- paste(output_dir, tool, sample, sep = "/")
match_performance_path <- paste(output_path_path, "match_performance.tsv", sep = "_")
pas_quantify_performance_path <- paste(output_path_path, "pas_quantify_performance.tsv", sep = "_")
de_apa_performance_path <- paste(output_path_path, "de_apa_performance.tsv", sep = "_")
# match_pas_path <- paste(output_path_path, "match_pas.csv", sep = "_")
# match_te_path <- paste(output_path_path, "match_te.csv", sep = "_")

fwrite(match_performance, match_performance_path, sep = "\t", quote = FALSE, col.names = TRUE)
fwrite(apa_detect_performance_full, de_apa_performance_path, sep = "\t", quote = FALSE, col.names = TRUE)
fwrite(pas_quantify_performance, pas_quantify_performance_path, sep = "\t", quote = FALSE, col.names = TRUE)
# fwrite(data.frame(value = pas_merged_dt[all_match == TRUE, unique(pas_name)], stringsAsFactors = FALSE), match_pas_path, sep = ",", quote = FALSE, col.names = TRUE)
# fwrite(data.frame(value = te_pas_match[pd_pas_count > 1, unique(gt_name)], stringsAsFactors = FALSE), match_te_path, sep = ",", quote = FALSE, col.names = TRUE)
