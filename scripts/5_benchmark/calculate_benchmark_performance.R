library(Seurat)
library(DEXSeq)
library(data.table)
library(magrittr)
library(dplyr)
library(optparse)

bedtools_path <- "/root/miniconda3/envs/dexseq/bin/"

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

dexseq_test_by_cell_type <- function(input_mtx, cell_type_map, split_num=6){
    mtx_for_dex <- mapply(function(mtx, te, cell_type_map){
            te <- gsub(":", "_", te)
            colnames(mtx) <- paste0(te, ":", seq_along(colnames(mtx)))
            rownames(mtx) <- cell_type_map[rownames(mtx)]
            mtx <- matrix(as.numeric(as.matrix(mtx)), nrow = nrow(mtx), ncol = ncol(mtx), dimnames = list(rownames(mtx), colnames(mtx)))
        },
        input_mtx,
        names(input_mtx),
        MoreArgs=list(cell_type_map = cell_type_map)
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
    mtx <- apply(mtx, 2, as.numeric)
    rownames(mtx) <- row_name
    ct_mtx <- rowsum(mtx, group = rownames(mtx))
    fisher_p <- fisher.test(ct_mtx, simulate.p.value=TRUE)$p
    return(fisher_p)
}

calculate_diff <- function(mtx, cell_type_map, index_type){

    row_name <- cell_type_map[rownames(mtx)]
    mtx <-  apply(mtx, 2, as.numeric)
    rownames(mtx) <- row_name
    ct_mtx <- rowsum(mtx, group = rownames(mtx))
    # calculate index
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
                         f1 = f1)
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
# Rscript calculate_benchmark_performance.R --tool scapatrap --sample SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2 --gt_pas /root/apabenchmark/data/sim_pas/mm10_sim_pas_gn5000_rep2.bed --pd_pas /root/apabenchmark/data/sim_bam_result/scapatrap/SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2/pas.bed --gt_mtx /root/apabenchmark/data/sim_bam/SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2.bam.expr.tsv --pd_mtx /root/apabenchmark/data/sim_bam_result/scapatrap/SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2/pas_counts.tsv --output_dir /root/apabenchmark/data/performance

# sample <- "SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2"
# tool <- "scapatrap"
# gt_pas_path <- "/root/apabenchmark/data/sim_pas/mm10_sim_pas_gn5000_rep2.bed"
# gt_mtx_path <- "/root/apabenchmark/data/sim_bam/SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2.bam.expr.tsv"
# pd_pas_path <- "/root/apabenchmark/data/sim_bam_result/scapatrap/SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2/pas.bed"
# pd_mtx_path <- "/root/apabenchmark/data/sim_bam_result/scapatrap/SpatialTranscriptomics_mouse_olfactorybulb_Rep11MOB_pas2_gn5000_rep2/pas_counts.tsv"
# output_dir <- "/root/apabenchmark/data/performance"

opt <- parse_args(OptionParser(option_list=option_list))

tool <- opt$tool
sample <- opt$sample
gt_pas_path <- opt$gt_pas
pd_pas_path <- opt$pd_pas
gt_mtx_path <- opt$gt_mtx
pd_mtx_path <- opt$pd_mtx
output_dir <- opt$output_dir

protocol <- unlist(strsplit(sample, "_"))[1]

# read pas counts
gt_mtx <- fread(gt_mtx_path) %>% as.matrix 
rownames(gt_mtx) <- gt_mtx[,1]
gt_mtx <- gt_mtx[,-1]

pd_mtx <- fread(pd_mtx_path) %>% as.matrix
rownames(pd_mtx) <- pd_mtx[,1]
pd_mtx <- pd_mtx[rownames(gt_mtx),-1]
pd_mtx <- pd_mtx[, colSums(pd_mtx != 0) > 0]

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

#match pas
window_size <- 200

gt_pas_bed_tmp <- tempfile(fileext = ".bed")
predict_pas_tmp <- tempfile(fileext = ".bed")
fwrite(gt_pas_dt, gt_pas_bed_tmp, sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(pd_pas_dt, predict_pas_tmp, sep = "\t", quote = FALSE, col.names = FALSE)
output_file <- tempfile()
bedtools_command <- sprintf("%sbedtools window -sw -a %s -b %s -w %d", bedtools_path, gt_pas_bed_tmp, predict_pas_tmp, window_size / 2)
system(paste(bedtools_command, ">", output_file))
match_dt <- fread(output_file)
unlink(gt_pas_bed_tmp)
unlink(predict_pas_tmp)
unlink(output_file)
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

pas_merged_dt <- merge(
  gt_pas_dt,
  match_dt,
  by.x = c("chr", "start", "end","strand"),
  by.y = c("gt_chr", "gt_start", "gt_end", "gt_strand"),
  all.x = TRUE
)

gt_cols <- names(pas_merged_dt)[grepl("^gt_", names(pas_merged_dt))]
pas_merged_dt[, (gt_cols) := NULL]
pas_merged_dt[, all_match := {
   ! is.na(best_match)
}]
pas_merged_dt[is.na(best_match), best_match := FALSE]
setnames(pas_merged_dt, paste0("res_", pd_col_num), "pd_pas_name")


# match te
gt_te_bed_tmp <- tempfile(fileext = ".bed")
predict_pas_tmp <- tempfile(fileext = ".bed")
fwrite(gt_te_dt, gt_te_bed_tmp, sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(pd_pas_dt, predict_pas_tmp, sep = "\t", quote = FALSE, col.names = FALSE)
output_file <- tempfile()
bedtools_command <- sprintf("%sbedtools intersect -s -a %s -b %s -wa -wb", bedtools_path, gt_te_bed_tmp, predict_pas_tmp)
system(paste(bedtools_command, ">", output_file))
te_pas_match <- fread(output_file)
unlink(gt_te_bed_tmp)
unlink(predict_pas_tmp)
unlink(output_file)
rename_matched_columns(te_pas_match, 6, pd_col_num)
te_pas_match[, pd_pas_count := .N, by = gt_name]
setnames(te_pas_match, paste0("res_", pd_col_num), "pd_pas_name")


# calculate match performance
pas_match_tp <- pas_merged_dt[all_match == TRUE, .N]
pas_match_fp <- pd_pas_dt[,.N] - pas_match_tp
pas_match_fn <- gt_pas_dt[,.N] - pas_match_tp

pas_match_precision <- pas_match_tp / (pas_match_tp + pas_match_fp)
pas_match_recall <- pas_match_tp / (pas_match_tp + pas_match_fn)
pas_match_f1 <- 2 * pas_match_precision * pas_match_recall / (pas_match_precision + pas_match_recall)

te_match_tp <- te_pas_match[(gt_score > 1) & (pd_pas_count > 1), length(unique(gt_name))]
te_match_fp <- te_pas_match[(gt_score == 1) & (pd_pas_count > 1), length(unique(gt_name))]
te_match_fn <- gt_te_dt[pas_num > 1, length(unique(te_name))] - te_match_tp

te_match_precision <- te_match_tp / (te_match_tp + te_match_fp)
te_match_recall <- te_match_tp / (te_match_tp + te_match_fn)
te_match_f1 <- 2 * te_match_precision * te_match_recall / (te_match_precision + te_match_recall)

match_performance <- data.table()
match_performance <- rbind(
  match_performance,
  data.table(
    sample = sample,
    tool = tool,
    protocol = protocol,
    match_type = "pas",
    tp = pas_match_tp,
    fp = pas_match_fp,
    fn = pas_match_fn,
    precision = pas_match_precision,
    recall = pas_match_recall,
    f1 = pas_match_f1
  )
)
match_performance <- rbind(
  match_performance,
  data.table(
    sample = sample,
    tool = tool,
    protocol = protocol,
    match_type = "te",
    tp = te_match_tp,
    fp = te_match_fp,
    fn = te_match_fn,
    precision = te_match_precision,
    recall = te_match_recall,
    f1 = te_match_f1
  )
)


# prepare pas-te maping
gt_te_pas_map <- split(colnames(gt_mtx[,-1]), sapply(colnames(gt_mtx[,-1]), function(col) {
  unlist(strsplit(col,split = "\\|"))[1]
}))
pd_te_pas_map <- split(te_pas_match[(!is.na(res_name)) & (pd_pas_count > 1), pd_pas_name], te_pas_match[(!is.na(res_name)) & (pd_pas_count > 1), gt_name])

# prepare cell type mapping
cell_type_map <- as.list(gt_mtx[,"cell_type"])

grouped_gt_mtx <- lapply(gt_te_pas_map, function(cols) {
  gt_mtx[, cols, drop=FALSE]
})
grouped_pd_mtx <- lapply(pd_te_pas_map, function(cols) {
  pd_mtx[, cols, drop=FALSE]
})

# statistical test
index_types <- c("PDUI", "PPUI", "RWUI", "DWUI")
gt_p_list <- list()
pd_p_list <- list()
suppressWarnings({
for (index in index_types){
    gt_index <- lapply(grouped_gt_mtx, calculate_index, index_type = index)
    pd_index <- lapply(grouped_pd_mtx, calculate_index, index_type = index)
    index <- paste0("wilcox_",index)
    gt_p_list[[index]] <- lapply(gt_index, wilcox_test_by_cell_type, cell_type_map = cell_type_map) %>% p.adjust(method="fdr") 
    pd_p_list[[index]] <- lapply(pd_index, wilcox_test_by_cell_type, cell_type_map = cell_type_map) %>% p.adjust(method="fdr")
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
    gt_filter_group[[paste0("dexseq_log2fc", "_", log2fc_threshold)]] <- gt_dxr[abs(log2fold_T2_T1) > log2fc_threshold, unique(groupID)] %>% as.vector %>% gsub(.,pattern = "_",replacement = ":")
    pd_filter_group[[paste0("dexseq_log2fc", "_", log2fc_threshold)]] <- pd_dxr[abs(log2fold_T2_T1) > log2fc_threshold, unique(groupID)] %>% as.vector %>% gsub(.,pattern = "_",replacement = ":")
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


# calculate quantification performance
mapped_pd_mtx <- pd_mtx %>% apply(., 2,as.numeric)
colnames(mapped_pd_mtx) <- pas_merged_dt[match(colnames(pd_mtx), pas_merged_dt$pd_pas_name), pas_name]
rownames(mapped_pd_mtx) <- rownames(pd_mtx)
mapped_pd_mtx <- mapped_pd_mtx %>% t %>% rowsum(., group=colnames(mapped_pd_mtx)) %>% t
mapped_pd_mtx <- mapped_pd_mtx[, !is.na(colnames(mapped_pd_mtx))]

mapped_gt_mtx <- gt_mtx %>% apply(., 2, as.numeric)
colnames(mapped_gt_mtx) <- colnames(gt_mtx)
rownames(mapped_gt_mtx) <- rownames(gt_mtx)
mapped_gt_mtx <- mapped_gt_mtx[,-1]
mapped_gt_mtx <- mapped_gt_mtx[rownames(mapped_pd_mtx), colnames(mapped_pd_mtx)]

mapped_ct_gt_mtx <- mapped_gt_mtx
mapped_ct_pd_mtx <- mapped_pd_mtx
rownames(mapped_ct_gt_mtx) <- cell_type_map[rownames(mapped_gt_mtx)]
rownames(mapped_ct_pd_mtx) <- cell_type_map[rownames(mapped_pd_mtx)]
mapped_ct_gt_mtx <- rowsum(mapped_ct_gt_mtx, group = rownames(mapped_ct_gt_mtx))
mapped_ct_pd_mtx <- rowsum(mapped_ct_pd_mtx, group = rownames(mapped_ct_pd_mtx))

cor_pas <- mapply(cor, as.data.frame(mapped_pd_mtx), as.data.frame(mapped_gt_mtx), MoreArgs = list(method="pearson")) %>% mean
rmse_pas <- sqrt(mean((mapped_gt_mtx - mapped_pd_mtx)^2))
mae_pas <- mean(abs(mapped_gt_mtx - mapped_pd_mtx))
mape_pas <- mean(ifelse(mapped_gt_mtx == 0, 0, abs(mapped_gt_mtx - mapped_pd_mtx) / mapped_gt_mtx))
mape_pas <- mean(mape_pas)
rmse_pas_ct <- sqrt(mean((mapped_ct_gt_mtx - mapped_ct_pd_mtx)^2))
mae_pas_ct <- mean(abs(mapped_ct_gt_mtx - mapped_ct_pd_mtx))
mape_pas_ct <- mean(ifelse(mapped_ct_gt_mtx == 0, 0, abs(mapped_ct_gt_mtx - mapped_ct_pd_mtx) / mapped_ct_gt_mtx))

pas_quantify_performance <- data.table(
  sample = sample,
  tool = tool,
  protocol = protocol,
  cor_pas = cor_pas,
  rmse_pas = rmse_pas,
  mae_pas = mae_pas,
  mape_pas = mape_pas,
  rmse_pas_ct = rmse_pas_ct,
  mae_pas_ct = mae_pas_ct,
  mape_pas_ct = mape_pas_ct
)

# save results

if (!file.exists(paste(output_dir, tool, sep = "/"))) {
  dir.create(paste(output_dir, tool, sep = "/"))
}

output_path_path <- paste(output_dir, tool, sample, sep = "/")
match_performance_path <- paste(output_path_path, "match_performance.tsv", sep = "_")
pas_quantify_performance_path <- paste(output_path_path, "pas_quantify_performance.tsv", sep = "_")
de_apa_performance_path <- paste(output_path_path, "de_apa_performance.tsv", sep = "_")
match_pas_path <- paste(output_path_path, "match_pas.csv", sep = "_")
match_te_path <- paste(output_path_path, "match_te.csv", sep = "_")

fwrite(match_performance, match_performance_path, sep = "\t", quote = FALSE, col.names = TRUE)
fwrite(apa_detect_performance, de_apa_performance_path, sep = "\t", quote = FALSE, col.names = TRUE)
fwrite(pas_quantify_performance, pas_quantify_performance_path, sep = "\t", quote = FALSE, col.names = TRUE)
fwrite(data.frame(value = pas_merged_dt[all_match == TRUE, unique(pas_name)], stringsAsFactors = FALSE), match_pas_path, sep = ",", quote = FALSE, col.names = TRUE)
fwrite(data.frame(value = te_pas_match[pd_pas_count > 1, unique(gt_name)], stringsAsFactors = FALSE), match_te_path, sep = ",", quote = FALSE, col.names = TRUE)