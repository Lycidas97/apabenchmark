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
  make_option(c("--stage10_bundle"), dest = "stage10_bundle", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$stage10_bundle) || opt$stage10_bundle == "") {
  stop("Missing required argument: --stage10_bundle")
}
if (!file.exists(opt$stage10_bundle)) {
  stop(sprintf("stage10 bundle not found: %s", opt$stage10_bundle))
}

bundle <- readRDS(opt$stage10_bundle)
cfg <- bundle$cfg

sample <- bundle$sample
tool <- bundle$tool
protocol <- bundle$protocol
bedtools_bin <- bundle$bedtools_bin
gt_te_dt <- bundle$gt_te_dt
pd_pas_dt <- bundle$pd_pas_dt
pd_mtx <- bundle$pd_mtx
pd_col_num <- bundle$pd_col_num
grouped_gt_mtx <- bundle$grouped_gt_mtx
cell_type_map <- bundle$cell_type_map
match_performance <- bundle$match_performance
pas_quantify_performance <- bundle$pas_quantify_performance
apa_inputs <- bundle$apa_inputs
quality_issues <- bundle$quality_issues %||% data.table(
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

apa_detect_performance <- data.table()
fisher_qc <- data.table()
empty_apa_row <- function(match_type_value) {
  data.table(
    sample = sample,
    tool = tool,
    protocol = protocol,
    match_type = match_type_value,
    filter_type_1 = NA_character_,
    filter_type_2 = NA_character_,
    precision = 0,
    recall = 0,
    f1 = 0,
    gt_count = 0,
    pd_count = 0
  )
}
append_fisher_qc <- function(grouped_mtx, source_label, match_type_label) {
  qc <- summarize_fisher_input_validity(
    grouped_mtx = grouped_mtx,
    cell_type_map = cell_type_map,
    source_label = source_label,
    match_type_label = match_type_label
  )
  qc[, `:=`(
    sample = sample,
    tool = tool,
    protocol = protocol
  )]
  fisher_qc <<- rbind(fisher_qc, qc, use.names = TRUE, fill = TRUE)
  if (nrow(qc) > 0) {
    for (i in seq_len(nrow(qc))) {
      row_i <- qc[i]
      match_type_value <- as.character(row_i$match_type)
      source_value <- as.character(row_i$source)
      if (row_i$total_te > 0 && row_i$tested_te_ncol_gt1 == 0) {
        append_quality_issue(
          stage = "stage20",
          issue_code = "fisher_no_testable_te",
          severity = "warning",
          detail = sprintf("No TE with >1 PAS for Fisher test (source=%s, match_type=%s).", source_value, match_type_value),
          metric_name = "match_type",
          metric_value = match_type_value
        )
      }
      if (row_i$fisher_skipped_missing_cell_type > 0) {
        append_quality_issue(
          stage = "stage20",
          issue_code = "fisher_missing_cell_type_mapping",
          severity = "warning",
          detail = sprintf(
            "Fisher skipped %d TE due to missing cell_type mapping (source=%s, match_type=%s).",
            as.integer(row_i$fisher_skipped_missing_cell_type), source_value, match_type_value
          ),
          metric_name = "match_type",
          metric_value = match_type_value
        )
      }
      if (row_i$fisher_skipped_invalid_contingency > 0) {
        append_quality_issue(
          stage = "stage20",
          issue_code = "fisher_invalid_contingency",
          severity = "info",
          detail = sprintf(
            "Fisher skipped %d TE due to invalid contingency tables (source=%s, match_type=%s).",
            as.integer(row_i$fisher_skipped_invalid_contingency), source_value, match_type_value
          ),
          metric_name = "match_type",
          metric_value = match_type_value
        )
      }
      if (row_i$tested_te_ncol_gt1 > 0 && row_i$fisher_valid == 0) {
        append_quality_issue(
          stage = "stage20",
          issue_code = "fisher_no_valid_tests",
          severity = "warning",
          detail = sprintf("No valid Fisher tests after filtering (source=%s, match_type=%s).", source_value, match_type_value),
          metric_name = "match_type",
          metric_value = match_type_value
        )
      }
    }
  }
}

append_fisher_qc(grouped_gt_mtx, "gt", "gt_precompute")

gt_apa_cache <- precompute_gt_apa_statistics(
  grouped_gt_mtx = grouped_gt_mtx,
  cell_type_map = cell_type_map,
  index_types = cfg$index_types,
  diff_types = cfg$diff_types,
  p_thresholds = cfg$p_thresholds,
  log2fc_thresholds = cfg$log2fc_thresholds,
  index_diff_thresholds = cfg$index_diff_thresholds,
  dexseq_split_num = cfg$dexseq_split_num,
  adjustment_method = cfg$adjustment_method,
  fisher_simulate_p_value = cfg$fisher_simulate_p_value,
  fisher_B = cfg$fisher_B,
  fisher_seed = cfg$fisher_seed
)

for (match_type in names(apa_inputs)) {
  grouped_pd <- apa_inputs[[match_type]]
  append_fisher_qc(grouped_pd, "pd", match_type)
  if (length(grouped_pd) == 0) {
    append_quality_issue(
      stage = "stage20",
      issue_code = "apa_input_empty",
      severity = "warning",
      detail = sprintf("No PD APA input TE for match_type=%s; using all-zero APA row.", match_type),
      metric_name = "match_type",
      metric_value = match_type
    )
    current <- empty_apa_row(match_type)
  } else {
    current <- perform_full_apa_analysis(
      grouped_gt_mtx = grouped_gt_mtx,
      grouped_pd_mtx = grouped_pd,
      cell_type_map = cell_type_map,
      sample = sample,
      tool = tool,
      protocol = protocol,
      index_types = cfg$index_types,
      diff_types = cfg$diff_types,
      p_thresholds = cfg$p_thresholds,
      log2fc_thresholds = cfg$log2fc_thresholds,
      index_diff_thresholds = cfg$index_diff_thresholds,
      dexseq_split_num = cfg$dexseq_split_num,
      adjustment_method = cfg$adjustment_method,
      fisher_simulate_p_value = cfg$fisher_simulate_p_value,
      fisher_B = cfg$fisher_B,
      fisher_seed = cfg$fisher_seed,
      precomputed_gt = gt_apa_cache
    )
    current$match_type <- match_type
  }
  apa_detect_performance <- rbind(apa_detect_performance, current, use.names = TRUE, fill = TRUE)
}

gt_te_bed_tmp <- tempfile(fileext = ".bed")
predict_pas_tmp <- tempfile(fileext = ".bed")
fwrite(gt_te_dt, gt_te_bed_tmp, sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(pd_pas_dt, predict_pas_tmp, sep = "\t", quote = FALSE, col.names = FALSE)
output_file <- tempfile()
intersect_stderr <- tempfile()
intersect_status <- system2(
  bedtools_bin,
  args = c("intersect", "-s", "-a", gt_te_bed_tmp, "-b", predict_pas_tmp, "-wa", "-wb"),
  stdout = output_file,
  stderr = intersect_stderr
)
if (!identical(intersect_status, 0L)) {
  err_msg <- if (file.exists(intersect_stderr)) paste(readLines(intersect_stderr, warn = FALSE), collapse = "\n") else ""
  stop(sprintf("bedtools intersect failed: %s", err_msg))
}
unlink(c(intersect_stderr))
te_pas_match <- read_table_if_exists(output_file)
unlink(c(gt_te_bed_tmp, predict_pas_tmp, output_file))

if (nrow(te_pas_match) > 0) {
  rename_matched_columns(te_pas_match, 6, pd_col_num)
  te_pas_match[, pd_pas_count := .N, by = gt_name]
  setnames(te_pas_match, paste0("res_", pd_col_num), "pd_pas_name")
}

if (nrow(te_pas_match) > 0) {
  te_match_tp <- te_pas_match[(gt_score > 1) & (pd_pas_count > 1), length(unique(gt_name))]
  te_match_fp <- te_pas_match[(gt_score == 1) & (pd_pas_count > 1), length(unique(gt_name))]
  pd_te_pas_map <- split(te_pas_match[(!is.na(res_name)) & (pd_pas_count > 1), pd_pas_name], te_pas_match[(!is.na(res_name)) & (pd_pas_count > 1), gt_name])
} else {
  append_quality_issue(
    stage = "stage20",
    issue_code = "te_intersect_empty",
    severity = "warning",
    detail = "No strand-aware TE/PAS intersections for TE-level evaluation."
  )
  te_match_tp <- 0
  te_match_fp <- 0
  pd_te_pas_map <- list()
}

te_match_fn <- gt_te_dt[pas_num > 1, length(unique(te_name))] - te_match_tp
te_match_precision <- safe_divide(te_match_tp, te_match_tp + te_match_fp)
te_match_recall <- safe_divide(te_match_tp, te_match_tp + te_match_fn)
te_match_f1 <- safe_f1(te_match_precision, te_match_recall)

match_performance <- rbind(match_performance, data.table(sample = sample, tool = tool, protocol = protocol, match_type = "te", tp = te_match_tp, fp = te_match_fp, fn = te_match_fn, precision = te_match_precision, recall = te_match_recall, f1 = te_match_f1))

grouped_pd_mtx_mapped_by_te <- lapply(pd_te_pas_map, function(cols) {
  pd_mtx[, cols, drop = FALSE]
})
if (length(grouped_pd_mtx_mapped_by_te) > 0) {
  keep_te <- vapply(grouped_pd_mtx_mapped_by_te, function(x) !is.null(x) && ncol(x) > 0, logical(1))
  grouped_pd_mtx_mapped_by_te <- grouped_pd_mtx_mapped_by_te[keep_te]
}
append_fisher_qc(grouped_pd_mtx_mapped_by_te, "pd", "te")

if (length(grouped_pd_mtx_mapped_by_te) == 0) {
  append_quality_issue(
    stage = "stage20",
    issue_code = "te_apa_input_empty",
    severity = "warning",
    detail = "No PD TE-level APA inputs after TE/PAS mapping; using all-zero TE APA row.",
    metric_name = "match_type",
    metric_value = "te"
  )
  apa_te <- empty_apa_row("te")
} else {
  apa_te <- perform_full_apa_analysis(
    grouped_gt_mtx = grouped_gt_mtx,
    grouped_pd_mtx = grouped_pd_mtx_mapped_by_te,
    cell_type_map = cell_type_map,
    sample = sample,
    tool = tool,
    protocol = protocol,
    index_types = cfg$index_types,
    diff_types = cfg$diff_types,
    p_thresholds = cfg$p_thresholds,
    log2fc_thresholds = cfg$log2fc_thresholds,
    index_diff_thresholds = cfg$index_diff_thresholds,
    dexseq_split_num = cfg$dexseq_split_num,
    adjustment_method = cfg$adjustment_method,
    fisher_simulate_p_value = cfg$fisher_simulate_p_value,
    fisher_B = cfg$fisher_B,
    fisher_seed = cfg$fisher_seed,
    precomputed_gt = gt_apa_cache
  )
  apa_te$match_type <- "te"
}

apa_detect_performance <- rbind(apa_detect_performance, apa_te, use.names = TRUE, fill = TRUE)

stage20_bundle <- list(
  sample = sample,
  tool = tool,
  protocol = protocol,
  output_dir = bundle$output_dir,
  cfg = cfg,
  window_sizes = cfg$window_sizes,
  match_performance = match_performance,
  pas_quantify_performance = pas_quantify_performance,
  apa_detect_performance = apa_detect_performance,
  fisher_qc = fisher_qc,
  quality_issues = quality_issues
)

prefix <- sub("_stage10_bundle.rds$", "", opt$stage10_bundle)
stage20_path <- paste0(prefix, "_stage20_bundle.rds")
saveRDS(stage20_bundle, stage20_path)
fwrite(match_performance, paste0(prefix, "_stage20_match_performance.tsv"), sep = "\t", quote = FALSE)
fwrite(apa_detect_performance, paste0(prefix, "_stage20_de_apa_raw.tsv"), sep = "\t", quote = FALSE)
fwrite(fisher_qc, paste0(prefix, "_stage20_fisher_qc.tsv"), sep = "\t", quote = FALSE)
fwrite(quality_issues, paste0(prefix, "_stage20_quality_issues.tsv"), sep = "\t", quote = FALSE)

cat(sprintf("stage20 bundle saved: %s\n", stage20_path))
