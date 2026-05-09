library(optparse)
library(data.table)

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

option_list <- list(
  make_option(c("--stage20_bundle"), dest = "stage20_bundle", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$stage20_bundle) || opt$stage20_bundle == "") {
  stop("Missing required argument: --stage20_bundle")
}
if (!file.exists(opt$stage20_bundle)) {
  stop(sprintf("stage20 bundle not found: %s", opt$stage20_bundle))
}

bundle <- readRDS(opt$stage20_bundle)
cfg <- bundle$cfg
sample <- bundle$sample
tool <- bundle$tool
protocol <- bundle$protocol

match_performance <- bundle$match_performance
pas_quantify_performance <- bundle$pas_quantify_performance
apa_detect_performance <- bundle$apa_detect_performance
quality_report <- bundle$quality_issues %||% data.table(
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
  quality_report <<- rbind(
    quality_report,
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

window_sizes_half <- cfg$window_sizes / 2
p_stats <- c(paste0("wilcox_", cfg$index_types), "fisher", "dexseq")

full_filter_types <- expand.grid(
  filter_type_1 = paste0(rep(p_stats, each = length(cfg$p_thresholds)), "_", rep(cfg$p_thresholds, length(p_stats))),
  filter_type_2 = c(
    paste0(rep(cfg$diff_types, length(cfg$index_diff_thresholds)), "_", rep(cfg$index_diff_thresholds, each = length(cfg$diff_types))),
    paste0("dexseq_log2fc_", cfg$log2fc_thresholds)
  ),
  stringsAsFactors = FALSE
)

window_based_meta <- expand.grid(
  sample = sample,
  tool = tool,
  protocol = protocol,
  match_type = paste0(rep(c("pd", "gt"), each = length(window_sizes_half)), "_", window_sizes_half),
  stringsAsFactors = FALSE
)

additional_te_block <- data.frame(
  sample = sample,
  tool = tool,
  protocol = protocol,
  match_type = "te",
  stringsAsFactors = FALSE
)

meta_combos <- rbind(window_based_meta, additional_te_block)
apa_complete_meta <- merge(meta_combos, full_filter_types, by = character(), allow.cartesian = TRUE)
apa_detect_performance_full <- merge(
  apa_complete_meta,
  apa_detect_performance,
  by = c("sample", "tool", "protocol", "match_type", "filter_type_1", "filter_type_2"),
  all.x = TRUE
)
setDT(apa_detect_performance_full)
required_metric_cols <- c("precision", "recall", "f1", "gt_count", "pd_count")
for (metric_col in required_metric_cols) {
  if (!(metric_col %in% names(apa_detect_performance_full))) {
    apa_detect_performance_full[, (metric_col) := 0]
  }
}
apa_detect_performance_full[is.na(apa_detect_performance_full)] <- 0

gt_nonzero_rows <- apa_detect_performance_full[gt_count > 0, .N]
pd_nonzero_rows <- apa_detect_performance_full[pd_count > 0, .N]
if (nrow(apa_detect_performance_full) > 0 && gt_nonzero_rows == 0) {
  append_quality_issue(
    stage = "stage30",
    issue_code = "all_gt_count_zero",
    severity = "warning",
    detail = "All rows in final de_apa_performance have gt_count=0.",
    metric_name = "rows_total",
    metric_value = nrow(apa_detect_performance_full)
  )
}
if (nrow(apa_detect_performance_full) > 0 && pd_nonzero_rows == 0) {
  append_quality_issue(
    stage = "stage30",
    issue_code = "all_pd_count_zero",
    severity = "warning",
    detail = "All rows in final de_apa_performance have pd_count=0.",
    metric_name = "rows_total",
    metric_value = nrow(apa_detect_performance_full)
  )
}
quality_report <- unique(quality_report)

output_dir <- bundle$output_dir
tool_dir <- file.path(output_dir, tool)
if (!dir.exists(tool_dir)) {
  dir.create(tool_dir, recursive = TRUE)
}

output_path_prefix <- file.path(tool_dir, sample)
fwrite(match_performance, paste0(output_path_prefix, "_match_performance.tsv"), sep = "\t", quote = FALSE)
fwrite(apa_detect_performance_full, paste0(output_path_prefix, "_de_apa_performance.tsv"), sep = "\t", quote = FALSE)
fwrite(pas_quantify_performance, paste0(output_path_prefix, "_pas_quantify_performance.tsv"), sep = "\t", quote = FALSE)
fwrite(quality_report, paste0(output_path_prefix, "_quality_report.tsv"), sep = "\t", quote = FALSE)

cat(sprintf("final outputs saved with prefix: %s\n", output_path_prefix))
