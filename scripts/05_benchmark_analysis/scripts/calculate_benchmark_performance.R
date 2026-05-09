library(optparse)

option_list <- list(
  make_option(c("--tool"), dest = "tool", type = "character"),
  make_option(c("--sample"), dest = "sample", type = "character"),
  make_option(c("--gt_pas"), dest = "gt_pas", type = "character"),
  make_option(c("--pd_pas"), dest = "pd_pas", type = "character"),
  make_option(c("--gt_mtx"), dest = "gt_mtx", type = "character"),
  make_option(c("--pd_mtx"), dest = "pd_mtx", type = "character"),
  make_option(c("--output_dir"), dest = "output_dir", type = "character"),
  make_option(c("--dexseq_split_num"), dest = "dexseq_split_num", type = "integer")
)

opt <- parse_args(OptionParser(option_list = option_list))

required_args <- c("tool", "sample", "gt_pas", "pd_pas", "gt_mtx", "pd_mtx", "output_dir", "dexseq_split_num")
for (arg in required_args) {
  if (is.null(opt[[arg]]) || opt[[arg]] == "" || (arg == "dexseq_split_num" && is.na(opt[[arg]]))) {
    stop(sprintf("Missing required argument: --%s", arg))
  }
}

if (opt$dexseq_split_num < 1) {
  stop("--dexseq_split_num must be >= 1")
}

input_paths <- c(opt$gt_pas, opt$pd_pas, opt$gt_mtx, opt$pd_mtx)
for (input_path in input_paths) {
  if (!file.exists(input_path)) {
    stop(sprintf("Input file not found: %s", input_path))
  }
}

script_arg <- commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]
script_path <- sub("--file=", "", script_arg[1])
script_dir <- dirname(normalizePath(script_path))

stage10_script <- file.path(script_dir, "10_stage_pas_match_quantify.R")
stage20_script <- file.path(script_dir, "20_stage_apa_and_te.R")
stage30_script <- file.path(script_dir, "30_stage_export_results.R")

required_scripts <- c(stage10_script, stage20_script, stage30_script)
for (stage_script in required_scripts) {
  if (!file.exists(stage_script)) {
    stop(sprintf("Stage script not found: %s", stage_script))
  }
}

rscript_bin <- file.path(R.home("bin"), "Rscript")
if (!file.exists(rscript_bin)) {
  rscript_bin <- "Rscript"
}

run_stage <- function(stage_script, stage_args, stage_name) {
  status <- system2(rscript_bin, args = c(stage_script, stage_args), stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop(sprintf("%s failed with exit code %d", stage_name, as.integer(status)))
  }
}

run_stage(
  stage_script = stage10_script,
  stage_args = c(
    "--tool", opt$tool,
    "--sample", opt$sample,
    "--gt_pas", opt$gt_pas,
    "--pd_pas", opt$pd_pas,
    "--gt_mtx", opt$gt_mtx,
    "--pd_mtx", opt$pd_mtx,
    "--output_dir", opt$output_dir,
    "--dexseq_split_num", as.character(opt$dexseq_split_num)
  ),
  stage_name = "stage10"
)

stage10_bundle <- file.path(opt$output_dir, opt$tool, paste0(opt$sample, "_stage10_bundle.rds"))
if (!file.exists(stage10_bundle)) {
  stop(sprintf("stage10 bundle not found after stage10: %s", stage10_bundle))
}

run_stage(
  stage_script = stage20_script,
  stage_args = c("--stage10_bundle", stage10_bundle),
  stage_name = "stage20"
)

stage20_bundle <- file.path(opt$output_dir, opt$tool, paste0(opt$sample, "_stage20_bundle.rds"))
if (!file.exists(stage20_bundle)) {
  stop(sprintf("stage20 bundle not found after stage20: %s", stage20_bundle))
}

run_stage(
  stage_script = stage30_script,
  stage_args = c("--stage20_bundle", stage20_bundle),
  stage_name = "stage30"
)

cat(sprintf("Done: unified wrapper finished for %s / %s\n", opt$tool, opt$sample))
