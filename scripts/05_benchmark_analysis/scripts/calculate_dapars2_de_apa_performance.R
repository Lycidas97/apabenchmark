library(optparse)
library(data.table)

safe_divide <- function(numerator, denominator) {
  if (is.na(numerator) || is.na(denominator) || denominator <= 0) {
    return(0)
  }
  numerator / denominator
}

safe_f1 <- function(precision, recall) {
  safe_divide(2 * precision * recall, precision + recall)
}

extract_gene_column <- function(dt, path) {
  candidates <- c("Genes", "Gene", "gene", "genes")
  matched <- candidates[candidates %in% colnames(dt)]
  if (length(matched) > 0) {
    return(matched[[1]])
  }
  stop(sprintf("Cannot find gene column in %s (expected one of: %s)", path, paste(candidates, collapse = ", ")))
}

extract_gene_id <- function(value) {
  parts <- strsplit(as.character(value), "|", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) {
    return("")
  }

  # Prefer gene-level Ensembl IDs when transcript and gene IDs coexist.
  gene_like <- parts[grepl("^ENS.*G[0-9]+", parts)]
  if (length(gene_like) > 0) {
    picked <- gene_like[[1]]
  } else {
    picked <- parts[[1]]
  }
  sub("\\.\\d+$", "", picked)
}

normalize_gene_ids <- function(values) {
  values <- as.character(values)
  values <- trimws(values)
  values <- values[!is.na(values) & nzchar(values)]
  ids <- vapply(values, extract_gene_id, FUN.VALUE = character(1), USE.NAMES = FALSE)
  unique(ids[nzchar(ids)])
}

read_gene_set <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Input file not found: %s", path))
  }
  dt <- fread(path)
  if (nrow(dt) == 0) {
    return(character(0))
  }
  gene_col <- extract_gene_column(dt, path)
  normalize_gene_ids(dt[[gene_col]])
}

option_list <- list(
  make_option(c("--pd_apa"), type = "character"),
  make_option(c("--gt_apa"), type = "character"),
  make_option(c("--tool"), type = "character"),
  make_option(c("--sample"), type = "character"),
  make_option(c("--match_type"), type = "character"),
  make_option(c("--filter_type_1"), type = "character"),
  make_option(c("--filter_type_2"), type = "character"),
  make_option(c("--output"), type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))
required_args <- c("pd_apa", "gt_apa", "tool", "sample", "match_type", "filter_type_1", "filter_type_2", "output")
for (arg in required_args) {
  if (is.null(opt[[arg]]) || opt[[arg]] == "") {
    stop(sprintf("Missing required argument: --%s", arg))
  }
}

pd_set <- read_gene_set(opt$pd_apa)
gt_set <- read_gene_set(opt$gt_apa)

tp <- length(intersect(pd_set, gt_set))
precision <- safe_divide(tp, length(pd_set))
recall <- safe_divide(tp, length(gt_set))
f1 <- safe_f1(precision, recall)

protocol <- strsplit(opt$sample, "_", fixed = TRUE)[[1]][1]

result <- data.table(
  sample = opt$sample,
  tool = opt$tool,
  protocol = protocol,
  match_type = opt$match_type,
  filter_type_1 = opt$filter_type_1,
  filter_type_2 = opt$filter_type_2,
  precision = precision,
  recall = recall,
  f1 = f1,
  gt_count = length(gt_set),
  pd_count = length(pd_set)
)

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
fwrite(result, opt$output, sep = "\t", quote = FALSE)

cat(sprintf("Saved dapars2 DE APA performance: %s\n", opt$output))
