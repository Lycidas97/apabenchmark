#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))

# Create command line arguments
parser <- ArgumentParser(description='Integrated PAS to DaPars2 converter with bedtools-based gene mapping')

parser$add_argument('--pas_counts', type='character', required=TRUE,
                    help='PAS count matrix file (e.g., pas_counts.tsv)')
parser$add_argument('--pas_coords', type='character', required=TRUE,
                    help='PAS coordinate file in BED format (e.g., pas.bed)')
parser$add_argument('--pas_annotation', type='character', required=TRUE,
                    help='PAS annotation file with exon_id column (e.g., mm10_sim_pas_gn2000_rep1.bed)')
parser$add_argument('--celltype_info', type='character', required=TRUE,
                    help='Cell type annotation file with cell barcode and type')
parser$add_argument('-o', '--output', type='character', default='integrated_dapars2_output.txt',
                    help='Output file name')
parser$add_argument('--min_exp', type='double', default=0,
                    help='Minimum expression threshold')
parser$add_argument('--min_cells', type='integer', default=0,
                    help='Minimum number of cells with expression above threshold')
parser$add_argument('--debug', action='store_true', default=FALSE,
                    help='Enable debug output')

args <- parser$parse_args()

cat("=== Integrated PAS to DaPars2 Converter ===\n")
cat("PAS counts file:", args$pas_counts, "\n")
cat("PAS coords file:", args$pas_coords, "\n")
cat("PAS annotation file:", args$pas_annotation, "\n")
cat("Cell type info file:", args$celltype_info, "\n")
cat("Output file:", args$output, "\n")
cat("Min expression:", args$min_exp, "\n")
cat("Min cells:", args$min_cells, "\n")
cat("Debug mode:", args$debug, "\n")
cat("\n")

#===============================================
# Module 1: Extract Gene Regions from PAS Annotation
#===============================================
extract_gene_regions <- function(annotation_file, debug = FALSE) {
  cat("Step 1: Extracting gene regions from PAS annotation...\n")

  pas_annotation <- read.table(annotation_file, sep='\t', header=TRUE, stringsAsFactors=FALSE)

  if (debug) cat("Original annotation dimensions:", nrow(pas_annotation), "x", ncol(pas_annotation), "\n")

  # Parse exon_id column to extract gene region information
  gene_regions <- data.frame()

  for (i in 1:nrow(pas_annotation)) {
    exon_id <- pas_annotation$exon_id[i]

    if (is.na(exon_id) || exon_id == "exon_id") {
      next
    }

    # Parse format: ENSMUSG00000104328.1:Gm37323:chr1:4583128:4585585:-:TE
    parts <- strsplit(exon_id, ":", fixed=TRUE)[[1]]

    if (length(parts) >= 6) {
      gene_id <- parts[1]
      gene_symbol <- parts[2]
      chr <- parts[3]
      start <- as.numeric(parts[4])
      end <- as.numeric(parts[5])
      strand <- parts[6]

      gene_regions <- rbind(gene_regions, data.frame(
        chr = chr,
        start = start,
        end = end,
        gene_id = gene_id,
        gene_symbol = gene_symbol,
        strand = strand,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Remove duplicates
  gene_regions <- gene_regions[!duplicated(gene_regions), ]

  cat("Unique gene regions extracted:", nrow(gene_regions), "\n")

  return(gene_regions)
}

#===============================================
# Module 2: Create PAS-Gene Mapping using bedtools
#===============================================
create_pas_gene_mapping <- function(pas_coords_file, gene_regions, debug = FALSE) {
  cat("Step 2: Creating PAS-gene mapping with bedtools...\n")

  empty_mapping <- data.frame(
    pas_chr = character(),
    pas_start = numeric(),
    pas_end = numeric(),
    pas_id = character(),
    pas_position = numeric(),
    pas_strand = character(),
    gene_chr = character(),
    gene_start = numeric(),
    gene_end = numeric(),
    gene_id = character(),
    gene_symbol = character(),
    gene_strand = character(),
    stringsAsFactors = FALSE
  )

  if (nrow(gene_regions) == 0) {
    cat("No gene regions available; returning empty PAS-gene mapping.\n")
    return(empty_mapping)
  }

  # Write gene regions to temporary BED file
  temp_gene_bed <- tempfile(fileext = ".bed")
  write.table(gene_regions[, c("chr", "start", "end", "gene_id", "gene_symbol", "strand")],
              temp_gene_bed, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

  # Run bedtools intersect
  system_check <- system("which bedtools > /dev/null 2>&1", ignore.stderr = TRUE)

  if (system_check != 0) {
    stop("Error: bedtools not found. Please install bedtools or ensure it's in your PATH.")
  }

  temp_mapping_file <- tempfile(fileext = ".bed")
  cmd <- paste("bedtools intersect -s -a", pas_coords_file, "-b", temp_gene_bed, "-wa -wb >", temp_mapping_file)

  if (debug) cat("Running command:", cmd, "\n")

  result <- system(cmd, ignore.stderr = TRUE)

  if (result != 0) {
    stop("Error running bedtools. Please check your input files.")
  }

  mapping_size <- file.info(temp_mapping_file)$size
  if (is.na(mapping_size) || mapping_size == 0) {
    unlink(c(temp_gene_bed, temp_mapping_file))
    cat("No overlapping PAS-gene mappings found.\n")
    return(empty_mapping)
  }

  # Read mapping results
  pas_mapping <- read.table(temp_mapping_file, sep='\t', stringsAsFactors=FALSE, check.names=FALSE)
  colnames(pas_mapping) <- c("pas_chr", "pas_start", "pas_end", "pas_id", "pas_position", "pas_strand",
                             "gene_chr", "gene_start", "gene_end", "gene_id", "gene_symbol", "gene_strand")

  # Clean up temporary files
  unlink(c(temp_gene_bed, temp_mapping_file))

  cat("PAS-gene mapping created:", nrow(pas_mapping), "mappings\n")

  return(pas_mapping)
}

#===============================================
# Module 3: Read Input Data
#===============================================
read_input_data <- function(pas_counts_file, celltype_file, debug = FALSE) {
  cat("Step 3: Reading input data...\n")

  # Read cell type information
  celltype_info <- read.table(celltype_file, sep='\t', header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
  if (debug) cat("Cell type info dimensions:", nrow(celltype_info), "x", ncol(celltype_info), "\n")
  if (nrow(celltype_info) == 0) {
    stop("Cell type information file is empty.")
  }

  cell_barcodes <- as.character(celltype_info[, 1])
  cell_types_all <- as.character(celltype_info[, 2])
  unique_cell_types <- sort(unique(cell_types_all))
  cat("Cell types found (", length(unique_cell_types), "):", paste(unique_cell_types, collapse=", "), "\n")

  return_empty_input <- function(reason) {
    cat(reason, "\n")
    return(list(
      pas_counts = data.frame(),
      cell_types = character(0),
      unique_cell_types = unique_cell_types
    ))
  }

  count_file_size <- file.info(pas_counts_file)$size
  if (is.na(count_file_size) || count_file_size == 0) {
    return(return_empty_input("PAS count matrix is empty. Continuing with empty matrix."))
  }

  # Read PAS count matrix
  pas_counts <- tryCatch(
    read.table(pas_counts_file, sep='\t', header=TRUE, stringsAsFactors=FALSE, check.names=FALSE),
    error = function(e) {
      if (debug) {
        cat("Failed to parse PAS count matrix:", e$message, "\n")
      }
      NULL
    }
  )
  if (is.null(pas_counts)) {
    return(return_empty_input("PAS count matrix is malformed and cannot be parsed. Continuing with empty matrix."))
  }
  if (debug) cat("PAS count matrix dimensions:", nrow(pas_counts), "x", ncol(pas_counts), "\n")
  if (nrow(pas_counts) == 0 || ncol(pas_counts) == 0) {
    return(return_empty_input("PAS count matrix has no usable entries. Continuing with empty matrix."))
  }

  # Align PAS count matrix rows with cell type metadata by barcode
  rowname_match_rate <- mean(rownames(pas_counts) %in% cell_barcodes)
  barcode_col_idx <- NA_integer_

  if (is.finite(rowname_match_rate) && rowname_match_rate >= 0.5) {
    pas_barcodes <- rownames(pas_counts)
    barcode_match_rate <- rowname_match_rate
  } else {
    candidate_cols <- which(vapply(pas_counts, is.character, logical(1)))
    if (length(candidate_cols) == 0) {
      candidate_cols <- 1
    }

    overlap_scores <- vapply(candidate_cols, function(col_idx) {
      mean(as.character(pas_counts[[col_idx]]) %in% cell_barcodes)
    }, numeric(1))

    barcode_col_idx <- candidate_cols[which.max(overlap_scores)]
    barcode_match_rate <- max(overlap_scores)

    if (!is.finite(barcode_match_rate) || barcode_match_rate < 0.5) {
      return(return_empty_input(
        "Unable to detect barcode source in PAS count matrix with sufficient overlap. Continuing with empty matrix."
      ))
    }

    pas_barcodes <- as.character(pas_counts[[barcode_col_idx]])
  }

  # Remove duplicated barcodes while keeping first occurrence
  pas_keep <- !duplicated(pas_barcodes)
  pas_counts <- pas_counts[pas_keep, , drop = FALSE]
  pas_barcodes <- pas_barcodes[pas_keep]

  cell_keep <- !duplicated(cell_barcodes)
  celltype_info <- celltype_info[cell_keep, , drop = FALSE]
  cell_barcodes <- as.character(celltype_info[, 1])

  matched_idx <- match(cell_barcodes, pas_barcodes)
  valid <- !is.na(matched_idx)

  if (!any(valid)) {
    return(return_empty_input(
      "No overlapping cell barcodes found between PAS counts and cell type file. Continuing with empty matrix."
    ))
  }

  if (debug) {
    cat("Detected barcode column index:", ifelse(is.na(barcode_col_idx), "ROWNAMES", barcode_col_idx), "\n")
    cat("Barcode overlap rate:", round(barcode_match_rate, 4), "\n")
    cat("Matched cells:", sum(valid), "/", length(valid), "\n")
  }

  celltype_info <- celltype_info[valid, , drop = FALSE]
  pas_counts <- pas_counts[matched_idx[valid], , drop = FALSE]
  rownames(pas_counts) <- as.character(celltype_info[, 1])

  # Drop barcode column (if barcodes were stored as a normal column)
  if (!is.na(barcode_col_idx)) {
    pas_counts <- pas_counts[, -barcode_col_idx, drop = FALSE]
  }
  if (ncol(pas_counts) == 0) {
    return(return_empty_input(
      "No PAS columns remain after barcode alignment. Continuing with empty matrix."
    ))
  }

  # Extract aligned cell types
  cell_types <- as.character(celltype_info[, 2])

  return(list(pas_counts = pas_counts, cell_types = cell_types, unique_cell_types = unique_cell_types))
}

#===============================================
# Module 4: Create PAS Info Table
#===============================================
create_pas_info_table <- function(pas_mapping, pas_counts, debug = FALSE) {
  cat("Step 4: Creating PAS info table...\n")

  pas_info <- data.frame(
    pas_id = character(),
    gene_id = character(),
    gene_symbol = character(),
    gene_chr = character(),
    gene_start = numeric(),
    gene_end = numeric(),
    gene_strand = character(),
    pas_position = numeric(),
    pas_strand = character(),
    stringsAsFactors = FALSE
  )

  if (nrow(pas_mapping) == 0) {
    cat("No PAS mappings available; PAS info table is empty.\n")
    return(pas_info)
  }

  for (i in 1:nrow(pas_mapping)) {
    mapping_row <- pas_mapping[i, ]
    pas_id <- mapping_row$pas_id

    # Check if this PAS exists in the count matrix
    if (pas_id %in% colnames(pas_counts)) {
      # Create PAS info entry in DaPars2 format
      pas_info <- rbind(pas_info, data.frame(
        pas_id = pas_id,
        gene_id = sub("\\..*$", "", mapping_row$gene_id),  # Remove version from gene ID
        gene_symbol = mapping_row$gene_symbol,
        gene_chr = mapping_row$gene_chr,
        gene_start = mapping_row$gene_start,
        gene_end = mapping_row$gene_end,
        gene_strand = mapping_row$gene_strand,
        pas_position = mapping_row$pas_position,
        pas_strand = mapping_row$pas_strand,
        stringsAsFactors = FALSE
      ))
    }
  }

  cat("PAS found in count matrix:", nrow(pas_info), "\n")
  return(pas_info)
}

#===============================================
# Module 5: Process Genes and Calculate PDUI
#===============================================
process_genes_and_calculate_pdui <- function(pas_info, pas_counts, cell_types, unique_cell_types,
                                            min_exp = 0, min_cells = 0, debug = FALSE) {
  cat("Step 5: Processing genes with multiple PAS...\n")

  unique_genes <- unique(pas_info$gene_id)
  cat("Number of unique genes:", length(unique_genes), "\n")

  results <- data.frame()

  for (gene_id in unique_genes) {
    if (debug) cat("\nProcessing gene:", gene_id, "\n")

    gene_pas_info <- pas_info[pas_info$gene_id == gene_id, ]

    if (nrow(gene_pas_info) < 2) {
      if (debug) cat("  Skipping: Less than 2 PAS sites\n")
      next
    }

    # Get PAS columns for this gene
    pas_columns <- gene_pas_info$pas_id
    missing_cols <- setdiff(pas_columns, colnames(pas_counts))
    if (length(missing_cols) > 0) {
      if (debug) cat("  Missing columns:", paste(missing_cols, collapse=", "), "\n")
      next
    }

    # Get expression data for this gene
    gene_exp <- pas_counts[, pas_columns, drop=FALSE]
    gene_exp$total_exp <- rowSums(gene_exp, na.rm=TRUE)

    # Filter cells with sufficient expression
    sufficient_cells <- gene_exp$total_exp >= min_exp
    n_sufficient <- sum(sufficient_cells)

    if (debug) cat("  Cells with sufficient expression:", n_sufficient, "/", length(sufficient_cells), "\n")

    if (n_sufficient < min_cells) {
      if (debug) cat("  Skipping: Insufficient cells with expression\n")
      next
    }

    # Get PAS positions and determine upstream/downstream
    pas_positions <- gene_pas_info$pas_position
    pas_strands <- gene_pas_info$gene_strand[1]  # All PAS for a gene should have same strand

    # Determine downstream PAS
    if (pas_strands == "+") {
      # Positive strand: larger position = downstream
      downstream_idx <- which.max(pas_positions)
      upstream_idx <- which.min(pas_positions)
    } else {
      # Negative strand: smaller position = downstream
      downstream_idx <- which.min(pas_positions)
      upstream_idx <- which.max(pas_positions)
    }

    # Get columns for upstream and downstream PAS
    upstream_col <- gene_pas_info$pas_id[upstream_idx]
    downstream_col <- gene_pas_info$pas_id[downstream_idx]

    if (debug) {
      cat("  Upstream PAS:", upstream_col, "Position:", pas_positions[upstream_idx], "\n")
      cat("  Downstream PAS:", downstream_col, "Position:", pas_positions[downstream_idx], "\n")
    }

    # Initialize result row with basic columns in correct DaPars2 order
    result_row <- data.frame(
      Gene = paste0(gene_id, "|", gene_pas_info$gene_symbol[1], "|", gene_pas_info$gene_chr[1], "|", gene_pas_info$gene_strand[1]),
      fit_value = 0,
      Predicted_Proximal_APA = pas_positions[upstream_idx],
      Loci = paste0(gene_pas_info$gene_chr[1], ":", gene_pas_info$gene_start[1], "-", gene_pas_info$gene_end[1]),
      stringsAsFactors = FALSE
    )

    # Calculate expression for each group using celltype naming
    valid_groups <- 0
    for (i in 1:length(unique_cell_types)) {
      cell_type <- unique_cell_types[i]
      group_cells <- cell_types == cell_type
      group_sufficient <- sufficient_cells & group_cells

      # Clean cell type name for use in column names
      clean_celltype <- gsub("[^a-zA-Z0-9_]", "_", cell_type)

      if (sum(group_sufficient) >= min_cells) {
        valid_groups <- valid_groups + 1

        down_exp_values <- gene_exp[group_sufficient, downstream_col]
        up_exp_values <- gene_exp[group_sufficient, upstream_col]

        # Calculate sum of expression values
        group_down_exp <- sum(down_exp_values, na.rm=TRUE)
        group_up_exp <- sum(up_exp_values, na.rm=TRUE)

        # Calculate PDUI: downstream sum / (downstream sum + upstream sum)
        group_pdui <- group_down_exp / (group_down_exp + group_up_exp)

        if (debug) {
          cat("  Group", cell_type, "- PDUI:", round(group_pdui, 3),
              "Down:", round(group_down_exp, 2), "Up:", round(group_up_exp, 2), "\n")
        }

        # Add columns using celltype naming
        long_col_name <- paste0(clean_celltype, "_long_exp")
        short_col_name <- paste0(clean_celltype, "_short_exp")
        pdui_col_name <- paste0(clean_celltype, "_PDUI")

        result_row[[long_col_name]] <- group_down_exp
        result_row[[short_col_name]] <- group_up_exp
        result_row[[pdui_col_name]] <- group_pdui
      } else {
        # Add NA columns for groups with insufficient cells
        long_col_name <- paste0(clean_celltype, "_long_exp")
        short_col_name <- paste0(clean_celltype, "_short_exp")
        pdui_col_name <- paste0(clean_celltype, "_PDUI")

        result_row[[long_col_name]] <- NA
        result_row[[short_col_name]] <- NA
        result_row[[pdui_col_name]] <- NA

        if (debug) cat("  Group", cell_type, "- Insufficient cells\n")
      }
    }

    if (valid_groups > 0) {
      results <- rbind(results, result_row)
      if (debug) cat("  Added gene to results\n")
    }
  }

  return(results)
}

#===============================================
# Module 6: Save Results
#===============================================
build_empty_dapars2_results <- function(unique_cell_types) {
  empty_results <- data.frame(
    Gene = character(),
    fit_value = numeric(),
    Predicted_Proximal_APA = numeric(),
    Loci = character(),
    stringsAsFactors = FALSE
  )

  for (cell_type in unique_cell_types) {
    clean_celltype <- gsub("[^a-zA-Z0-9_]", "_", cell_type)
    empty_results[[paste0(clean_celltype, "_long_exp")]] <- numeric()
    empty_results[[paste0(clean_celltype, "_short_exp")]] <- numeric()
    empty_results[[paste0(clean_celltype, "_PDUI")]] <- numeric()
  }

  empty_results
}

save_results <- function(results, output_file, unique_cell_types, debug = FALSE) {
  cat("\nStep 6: Final results...\n")
  cat("Results dimensions:", nrow(results), "x", ncol(results), "\n")

  results_to_write <- results
  if (nrow(results) == 0 || ncol(results) == 0) {
    cat("No results generated; writing an empty DaPars2-format table.\n")
    results_to_write <- build_empty_dapars2_results(unique_cell_types)
  } else if (debug) {
    cat("Results preview:\n")
    print(head(results[, 1:min(6, ncol(results))]))
    if (ncol(results) > 6) cat("... (", ncol(results)-6, "more columns)\n")
  }

  # Save results
  write.table(results_to_write, output_file, sep='\t', quote=FALSE, row.names=FALSE, na='NA')
  cat("\nIntegrated conversion completed successfully!\n")
  cat("Output saved to:", output_file, "\n")
  cat("Number of genes processed:", nrow(results_to_write), "\n")
  cat("Number of groups:", length(unique_cell_types), "\n")
  cat("Groups:", paste(unique_cell_types, collapse=", "), "\n")

  # Show column names for verification
  cat("\nOutput columns:\n")
  for (i in 1:ncol(results_to_write)) {
    cat(i, ":", colnames(results_to_write)[i], "\n")
  }
}

#===============================================
# Main Execution Pipeline
#===============================================
cat("Starting integrated analysis pipeline...\n\n")

tryCatch({
  # Step 1: Extract gene regions
  gene_regions <- extract_gene_regions(args$pas_annotation, args$debug)

  # Step 2: Create PAS-gene mapping
  pas_mapping <- create_pas_gene_mapping(args$pas_coords, gene_regions, args$debug)

  # Step 3: Read input data
  input_data <- read_input_data(args$pas_counts, args$celltype_info, args$debug)

  # Step 4: Create PAS info table
  pas_info <- create_pas_info_table(pas_mapping, input_data$pas_counts, args$debug)

  # Step 5: Process genes and calculate PDUI
  results <- process_genes_and_calculate_pdui(
    pas_info, input_data$pas_counts, input_data$cell_types,
    input_data$unique_cell_types, args$min_exp, args$min_cells, args$debug
  )

  # Step 6: Save results
  save_results(results, args$output, input_data$unique_cell_types, args$debug)

}, error = function(e) {
  cat("Error occurred during execution:\n")
  cat(e$message, "\n")
  quit(status = 1)
})

cat("\nPipeline completed successfully!\n")
