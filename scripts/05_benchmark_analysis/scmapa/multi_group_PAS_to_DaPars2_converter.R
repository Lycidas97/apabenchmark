#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))

# Create command line arguments
parser <- ArgumentParser(description='Convert PAS expression matrix to scMAPA DaPars2 format supporting multiple groups')

parser$add_argument('-i', '--input', type='character', required=TRUE,
                    help='Input PAS expression matrix file (e.g., from DaPars2)')
parser$add_argument('-o', '--output', type='character', default='multi_group_output.txt',
                    help='Output file name')
parser$add_argument('--min_exp', type='double', default=0,
                    help='Minimum expression threshold')
parser$add_argument('--min_cells', type='integer', default=0,
                    help='Minimum number of cells with expression above threshold (default: 0)')
parser$add_argument('--debug', action='store_true', default=FALSE,
                    help='Enable debug output')

args <- parser$parse_args()

cat("=== Multi-Group PAS to DaPars2 Converter ===\n")
cat("Input file:", args$input, "\n")
cat("Output file:", args$output, "\n")
cat("Min expression:", args$min_exp, "\n")
cat("Min cells:", args$min_cells, "\n")
cat("Debug mode:", args$debug, "\n")
cat("\n")

# Read data
cat("Step 1: Reading data...\n")
data <- read.table(args$input, sep='\t', header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
cat("Original data dimensions:", nrow(data), "x", ncol(data), "\n")

# Check if there's a cell_type column
if ("cell_type" %in% colnames(data)) {
  cat("Found cell_type column\n")
  cell_types <- sort(unique(data$cell_type))
  cat("Cell types found (", length(cell_types), "):", paste(cell_types, collapse=", "), "\n")
} else {
  cat("Error: No cell_type column found. Input file must have a 'cell_type' column specifying group assignments.\n")
  quit(status=1)
}

# Extract PAS column information
pas_columns <- colnames(data)
pas_columns <- pas_columns[pas_columns != "cell_type"]
cat("Number of PAS columns:", length(pas_columns), "\n")

# Remove empty column names (caused by leading tab in file)
pas_columns <- pas_columns[pas_columns != "" & !is.na(pas_columns)]
cat("Number of PAS columns after removing empty names:", length(pas_columns), "\n")

# Parse PAS column names
cat("Step 2: Parsing PAS column names...\n")
pas_info_list <- list()

for (i in 1:length(pas_columns)) {
  col_name <- pas_columns[i]

  # Skip empty column names
  if (col_name == "" || is.na(col_name)) {
    next
  }

  # Parse: GeneID:GeneSymbol:Chr:Start:End:Strand:TE|Chr:Position:Position:Strand
  parts <- strsplit(col_name, "|", fixed=TRUE)[[1]]

  if (length(parts) == 2) {
    gene_part <- parts[1]
    pas_part <- parts[2]

    # Parse gene info
    gene_parts <- strsplit(gene_part, ":", fixed=TRUE)[[1]]
    if (length(gene_parts) >= 6) {
      gene_id <- gene_parts[1]
      gene_symbol <- gene_parts[2]
      chr <- gene_parts[3]
      start <- gene_parts[4]
      end <- gene_parts[5]
      strand <- gene_parts[6]

      # Parse PAS position
      pas_parts <- strsplit(pas_part, ":", fixed=TRUE)[[1]]
      if (length(pas_parts) >= 4) {
        pas_chr <- pas_parts[1]
        pas_position <- as.numeric(pas_parts[2])
        pas_strand <- pas_parts[4]

        # Clean gene ID to remove version number
        gene_id_clean <- sub("\\..*$", "", gene_id)

        pas_info_list[[length(pas_info_list) + 1]] <- list(
          original_column = col_name,
          gene_id = gene_id_clean,
          gene_symbol = gene_symbol,
          chr = chr,
          start = as.numeric(start),
          end = as.numeric(end),
          strand = strand,
          pas_position = pas_position,
          pas_strand = pas_strand
        )
      }
    }
  }
}

cat("Successfully parsed", length(pas_info_list), "PAS columns\n")

if (args$debug && length(pas_info_list) > 0) {
  cat("First 3 parsed PAS entries:\n")
  for (i in 1:min(3, length(pas_info_list))) {
    pas_info <- pas_info_list[[i]]
    cat("Column", i, ":", pas_info$gene_id, "|", pas_info$gene_symbol,
        "|", pas_info$chr, "|", pas_info$start, "|", pas_info$end, "|", pas_info$strand,
        "| Position:", pas_info$pas_position, "\n")
  }
  cat("\n")
}

# Convert to data frame
cat("Step 3: Converting to data frame...\n")
if (length(pas_info_list) > 0) {
  pas_info_df <- do.call(rbind, lapply(pas_info_list, function(x) as.data.frame(x)))
  cat("PAS info dataframe dimensions:", nrow(pas_info_df), "x", ncol(pas_info_df), "\n")

  # Group by gene and analyze
  cat("Step 4: Analyzing PAS by gene...\n")
  unique_genes <- unique(pas_info_df$gene_id)
  cat("Number of unique genes:", length(unique_genes), "\n")

  results <- data.frame()

  for (gene_id in unique_genes) {
    if (args$debug) cat("\nProcessing gene:", gene_id, "\n")

    gene_pas <- pas_info_df[pas_info_df$gene_id == gene_id, ]

    if (nrow(gene_pas) < 2) {
      if (args$debug) cat("  Skipping: Less than 2 PAS sites\n")
      next
    }

    # Get expression data for this gene
    gene_columns <- gene_pas$original_column
    missing_cols <- setdiff(gene_columns, colnames(data))
    if (length(missing_cols) > 0) {
      if (args$debug) cat("  Missing columns:", paste(missing_cols, collapse=", "), "\n")
      next
    }

    gene_exp <- data[, gene_columns, drop=FALSE]
    gene_exp$total_exp <- rowSums(gene_exp, na.rm=TRUE)

    # Filter cells with sufficient expression
    sufficient_cells <- gene_exp$total_exp >= args$min_exp
    n_sufficient <- sum(sufficient_cells)

    if (args$debug) cat("  Cells with sufficient expression:", n_sufficient, "/", length(sufficient_cells), "\n")

    if (n_sufficient < args$min_cells) {
      if (args$debug) cat("  Skipping: Insufficient cells with expression\n")
      next
    }

    # Get PAS positions and determine upstream/downstream
    pas_positions <- gene_pas$pas_position
    pas_strands <- gene_pas$strand[1]  # All PAS for a gene should have same strand

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
    upstream_col <- gene_pas$original_column[upstream_idx]
    downstream_col <- gene_pas$original_column[downstream_idx]

    if (args$debug) {
      cat("  Upstream PAS index:", upstream_idx, "Position:", pas_positions[upstream_idx], "\n")
      cat("  Downstream PAS index:", downstream_idx, "Position:", pas_positions[downstream_idx], "\n")
    }

    # Initialize result row with basic columns in correct DaPars2 order
    result_row <- data.frame(
      Gene = paste0(gene_id, "|", gene_pas$gene_symbol[1], "|", gene_pas$chr[1], "|", gene_pas$strand[1]),
      fit_value = 0,
      Predicted_Proximal_APA = pas_positions[upstream_idx],
      Loci = paste0(gene_pas$chr[1], ":", gene_pas$start[1], "-", gene_pas$end[1]),
      stringsAsFactors = FALSE
    )

    # Calculate expression for each group using standard DaPars2 naming
    valid_groups <- 0
    for (i in 1:length(cell_types)) {
      cell_type <- cell_types[i]
      group_cells <- data$cell_type == cell_type
      group_sufficient <- sufficient_cells & group_cells

      # Use standard DaPars2 naming: Sample_1, Sample_2, Sample_3...
      # sample_name <- paste0("Sample_", i)
      sample_name <- gsub("[^a-zA-Z0-9_]", "_", cell_type)

      if (sum(group_sufficient) >= args$min_cells) {
        valid_groups <- valid_groups + 1

        down_exp_values <- gene_exp[group_sufficient, downstream_col]
        up_exp_values <- gene_exp[group_sufficient, upstream_col]

        # Calculate sum of expression values
        group_down_exp <- sum(down_exp_values, na.rm=TRUE)
        group_up_exp <- sum(up_exp_values, na.rm=TRUE)

        # Calculate PDUI: downstream sum / (downstream sum + upstream sum)
        group_pdui <- group_down_exp / (group_down_exp + group_up_exp)

        if (args$debug) {
          cat("  Group", cell_type, "(", sample_name, ") - PDUI:", round(group_pdui, 3),
              "Down:", round(group_down_exp, 2), "Up:", round(group_up_exp, 2), "\n")
        }

        # Add columns using standard DaPars2 naming
        long_col_name <- paste0(sample_name, "_long_exp")
        short_col_name <- paste0(sample_name, "_short_exp")
        pdui_col_name <- paste0(sample_name, "_PDUI")

        result_row[[long_col_name]] <- group_down_exp
        result_row[[short_col_name]] <- group_up_exp
        result_row[[pdui_col_name]] <- group_pdui
      } else {
        # Add NA columns for groups with insufficient cells
        long_col_name <- paste0(sample_name, "_long_exp")
        short_col_name <- paste0(sample_name, "_short_exp")
        pdui_col_name <- paste0(sample_name, "_PDUI")

        result_row[[long_col_name]] <- NA
        result_row[[short_col_name]] <- NA
        result_row[[pdui_col_name]] <- NA

        if (args$debug) cat("  Group", cell_type, "(", sample_name, ") - Insufficient cells\n")
      }
    }

    if (valid_groups > 0) {
      results <- rbind(results, result_row)
      if (args$debug) cat("  Added gene to results\n")
    }
  }

  cat("\nStep 5: Final results...\n")
  cat("Results dimensions:", nrow(results), "x", ncol(results), "\n")

  if (nrow(results) > 0) {
    if (args$debug) {
      cat("Results preview:\n")
      print(head(results[, 1:min(6, ncol(results))]))
      if (ncol(results) > 6) cat("... (", ncol(results)-6, "more columns)\n")
    }

    # Save results
    write.table(results, args$output, sep='\t', quote=FALSE, row.names=FALSE, na='NA')
    cat("\nConversion completed successfully!\n")
    cat("Output saved to:", args$output, "\n")
    cat("Number of genes processed:", nrow(results), "\n")
    cat("Number of groups:", length(cell_types), "\n")
    cat("Groups:", paste(cell_types, collapse=", "), "\n")

    # Show column names for verification
    cat("\nOutput columns:\n")
    for (i in 1:ncol(results)) {
      cat(i, ":", colnames(results)[i], "\n")
    }
  } else {
    cat("No results generated\n")
  }
} else {
  cat("No PAS columns parsed successfully\n")
}
