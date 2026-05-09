#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(argparse)
})

# Create command line arguments
parser <- ArgumentParser(description='Parameterized scMAPA analysis for differential APA detection')

# Input/Output parameters
parser$add_argument('-i', '--input_file', type='character', required=TRUE,
                    help='Input file (DaPars2 format or PAS converter output)')
parser$add_argument('-o', '--output_dir', type='character', default='scmapa_results',
                    help='Output directory for results')


# Group comparison parameters
parser$add_argument('--comparison_mode', type='character', default='single',
                    choices=c('single', 'batch'),
                    help='Comparison mode: single (one pair) or batch (multiple pairs)')
parser$add_argument('--group1', type='character', default=NULL,
                    help='First group name for single comparison')
parser$add_argument('--group2', type='character', default=NULL,
                    help='Second group name for single comparison')
parser$add_argument('--comparison_file', type='character', default=NULL,
                    help='Tab-separated file with group pairs for batch comparison (group1\\tgroup2)')

# Data filtering parameters
parser$add_argument('--NAcutoff', type='integer', default=2,
                    help='Minimum number of groups with non-NA expression (default: 2)')
parser$add_argument('--CPMcutoff_L', type='double', default=5,
                    help='CPM cutoff for long isoforms (default: 5)')
parser$add_argument('--CPMcutoff_S', type='double', default=5,
                    help='CPM cutoff for short isoforms (default: 5)')

# APA analysis parameters
parser$add_argument('--coverageCutoff', type='double', default=10,
                    help='Coverage cutoff for APA analysis (default: 10)')
parser$add_argument('--ORcutoff', type='double', default=0.3,
                    help='Odds ratio cutoff for APA analysis (default: 0.3)')
parser$add_argument('--adPval', type='double', default=0.1,
                    help='FDR cutoff for APA significance (default: 0.1)')

args <- parser$parse_args()

cat("=== Simplified scMAPA Analysis ===\n")
cat("Input file:", args$input_file, "\n")
cat("Output directory:", args$output_dir, "\n")

# Validate input parameters
if (args$comparison_mode == 'single') {
    if (is.null(args$group1) || is.null(args$group2)) {
        cat("Error: For single comparison mode, both --group1 and --group2 must be specified\n")
        quit(status=1)
    }
    if (args$group1 == args$group2) {
        cat("Error: group1 and group2 must be different\n")
        quit(status=1)
    }
} else if (args$comparison_mode == 'batch') {
    if (is.null(args$comparison_file)) {
        cat("Error: For batch comparison mode, --comparison_file must be specified\n")
        quit(status=1)
    }
    if (!file.exists(args$comparison_file)) {
        cat("Error: Comparison file does not exist:", args$comparison_file, "\n")
        quit(status=1)
    }
}

# Check input file
if (!file.exists(args$input_file)) {
    cat("Error: Input file does not exist:", args$input_file, "\n")
    quit(status=1)
}

# Create output directory
if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive=TRUE)
    cat("Created output directory:", args$output_dir, "\n")
}

# Load required packages (fail fast instead of runtime installation)
cat("\nChecking required packages...\n")
required_pkgs <- c("multcomp", "lmtest", "wec", "Matrix.utils", "nnet", "stringr", "scMAPA")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
    cat("Error: Missing required R packages:", paste(missing_pkgs, collapse = ", "), "\n")
    cat("Please pre-install required packages in the execution environment before running this script.\n")
    quit(status = 1)
}

suppressPackageStartupMessages({
    library(stringr)
    library(multcomp)
    library(lmtest)
    library(wec)
    library(Matrix.utils)
    library(scMAPA)
})

cat("Packages loaded successfully\n")

# Function to read and validate input file
create_empty_input_result <- function() {
    groups <- unique(c(args$group1, args$group2))
    groups <- groups[!is.na(groups) & nzchar(groups)]

    empty_data <- data.frame(
        Gene = character(),
        fit_value = numeric(),
        Predicted_Proximal_APA = numeric(),
        Loci = character(),
        stringsAsFactors = FALSE
    )

    for (group in groups) {
        empty_data[[paste0(group, "_long_exp")]] <- numeric()
        empty_data[[paste0(group, "_short_exp")]] <- numeric()
        empty_data[[paste0(group, "_PDUI")]] <- numeric()
    }

    list(data = empty_data, format = "dapars2")
}

read_input_data <- function(file_path) {
    cat("\nReading input data...\n")
    file_size <- file.info(file_path)$size
    if (is.na(file_size) || file_size == 0) {
        cat("Input file is empty; using an empty template for downstream zero-output handling.\n")
        return(create_empty_input_result())
    }

    data <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    cat("Input data dimensions:", nrow(data), "x", ncol(data), "\n")

    # Check if this is a PAS converter output (has cell_type column) or DaPars2 output
    if ("cell_type" %in% colnames(data)) {
        cat("Detected PAS converter output format\n")
        return(list(data=data, format="pas_converter"))
    } else {
        cat("Detected DaPars2 standard output format\n")
        return(list(data=data, format="dapars2"))
    }
}

# Function to extract group names from column names
extract_groups_from_columns <- function(col_names) {
    # Extract group names from columns like {group}_long_exp, {group}_short_exp
    long_cols <- col_names[str_detect(col_names, "_long_exp$")]
    groups <- str_extract(long_cols, "^(.+)_long_exp$")
    groups <- str_remove(groups, "_long_exp$")
    return(groups)
}

# Function to prepare comparison pairs
prepare_comparison_pairs <- function(groups) {
    if (args$comparison_mode == 'single') {
        # Single comparison
        if (!args$group1 %in% groups || !args$group2 %in% groups) {
            cat("Error: Specified groups not found in data\n")
            cat("Available groups:", paste(groups, collapse=", "), "\n")
            quit(status=1)
        }
        return(data.frame(group1=args$group1, group2=args$group2, stringsAsFactors=FALSE))
    } else {
        # Batch comparison from file
        comp_data <- read.delim(args$comparison_file, header=FALSE, stringsAsFactors=FALSE)
        colnames(comp_data) <- c("group1", "group2")

        # Validate groups
        invalid_groups <- c(comp_data$group1, comp_data$group2)[!c(comp_data$group1, comp_data$group2) %in% groups]
        if (length(invalid_groups) > 0) {
            cat("Error: Invalid groups in comparison file:", paste(unique(invalid_groups), collapse=", "), "\n")
            cat("Available groups:", paste(groups, collapse=", "), "\n")
            quit(status=1)
        }

        return(comp_data)
    }
}

# Function to format data for scMAPA
format_data_for_scmapa <- function(data, format_type, groups) {
    cat("Formatting data for scMAPA...\n")

    if (format_type == "pas_converter") {
        # Convert PAS converter output to DaPars2 format
        # PAS converter already has the right format: {group}_long_exp, {group}_short_exp
        formatted_data <- data
    } else {
        # Standard DaPars2 format - need to rename columns
        formatted_data <- data

        # Check if columns need renaming (common DaPars2 formats)
        for (group in groups) {
            # Handle various possible DaPars2 column naming patterns
            patterns <- list(
                paste0("^", group, "_long_exp$"),
                paste0("^Sample_", group, "_long_exp$"),
                paste0("^cluster", group, "_long_exp$")
            )

            for (pattern in patterns) {
                matching_cols <- grep(pattern, colnames(formatted_data), value=TRUE)
                if (length(matching_cols) > 0) {
                    # Rename to standard format
                    new_name <- paste0(group, "_long_exp")
                    if (matching_cols[1] != new_name) {
                        colnames(formatted_data)[colnames(formatted_data) == matching_cols[1]] <- new_name
                    }
                    break
                }
            }

            # Same for short_exp
            patterns_short <- list(
                paste0("^", group, "_short_exp$"),
                paste0("^Sample_", group, "_short_exp$"),
                paste0("^cluster", group, "_short_exp$")
            )

            for (pattern in patterns_short) {
                matching_cols <- grep(pattern, colnames(formatted_data), value=TRUE)
                if (length(matching_cols) > 0) {
                    new_name <- paste0(group, "_short_exp")
                    if (matching_cols[1] != new_name) {
                        colnames(formatted_data)[colnames(formatted_data) == matching_cols[1]] <- new_name
                    }
                    break
                }
            }
        }
    }

    # Verify all required columns exist
    required_cols <- c()
    for (group in groups) {
        required_cols <- c(required_cols, paste0(group, "_long_exp"), paste0(group, "_short_exp"))
    }

    missing_cols <- setdiff(required_cols, colnames(formatted_data))
    if (length(missing_cols) > 0) {
        cat("Error: Missing required columns:", paste(missing_cols, collapse=", "), "\n")
        quit(status=1)
    }

    return(formatted_data)
}

# Function to run scMAPA analysis for a single comparison
compute_pdui <- function(long_exp, short_exp, fallback = 0.5) {
    denominator <- long_exp + short_exp
    pdui <- long_exp / denominator
    pdui[is.na(pdui) | is.infinite(pdui)] <- fallback
    return(pdui)
}

run_scmapa_comparison <- function(data, group1, group2, comparison_name, args_obj) {
    cat("\nRunning scMAPA analysis for", comparison_name, "(", group1, "vs", group2, ")...\n")

    # Select only the two groups for comparison
    comparison_cols <- c("Gene", paste0(group1, "_long_exp"), paste0(group1, "_short_exp"),
                        paste0(group2, "_long_exp"), paste0(group2, "_short_exp"))

    # Check if all required columns exist
    missing_cols <- setdiff(comparison_cols, colnames(data))
    if (length(missing_cols) > 0) {
        stop(sprintf(
            "Missing columns for %s: %s",
            comparison_name,
            paste(missing_cols, collapse = ", ")
        ))
    }

    comparison_data <- data[, comparison_cols, drop=FALSE]
    if (nrow(comparison_data) == 0) {
        cat("No candidate genes available for", comparison_name, "- writing empty APA result.\n")
        return(list(
            APA_results = data.frame(Gene = character(), stringsAsFactors = FALSE),
            comparison_name = comparison_name
        ))
    }

    # Create temporary directory for this comparison
    temp_dir <- file.path(args_obj$output_dir, "temp", comparison_name)
    if (!dir.exists(temp_dir)) {
        success <- dir.create(temp_dir, recursive=TRUE, showWarnings=FALSE)
        if (!success) {
            stop("Failed to create directory: ", temp_dir)
        }
    }

    # Prepare data in the format expected by scMAPA
    scmapa_format_data <- comparison_data
    scmapa_format_data$fit_value <- 1.0  # Placeholder
    scmapa_format_data$Predicted_Proximal_APA <- 1.0  # Placeholder
    scmapa_format_data$Loci <- "chr1:1-1000"  # Placeholder

    # Calculate PDUI (long / (long + short))
    scmapa_format_data[[paste0(group1, "_PDUI")]] <- compute_pdui(
        scmapa_format_data[[paste0(group1, "_long_exp")]],
        scmapa_format_data[[paste0(group1, "_short_exp")]]
    )
    scmapa_format_data[[paste0(group2, "_PDUI")]] <- compute_pdui(
        scmapa_format_data[[paste0(group2, "_long_exp")]],
        scmapa_format_data[[paste0(group2, "_short_exp")]]
    )

    # Reorder columns to match scMAPA expected format
    scmapa_cols <- c("Gene", "fit_value", "Predicted_Proximal_APA", "Loci",
                    paste0(group1, "_long_exp"), paste0(group1, "_short_exp"), paste0(group1, "_PDUI"),
                    paste0(group2, "_long_exp"), paste0(group2, "_short_exp"), paste0(group2, "_PDUI"))

    scmapa_format_data <- scmapa_format_data[, scmapa_cols]

    # Make sure gene names are unique
    if (any(duplicated(scmapa_format_data$Gene))) {
        scmapa_format_data$Gene <- paste0(scmapa_format_data$Gene, "_", 1:nrow(scmapa_format_data))
    }

    # Save comparison data
    temp_file <- file.path(temp_dir, "PDUI_Comparison_1_All_Prediction_Results.txt")
    write.table(scmapa_format_data, temp_file, sep="\t", quote=FALSE, row.names=FALSE)

    # Run scMAPA analysis
    original_wd_local <- getwd()
    setwd(temp_dir)
    on.exit(setwd(original_wd_local), add=TRUE)

    tryCatch({
        ISOMatrix <- readinPAsites(
            path = ".",
            NAcutoff = args_obj$NAcutoff,
            CPMcutoff_L = args_obj$CPMcutoff_L,
            CPMcutoff_S = args_obj$CPMcutoff_S,
            clusterOfInterests = "all"
        )

        cat("Filtered genes:", nrow(ISOMatrix), "\n")

        if (nrow(ISOMatrix) >= 1) {
            # APA significance test
            APA_results <- APAtest(
                countMatrix = ISOMatrix,
                coverageCutoff = args_obj$coverageCutoff,
                ORcutoff = args_obj$ORcutoff,
                adPval = args_obj$adPval
            )

            cat("Significant APA genes:", nrow(APA_results), "\n")

            if (nrow(APA_results) > 0) {
                # Return only APA_results
                return(list(
                    APA_results = APA_results,
                    comparison_name = comparison_name
                ))
            } else {
                cat("No significant APA genes found\n")
                return(list(
                    APA_results = data.frame(Gene = character(), stringsAsFactors = FALSE),
                    comparison_name = comparison_name
                ))
            }
        } else {
            cat("Insufficient genes for analysis\n")
            return(list(
                APA_results = data.frame(Gene = character(), stringsAsFactors = FALSE),
                comparison_name = comparison_name
            ))
        }

    }, error = function(e) {
        stop("Analysis failed for ", comparison_name, ": ", e$message)
    })

}

# Function to save results
save_results <- function(results, comparison_name, args_obj) {
    cat("Saving APA results for", comparison_name, "...\n")

    result_dir <- file.path(args_obj$output_dir, comparison_name)
    if (!dir.exists(result_dir)) {
        dir.create(result_dir, recursive=TRUE, showWarnings=FALSE)
    }

    apa_results <- NULL
    if (!is.null(results) && !is.null(results$APA_results) && is.data.frame(results$APA_results)) {
        apa_results <- results$APA_results
    }
    if (is.null(apa_results) || ncol(apa_results) == 0) {
        apa_results <- data.frame(Gene = character(), stringsAsFactors = FALSE)
    }

    write.table(
        apa_results,
        file.path(result_dir, "APA_results.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    cat("Saved APA results for", comparison_name, "(rows:", nrow(apa_results), ")\n")

    cat("Results saved to:", result_dir, "\n")
}

# Main analysis pipeline
main_analysis <- function() {
    # Read input data
    input_result <- read_input_data(args$input_file)
    data <- input_result$data
    format_type <- input_result$format

    # Extract group names from columns
    groups <- extract_groups_from_columns(colnames(data))
    cat("Detected groups:", paste(groups, collapse=", "), "\n")

    # Prepare comparison pairs
    comparison_pairs <- prepare_comparison_pairs(groups)
    cat("Number of comparisons:", nrow(comparison_pairs), "\n")

    # Format data for scMAPA
    formatted_data <- format_data_for_scmapa(data, format_type, groups)

    # Run analyses
    for (i in 1:nrow(comparison_pairs)) {
        group1 <- comparison_pairs$group1[i]
        group2 <- comparison_pairs$group2[i]
        comparison_name <- paste0(group1, "_vs_", group2)

        results <- run_scmapa_comparison(formatted_data, group1, group2, comparison_name, args)

        save_results(results, comparison_name, args)
    }

    cat("\n=== Analysis completed ===\n")
    cat("APA results saved in:", args$output_dir, "\n")
}

# Run the analysis
main_analysis()
