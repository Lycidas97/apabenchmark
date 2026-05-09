library(optparse)

# Define the options
option_list <- list(
    make_option(c("-b", "--bamfile"), type="character", default=NULL,
                help="Path to the BAM file", metavar="path"),
    make_option(c("-c", "--cellbarcode"), type="character", default=NULL,
                help="Path to the cellbarcode file", metavar="path"),
    make_option(c("-w", "--workdir"), type="character", default=NULL,
                help="Working directory path", metavar="path"),
    make_option(c("-r", "--reference"), type="character", default=NULL,
                help="Path to the reference gtf file", metavar="path"),
    make_option(c("-t", "--thread"), type="numeric", default=1,
                help="Number of threads [default %default]", metavar="number"),
    make_option(c("--cb_tag"), type="character", default="CB",
                help="Cell barcode tag in the BAM file [default %default]", metavar="tag"),
    make_option(c("--umi_tag"), type="character", default="UB",
                help="UMI tag in the BAM file [default %default]", metavar="tag")
)

# Parse the options
parser <- OptionParser(
  option_list=option_list,
  description = "This script runs the Sierra pipeline. Save PAS site file to workdir/pas.bed and PAS count file to workdir/pas_counts.tsv."
  )
args <- parse_args(parser)

# Assign the parsed options to variables
bamfile <- args$bamfile
cellbarcode <- args$cellbarcode
workdir <- args$workdir
reference <- args$reference
thread <- args$thread
cb_tag <- args$cb_tag
umi_tag <- args$umi_tag


# # test paramaters
# bamfile <- "/path/to/apabenchmark_final/simulation/bam/Chromium_mouse_olfactorybulb_GSM3449591_pas1_rep1_1.bam"
# workdir <- "/path/to/apabenchmark_final/benchmark/result/sierra/test"
# reference <- "/path/to/gencode.vM25.annotation.gtf"
# thread <- 8
# junctions.file <- "junctions.bed"
# cellbarcode <- "/path/to/apabenchmark_final/simulation/bam/Chromium_mouse_olfactorybulb_GSM3449591_pas1_rep1_1.barcode_list.txt"

# Optional: Check if any of the mandatory options are missing and handle accordingly
if (is.null(bamfile) || is.null(workdir) || is.null(reference) || is.null(cellbarcode)){
    stop("Mandatory options missing. Please check the input arguments.")
}

if (!file.exists(workdir)) {
  dir.create(workdir, recursive = TRUE)
}

setwd(workdir)

library(Rsamtools)
bam <- scanBamHeader(bamfile)
header <- bam[[1]]$header
chroms <- labels(bam[[1]]$targets)
chroms <- chroms[grepl("^chr", chroms)]

temp_files <- vector("list", length(chroms))

for (i in seq_along(chroms)) {
    chrom <- chroms[i]
    temp_file <- paste0("temp_", chrom, ".bed")
    temp_files[[i]] <- temp_file

    command <- glue::glue('regtools junctions extract -s RF {bamfile} -o {temp_file} -r {chrom}')
    system(command, wait = TRUE)
}
# wait for files written to disk
Sys.sleep(10)

junctions.file <- "junctions.bed"
command <- glue::glue("cat ", paste(temp_files, collapse = " "), " > ", junctions.file)
system(command, wait = TRUE)
sapply(temp_files, unlink)
Sys.sleep(5)
if (!file.exists(junctions.file)) {
  stop("junctions.bed is not created. Please check the input arguments.")
}

# white_list <- "whitelist.tsv"
# if (!file.exists(white_list)) {
#   command <- glue::glue("samtools view", bamfile, "| awk -F '\t' '{{for(i=12;i<=NF;i++) if($i ~ /^CB:Z:/) print substr($i,6)}}' | sort | uniq > ", white_list)
#   system(command, wait = TRUE)
# }

# generate cellbarcode if not provided
# if [ -z "$cellbarcode" ]; then
#     cellbarcode=$outputdir/cellbarcode.tsv
#     samtools view $bamfile | awk -F '\t' '{for(i=12;i<=NF;i++) if($i ~ /^CB:Z:/) print substr($i,6)}' | sort | uniq > $cellbarcode
# fi

# if [[ $cellbarcode != /* ]]; then
#     cellbarcode=$PWD/$cellbarcode
# fi



library(Sierra)

peakfile <- "tmp.peaks.tsv"
output <- "pas.bed"

FindPeaks(output.file = peakfile,   
          gtf.file = reference,           
          bamfile = bamfile,                
          junctions.file = junctions.file,     
          ncores = thread)

CountPeaks(peak.sites.file = peakfile, 
           gtf.file = reference,
           bamfile = bamfile, 
           whitelist.file = cellbarcode,
           CBtag = cb_tag,
           UMItag = umi_tag,
           output.dir = workdir, 
           ncores = thread,
           countUMI=FALSE)

mtx <- Matrix::readMM("matrix.mtx.gz")
barcodes <- read.delim("barcodes.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
sitenames <- read.delim("sitenames.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
rownames(mtx) <- sitenames[,1]
colnames(mtx) <- barcodes[,1]
mtx_transposed <- t(as.matrix(mtx))
mtx_pas_df <- as.data.frame(mtx_transposed)
write.table(mtx_pas_df, file = "pas_counts.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


tryCatch({
  pas_df <- read.table(peakfile,header = T)
  pas_df$strand<-ifelse(pas_df$Strand=="1","+","-")
  pas_df$start <- ifelse(pas_df$strand == "+", pas_df$Fit.end, pas_df$Fit.start - 1)
  pas_df$end <- ifelse(pas_df$strand == "+", pas_df$start + 1, pas_df$Fit.start)
  pas_df$score = ifelse(pas_df$strand == "+", pas_df$start, pas_df$end)
  pas_df$chr <- pas_df$Chr
  pas_bed <- pas_df[, c("chr", "start", "end", "polyA_ID", "score", "strand")]
  pas_bed <- pas_bed[order(pas_bed$chr, pas_bed$start), ]
},
  error = function(e) {
    print(e)
    assign('pas_bed', 
           data.frame(chr = character(),
                      start = integer(),
                      end = integer(),
                      PeakID = character(),
                      pas_position = integer(),
                      strand = character()),
           envir = .GlobalEnv)
  })

write.table(pas_bed, file = output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

