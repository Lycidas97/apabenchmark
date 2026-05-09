library(optparse)

# Define the options
option_list <- list(
    make_option(c("-b", "--bamfile"), type="character", default=NULL,
                help="Path to the BAM file", metavar="path"),
    make_option(c("-w", "--workdir"), type="character", default=NULL,
                help="Working directory path", metavar="path"),
    make_option(c("-r", "--readlength"), type="numeric", default=NULL,
                help="Length of the read", metavar="number"),
    make_option(c("-t", "--thread"), type="numeric", default=1,
                help="Number of threads [default %default]", metavar="number"),
    make_option(c("-c", "--celltag"), type="character", default=NULL,
                help="Cell tag in BAM file", metavar="character"),
    make_option(c("-u", "--umitag"), type="character", default=NULL,
                help="UMI tag in BAM file", metavar="character")
)

# Parse the options
parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

# Assign the parsed options to variables
bamfile <- args$bamfile
workdir <- args$workdir
readlength <- args$readlength
thread <- args$thread
celltag <- args$celltag
umitag <- args$umitag

# bamfile <- "/path/to/apabenchmark_final/simulation/bam/Chromium_mouse_olfactorybulb_GSM3449591_pas1_rep1_1.bam"
# cellbarcode <- "/path/to/apabenchmark_final/simulation/bam/Chromium_mouse_olfactorybulb_GSM3449591_pas1_rep1_1.barcode_list.txt"
# workdir <- "/path/to/apabenchmark_final/benchmark/result/scapatrap/test"
# readlength <- 120
# thread <- 8
# celltag <- "CB"
# umitag <- "UB"

# Optional: Check if any of the mandatory options are missing and handle accordingly
if (is.null(bamfile) || is.null(workdir) || is.null(readlength)) {
    stop("Mandatory options missing. Please check the input arguments.")
}


if (!file.exists(workdir)) {
  dir.create(workdir, recursive = TRUE)
}

setwd(workdir)

library(scAPAtrap)
library(Rsamtools)

countPeaks <- function(umitools.path, input, outputdir, celltag=NULL, umitag=NULL, ...) {

  fn='countPeaks'
  scAPAtrap:::.msg(sprintf("%s: start.", fn), preLine=FALSE, ...)

  output <- paste0(outputdir, "/counts.tsv.gz")
  if (is.null(celltag) | is.null(umitag)) {
    command <- paste0(umitools.path, " count --per-gene --gene-tag=XT --assigned-status-tag=XS --method=unique --per-cell -I ",
        input, " -S ", output)
  } else {
    command <- paste0(umitools.path, " count --per-gene --gene-tag=XT --assigned-status-tag=XS --method=unique --extract-umi-method=tag --umi-tag=", umitag, " --cell-tag=", celltag, " --per-cell -I ",
        input, " -S ", output)
  }
  scAPAtrap:::.runCommand(fn, command, ofile=output, ocheck=TRUE, ...)
  scAPAtrap:::.addTrapFiles(output, 'countPeaks: umitools counts')
  scAPAtrap:::.msg(sprintf("%s: finish.", fn), ...)
  return(output)
}


maxwidth <- 1000
cov.cutoff <- 10
min.cells <- 10

logf='process.log'

chrs=scAPAtrap:::.getBAMchrs(bamfile)
if (file.exists(logf)) {
  file.remove(logf)
}
unlink(logf)

samtools.path <- "samtools"
featureCounts.path <- 'featureCounts'
umi_tools.path <- 'umi_tools'
featureCounts.path <- 'featureCounts'

unlink('merged.bam')
unlink('merged.bam.bai')
system(command = glue::glue('ln -s {bamfile} merged.bam'))
system(command = glue::glue('ln -s {bamfile}.bai merged.bam.bai'))

input = 'merged.bam'

input.seps <- separateBamBystrand(samtools.path,
                input=input,
                thread=thread,
                logf=logf)

forwardPeaksFile=gsub(".bam", ".bed", input.seps[1])
forwardPeaksFile <-findPeaksByStrand(bamFile=input.seps[1],
    chrs=chrs,
    strand='+',
    L=readlength,
    maxwidth=maxwidth,
    cutoff=cov.cutoff,
    ofile=forwardPeaksFile,
    logf=logf)

reversePeaksFile=gsub(".bam", ".bed", input.seps[2])
reversePeaksFile <-findPeaksByStrand(bamFile=input.seps[2],
    chrs=chrs,
    strand='-',
    L=readlength,
    maxwidth=maxwidth,
    cutoff=cov.cutoff,
    ofile=reversePeaksFile,
    logf=logf)

peaksfile <- generateSAF(forwardPeaks=forwardPeaksFile,
      reversePeaks=reversePeaksFile,
      outputdir=workdir,
      logf=logf)



tryCatch({
    peaks <- read.table(peaksfile, sep="\t", header=FALSE)
    colnames(peaks) <- c("PeakID", "chr", "peak_start", "peak_end", "strand")

    peaks$start <- ifelse(peaks$strand == "+", peaks$peak_end, peaks$peak_start - 1)
    peaks$end <- ifelse(peaks$strand == "+", peaks$start + 1, peaks$peak_start)
    peaks$score = ifelse(peaks$strand == "+", peaks$start, peaks$end)

    pas_bed <- peaks[, c("chr", "start", "end", "PeakID", "score", "strand")]
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

pas_output <- "pas.bed"
write.table(pas_bed, file = pas_output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

final.bam <- generateFinalBam(featureCounts.path,
      samtools.path,
      input=input,
      peakfile=peaksfile,
      thread=thread,
      logf=logf)

scAPAtrap:::.checkAndPrintFiles(varnames='final.bam',
    filenames=final.bam, logf=logf)

countsfile <- countPeaks(umi_tools.path,
    input=final.bam,
    outputdir=workdir,
    celltag=celltag,
    umitag=umitag,
    logf=logf)

# reshape matrixs
library(data.table)
library(R.utils)
library(Matrix)
countslist <- fread(gunzip(countsfile, temporary = TRUE), sep="\t", header=TRUE, check.names=TRUE)

cells <- unique(countslist$cell)
genes <- unique(countslist$gene)

row_indices <- match(countslist$cell, cells)
col_indices <- match(countslist$gene, genes)
values <- countslist$count

pas_mtx <- sparseMatrix(i = row_indices, j = col_indices, x = values, 
                        dims = c(length(cells), length(genes)),
                        dimnames = list(cells, genes))



num_rows <- nrow(pas_mtx)

block_size <- min(1000, num_rows)

num_row_blocks <- ceiling(num_rows / block_size)


output_file <- "pas_counts.tsv"
cat(paste(colnames(pas_mtx), collapse = "\t"), file = output_file, append = FALSE)
cat("\n", file = output_file, append = TRUE)


for (i in 1:num_row_blocks) {
  row_start <- (i - 1) * block_size + 1
  row_end <- min(i * block_size, num_rows)
  
  block <- pas_mtx[row_start:row_end, ]
  dense_block <- as.matrix(block)
  dense_block <- cbind(rownames(block), dense_block)
  write.table(dense_block, file = output_file, sep = "\t", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)
}


# countslist <- read.table(countsfile, sep="\t", header=TRUE, check.names=TRUE)
# library(reshape2)
# pas_mtx <- acast(countslist, formula = cell ~ gene, value.var = "count", fill = 0)
# write.table(pas_mtx, file = "pas_counts.tsv", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

bam_files <- list.files(path = workdir, pattern = "\\.bam$", full.names = TRUE)
bai_files <- list.files(path = workdir, pattern = "\\.bai$", full.names = TRUE)
file_to_remove <- c(bam_files, bai_files, forwardPeaksFile, reversePeaksFile, peaksfile)
file.remove(file_to_remove)
# outputF <- paste0(gsub(".bam", "", input), ".forward.bam")
# outputR <- paste0(gsub(".bam", "", input), ".reverse.bam")

# samtools.Fcommand <- paste0(samtools.path, " view -@ ", thread,
#                             " -h -F 16 -bS ", input, " > ", outputF)
# system(command = samtools.Fcommand, wait = T)
# samtools.Rcommand <- paste0(samtools.path, " view -@ ", thread,
#                             " -h -f 16 -bS ", input, " > ", outputR)
# system(command = samtools.Rcommand, wait = T)

# command <- paste0(samtools.path, " index -@ ", thread, " ",
#                   outputF)
# system(command = command, wait = T)
# command <- paste0(samtools.path, " index -@ ", thread, " ",
#                   outputR)
# system(command = command, wait = T)

# nextinput <- c('merged.forward.bam', 'merged.reverse.bam')



# maxwidth <- 1000

# fullcovF <- loadBpCoverages(nextinput[1],chrs)
# fullcovR <- loadBpCoverages(nextinput[2],chrs)

# fullcovF <- fullcovF[!is.na(names(fullcovF))] # remove NA name columns in list
# fullcovR <- fullcovR[!is.na(names(fullcovR))]

# forwardPeaks <-findPeaks(fullcovF, '+', read_length, maxwidth)
# reversePeaks <-findPeaks(fullcovR, '-', read_length, maxwidth)

# forwardPeaks$start <- forwardPeaks$end
# forwardPeaks$end <- forwardPeaks$end + 1
# reversePeaks$end <- reversePeaks$start
# reversePeaks$start <- reversePeaks$start - 1
# forwardPeaks$pas_position <- forwardPeaks$start
# reversePeaks$pas_position <- reversePeaks$end

# peaks <- rbind(forwardPeaks, reversePeaks)

# tryCatch({
#   peaks$PeakID <- paste0("peak", "_", 1:length(peaks$chr))
#   peaks.bed <- peaks[, c("chr", "start", "end", "PeakID", "pas_position","strand")]
#   peaks.bed <- peaks.bed[order(peaks.bed$chr, peaks.bed$start), ]
# },
#   error = function(e) {
#     print(e)
#     assign('peaks.bed', 
#            data.frame(chr = character(),
#                       start = integer(),
#                       end = integer(),
#                       PeakID = character(),
#                       pas_position = integer(),
#                       strand = character()),
#            envir = .GlobalEnv)
#   })


# output <- "pas.bed"
# write.bed(peaks.bed, output)




