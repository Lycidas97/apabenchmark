library(optparse)
option_list <- list(
    make_option(c("-b", "--bamfile"), type="character", default=NULL,
                help="Path to the BAM file", metavar="path"),
    make_option(c("-c", "--cellbarcode"), type="character", default=NULL,
                help="Path to the cellbarcode file", metavar="path"),
    make_option(c("-w", "--workdir"), type="character", default=NULL,
                help="Working directory path", metavar="path"),
    make_option(c("-g", "--genome"), type="character", default=NULL,
                help="Genome version, only support mm10, mm9, hg38, hg19", metavar="character"),
    make_option(c("-t", "--thread"), type="numeric", default=4,
                help="Number of threads [default %default]", metavar="number"),
    make_option(c("-l", "--read_length"), type="numeric", default=150,
                help="Read length [default %default]", metavar="number"),
    make_option(c("--cb_tag"), type="character", default="CB",
                help="Cell barcode tag [default %default]", metavar="character")
)

options <- parse_args(OptionParser(option_list=option_list))

genome <- options$genome
bamfile <- options$bamfile
workdir <- options$workdir
thread <- options$thread
cellbarcode <- options$cellbarcode
read_length <- options$read_length
cb_tag <- options$cb_tag

require(package = "dplyr", warn.conflicts = F)
require(package = "tidyr", warn.conflicts = F)
require(package = "ggplot2", warn.conflicts = F)
require(package = "EnvStats", warn.conflicts = F)
require(package = "parallel", warn.conflicts = F)
require(package = "mclust", warn.conflicts = F)
require("Rsubread")
require("scAPA")
require("org.Hs.eg.db")

# check genome
if(genome %in% c("mm10", "mm9", "hg38", "hg19")) {
    print(paste0("Valid genome version ", genome))
} else {
    stop("Invalid genome version. Only support mm10, mm9, hg38, hg19")
}

# bamfile <- "/path/to/apabenchmark_final/simulation/bam/Chromium_mouse_olfactorybulb_GSM3449591_pas1_rep1_1.bam"
# workdir <- "/path/to/apabenchmark_final/benchmark/result/scapa/test1"
# threads <- 8
# genome <- "Mm10"
# cellbarcode <- "/path/to/apabenchmark_final/simulation/bam/Chromium_mouse_olfactorybulb_GSM3449591_pas1_rep1_1.barcode_list.txt"
# read_length <- 120
if (!file.exists(workdir)) {
    dir.create(workdir)
    print("Warning: working directory does not exist. Creating a new one.")
}
setwd(workdir)

# liftover_chain_path <- "/path/to/apabenchmark_final/benchmark/bin/liftover"
liftover_chain_path <- "/"

# generate genome utr annotation file
if (genome == "mm10") {
  data(mm10.utr)
  write.bed(.x = mm10.utr, f = paste0("./","threeprimeUTRs_new.bed"))
  data(mm10.introns)
  write.bed(.x = mm10.introns, f = paste0("./","introns_new.bed"))
  utrs <- mm10.utr
}
if (genome == "mm9") {
  data(mm10.utr)
  mm10.utr$V1 <- paste0("chr", mm10.utr$V1)
  write.bed(.x = mm10.utr, f = paste0("./","threeprimeUTRs.bed"))
  liftover.command <- paste0("liftOver ./","threeprimeUTRs.bed ",liftover_chain_path,"/mm10ToMm9.over.chain"," ./","threeprimeUTRs_new.bed ./","umap.bed")
  system(command = liftover.command, wait = T)
  library(rtracklayer)
  utrs <- import.bed(paste0("./","threeprimeUTRs_new.bed"))

  data(mm10.introns)
  mm10.introns$V1 <- paste0("chr", mm10.introns$V1)
  write.bed(.x = mm10.introns, f = paste0("./","introns.bed"))
  liftover.command.introns <- paste0("liftOver ./","introns.bed ",liftover_chain_path,"/mm10ToMm9.over.chain"," ./","introns_new.bed ./","umap_introns.bed")
  system(command = liftover.command.introns, wait = T)

  threeprimeUTRs_new_bed <- read.table("threeprimeUTRs_new.bed", header = F, stringsAsFactors = F)
  threeprimeUTRs_new_bed$V1 <- as.character(threeprimeUTRs_new_bed$V1)
  threeprimeUTRs_new_bed$V1 <- sub("chr", "", threeprimeUTRs_new_bed$V1)
  write.table(threeprimeUTRs_new_bed, file = "threeprimeUTRs_new.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
 
  introns_new_bed <- read.table("introns_new.bed", header = F, stringsAsFactors = F)
  introns_new_bed$V1 <- as.character(introns_new_bed$V1)
  introns_new_bed$V1 <- sub("chr", "", introns_new_bed$V1)
  write.table(introns_new_bed, file = "introns_new.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}
if (genome == "hg19") {
  data(hg19.utr)
  hg19.utr$V1 <- paste0("chr", hg19.utr$V1)
  write.table(x = hg19.utr, file = paste0("./","threeprimeUTRs_new.bed"), 
              sep = "\t", quote = F, row.names = F, col.names = F)
  utrs <- hg19.utr
  data(hg19.introns)
  write.bed(.x = hg19.introns, f = paste0("./","introns_new.bed"))
}
if (genome == "hg38") {
  data(hg19.utr)
  hg19.utr$V1 <- paste0("chr", hg19.utr$V1)
  write.table(x = hg19.utr, file = paste0("./","threeprimeUTRs.bed"), 
              sep = "\t", quote = F, row.names = F, col.names = F)
  liftover.command <- paste0("liftOver ./","threeprimeUTRs.bed ",liftover_chain_path,"/hg19ToHg38.over.chain"," ./","threeprimeUTRs_new.bed ./","umap.bed")
  system(command = liftover.command, wait = T)
  library(rtracklayer)
  utrs <- import.bed(paste0("./","threeprimeUTRs_new.bed"))
  data(hg19.introns)
  hg19.introns$V1 <- paste0("chr", hg19.introns$V1)
  write.bed(.x = hg19.introns, f = paste0("./","introns.bed"))
  liftover.command.introns <- paste0("liftOver ./","introns.bed ",liftover_chain_path,"/hg19ToHg38.over.chain"," ./","introns_new.bed ./","umap_introns.bed")
  system(command = liftover.command.introns, wait = T)
  
  threeprimeUTRs_new_bed <- read.table("threeprimeUTRs_new.bed", header = F, stringsAsFactors = F)
  threeprimeUTRs_new_bed$V1 <- as.character(threeprimeUTRs_new_bed$V1)
  threeprimeUTRs_new_bed$V1 <- sub("chr", "", threeprimeUTRs_new_bed$V1)
  write.table(threeprimeUTRs_new_bed, file = "threeprimeUTRs_new.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
 
  introns_new_bed <- read.table("introns_new.bed", header = F, stringsAsFactors = F)
  introns_new_bed$V1 <- as.character(introns_new_bed$V1)
  introns_new_bed$V1 <- sub("chr", "", introns_new_bed$V1)
  write.table(introns_new_bed, file = "introns_new.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}



unlink('merged.bam')
unlink('merged.bam.bai')
system(command = glue::glue('ln -s {bamfile} merged.bam'))
system(command = glue::glue('ln -s {bamfile}.bai merged.bam.bai'))

write_log_start(f = "./scAPA.script.log", "Splitting bams", command = NA)
if (!dir.exists("CellBams"))
    dir.create("CellBams")
list.cells <- read.delim(file = cellbarcode, header = F)
list.cells <- data.frame(cell = as.character(list.cells[,1]), bam = "merged.bam")
list.cells <- split(x = list.cells, f = list.cells$cell, drop = T)

split_bams <- function(x) {
    cell <- x["cell"]
    input_bam <- x["bam"]
    output_bam <- paste0(cell, ".bam")
    FilterBAMbyTag.command <- paste0("FilterBamByTag ",
                                        "TAG=",cb_tag," I=",
                                        input_bam, " O=",
                                        "./CellBams/", output_bam,
                                        " TAG_VALUE=", cell)
    # print(FilterBAMbyTag.command)
    system(command = FilterBAMbyTag.command, wait = T)
}
if(thread > 1){
parallel::mclapply(X = list.cells, FUN = split_bams, mc.cores = thread,
                    mc.preschedule = T)
}
if (thread ==1){
lapply(X = list.cells, FUN = split_bams)
}
write_log(f = "./scAPA.script.log", stage = "Splitting bams")

# find peaks
write_log_start(f = "./scAPA.script.log", "Find Peaks", command = NA)
system(command = paste0("echo \'makeTagDirectory ","Tagdirectory <(samtools view -h ", "merged.bam)\'|bash"),wait = T)
system(command = paste0("findPeaks ./","Tagdirectory -size 50 -fragLength {read_length} -minDist 1 -strand separate -o ","Peakfile"),wait = T)
write_log(f = "./scAPA.script.log", stage = "Find Peaks")

write_log_start(f = "./scAPA.script.log", "Merge Peaks", command = NA)

scAPA::merge_peaks(bedtools.path="", peaks.file="Peakfile", path= "./")
peak_bed <- read.table("merge.peakfile.bed", header = F, stringsAsFactors = F)
peak_bed$V1 <- as.character(peak_bed$V1)
data_starts_with_chr <- any(startsWith(peak_bed$V1, "chr"))
if (data_starts_with_chr){
    peak_bed$V1 <- sub("chr", "", peak_bed$V1)
}
write.table(peak_bed, file = "merge.peakfile.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write_log(f = "./scAPA.script.log", stage = "Merge Peaks")


write_log_start(f = "./scAPA.script.log", "Intersect Peaks", command = NA)
intersect_peaks <- function (path, bed.name, bedtools.path, introns = F) 
{
    if (!introns) {
        intersect.command <- paste0(bedtools.path, "bedtools intersect -wa ", 
            "-wb -s -a ./", bed.name, " -b ./","threeprimeUTRs_new.bed >", 
            " ./","annotatedpeaks")
        system(command = intersect.command, wait = T)
        peakstosort <- read.delim(file = paste0("./","annotatedpeaks"), 
            header = F)
        peakstosort <- peakstosort[, c(1, 2, 3, 10, 4, 6)]
        colnames(peakstosort) <- paste0("V", 1:6)
        peakstosort <- dplyr::group_by(peakstosort, V4)
        peakstosort <- dplyr::group_by(peakstosort, V4)
        peakstosort <- dplyr::mutate(peakstosort, desire = ifelse(V6 == 
            "+", dplyr::dense_rank(x = V3), dplyr::dense_rank(dplyr::desc(x = V3))))
        peakstosort$V4 <- paste0(as.character(peakstosort$V4), 
            "_", peakstosort$desire)
        peakstosort <- peakstosort[, -7]
        system(command = paste0("rm ./","annotatedpeaks"), wait = T)
    }
    if (introns) {
        intersect.command <- paste0(bedtools.path, "bedtools intersect -wa ", 
            "-wb -s -a ./", bed.name, " -b ./","introns_new.bed >", 
            " ./","annotatedpeaks")
        system(command = intersect.command, wait = T)
        peakstosort <- read.delim(file = paste0("./","annotatedpeaks"), 
            header = F)
        peakstosort <- peakstosort[, c(1, 2, 3, 10, 4, 6)]
        colnames(peakstosort) <- paste0("V", 1:6)
        peakstosort <- dplyr::group_by(peakstosort, V4)
        peakstosort <- dplyr::group_by(peakstosort, V4)
        peakstosort <- dplyr::mutate(peakstosort, desire = ifelse(V6 == 
            "+", dplyr::dense_rank(x = V3), dplyr::dense_rank(dplyr::desc(x = V3))))
        peakstosort$V4 <- paste0(as.character(peakstosort$V4), 
            "_", peakstosort$desire)
        peakstosort <- peakstosort[, -7]
        peakstosort <- rbind.data.frame(peakstosort, utrs)
        system(command = paste0("rm ./","annotatedpeaks"), wait = T)
    }
    peakstosort
}

peaks.bed <- intersect_peaks(bed.name = "merge.peakfile.bed",
                             path = "", bedtools.path = '',
                             introns = F)

if (data_starts_with_chr){
    peaks.bed$V1 <- paste0("chr", peaks.bed$V1)
}
write.bed(peaks.bed, "peaks.bed")

system(command = paste0("bedtools genomecov -ibam ", bamfile,  " -bg -strand + | awk 'BEGIN {OFS = \"\t\"}{print $1, $2, $3, $4, \".\", \"+\"}' > merged.wig"), wait=T)
system(command = paste0("bedtools genomecov -ibam ", bamfile,  " -bg -strand - | awk 'BEGIN {OFS = \"\t\"}{print $1, $2, $3, $4, \".\", \"-\"}' > merged.wig"), wait=T)
system(command = paste0("bedtools intersect -s -wb -b peaks.bed -a merged.wig > intersected.wig"), wait=T)

peaks.wig <- read.delim(file = "intersected.wig", header = F)
peaks.wig <- split(x = peaks.wig, f = peaks.wig$V10, drop = T)
bed <- plyr::rbind.fill(lapply(1:length(peaks.wig),
                                           FUN = scAPA::creat_mclus))                        
write.bed(bed, "intersected.bed")

tryCatch({
    df <- read.table(paste0("./","intersected.bed"),
        header = F,
        col.names = c("chr", "S", "E", "polyA_ID", "score", "strand"))
    df$start <- ifelse(df$strand == "+", df$E, df$S - 1)
    df$end <- ifelse(df$strand == "+", df$start + 1, df$S)
    df$score = ifelse(df$strand == "+", df$start, df$end)
    pas_bed <- df[, c("chr", "start", "end", "polyA_ID", "score", "strand")]
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
output <- "pas.bed"
write.table(pas_bed, file = output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write_log(f = "./scAPA.script.log", stage = "Intersect Peaks")







# Count reads --------------------------------------------------------------
utr.saf <- bed[, c(4, 1, 2, 3, 6)]
colnames(utr.saf) <- c("GeneID", "Chr", "Start", "End", "Strand")
write_log_start("Stage 2: Quantifying peak usage\n\n", command = NA,
                f = "./scAPA.script.log")
write_log_start(f = "./scAPA.script.log", "featureCounts for 3UTRs",
                command = NA)
bam.cluster.files <- list.files(path = "./CellBams", pattern = ".bam$",
                                full.names = F)
cellnames <- gsub(x = bam.cluster.files, pattern = ".bam", replacement = "")
bam.cluster.files <- list.files(path = "./CellBams", pattern = ".bam$",
                                full.names = T)
counts <- Rsubread::featureCounts(files = bam.cluster.files, isGTFAnnotationFile = F,
                        strandSpecific = 1, annot.ext = utr.saf,
                        largestOverlap = T, nthreads = thread)

write_log(f = "./scAPA.script.log", stage = "featureCounts for 3UTRs")


counts_transposed <- t(as.matrix(counts$counts))
rownames(counts_transposed) <- gsub(".bam", "", rownames(counts_transposed))
counts_transposed_df <- as.data.frame(counts_transposed)
write.table(counts_transposed_df, file = paste0(workdir, "/pas_counts.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

system(command = "rm -r ./CellBams")

write_log(f = "./scAPA.script.log", stage = "Write pas_counts.tsv")

