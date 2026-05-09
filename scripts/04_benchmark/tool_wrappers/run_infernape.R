library(optparse)

# Define the options
option_list <- list(
    make_option(c("-b", "--bamfile"), type="character", default=NULL,
                help="Path to the BAM file", metavar="path"),
    make_option(c("-c", "--cellbarcode"), type="character", default=NULL,
                help="Path to the cellbarcode file", metavar="path"),
    make_option(c("-w", "--workdir"), type="character", default=NULL,
                help="Working directory path", metavar="path"),
    make_option(c("-g", "--genome"), type="character", default=NULL,
                help="Genome version, only support mm10, or hg38", metavar="character"),
    make_option(c("-t", "--thread"), type="numeric", default=4,
                help="Number of threads [default %default]", metavar="number")
)

# Usage

# Parse options
options <- parse_args(OptionParser(option_list=option_list))

genome <- options$genome
bamfile <- options$bamfile
workdir <- options$workdir
thread <- options$thread
cellbarcode <- options$cellbarcode


if (genome == "mm10") {
    genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    # genome.ref <- "/path/to/apabenchmark_final/benchmark/env/infernape/gencode.vM25.pcg.Infernape.csv"
    # pas_annotation <- "/path/to/apabenchmark_final/benchmark/env/infernape/mm10.PAS.csv"
    genome.ref <- "/gencode.vM25.pcg.Infernape.csv"
    pas_annotation <- "/mm10.PAS.csv"
} else if (genome == "hg38") {
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    # genome.ref <- "/path/to/apabenchmark_final/benchmark/env/infernape/gencode.v40.pcg.Infernape.csv"
    # pas_annotation <- "/path/to/apabenchmark_final/benchmark/env/infernape/hg38.PAS.csv"
    genome.ref <- "/gencode.v40.pcg.Infernape.csv"
    pas_annotation <- "/hg38.PAS.csv"
} else {
    stop("Invalid genome version. Only support mm9, mm10, hg19, hg38, or others.")
}

library(Infernape)
# genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
# genome.ref <- "/path/to/apabenchmark_final/benchmark/env/infernape/gencode.vM25.pcg.Infernape.csv"
# bamfile <- "/path/to/apabenchmark_final/simulation/bam/Chromium_mouse_olfactorybulb_GSM3449591_pas1_rep1_1.bam"
# pas <- "/path/to/apabenchmark_final/benchmark/env/infernape/mm10.PAS.csv"
# workdir <- "/path/to/apabenchmark_final/benchmark/result/infernape/test"
# threads <- 8
# cellbarcode <- "/path/to/apabenchmark_final/simulation/bam/Chromium_mouse_olfactorybulb_GSM3449591_pas1_rep1_1.barcode_list.txt"

if (!file.exists(workdir)) {
    dir.create(workdir)
    print("Warning: working directory does not exist. Creating a new one.")
}
setwd(workdir)
white_list <- read.table(cellbarcode, header = FALSE, stringsAsFactors = FALSE)
white_list$x <- white_list$V1
write.table(white_list, file = "white_list.txt", quote = FALSE, sep=",", row.names=FALSE)



ref.df <- utils::read.csv(genome.ref)
feature = chromosome = start = end = strand = gene_symbol = NULL
gene.ref.df <- ref.df %>% dplyr::filter(feature == 'gene') %>% dplyr::select(chromosome, start, end, strand, gene_symbol)
n.genes <- nrow(gene.ref.df)

Infernape_cnt(genome.ref = genome.ref,
              bam = bamfile,
              batch.start = 1,
              batch.end = n.genes,
              ncores = thread,
              d = 31,
              h = 5,
              d.cut = 50,
              hr = 160,
              min.mode.prop = 0.05,
              min.mode.cutoff = 5,
              output.path = './',
              pas.reference.file = pas_annotation,
              genome = genome,
              pas.search.cut.1 = 0,
              pas.search.cut.2 = 300,
              polystretch_length = 13,
              max_mismatch = 1,
              motif.search.cut = 300,
              invert_strand = FALSE,
              q = c(110, 200),
              whitelist.file = "white_list.txt",
              start.cid = NULL,
              end.cid = NULL
)

mtx <- Matrix::readMM(paste0(workdir, "/cnt_mat/matrix.mtx"))
barcodes <- read.delim(paste0(workdir,"/cnt_mat/barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE)
sitenames <- read.delim(paste0(workdir,"/cnt_mat/sitenames.tsv"), header = FALSE, stringsAsFactors = FALSE)
rownames(mtx) <- sitenames[,1]
colnames(mtx) <- barcodes[,1]
mtx_transposed <- t(as.matrix(mtx))
mtx_pas_df <- as.data.frame(mtx_transposed)
write.table(mtx_pas_df, file = paste0(workdir, "/pas_counts.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

tryCatch({
pas.df <- read.csv(paste0(workdir, "/anno_filtered.csv"))
pas.df$start <- ifelse(pas.df$strand == "+", pas.df$to, pas.df$from - 1)
pas.df$end <- ifelse(pas.df$strand == "+", pas.df$start + 1, pas.df$from)
pas.df$score = ifelse(pas.df$strand == "+", pas.df$start, pas.df$end)
pas.df$chr <- pas.df$seq
# pas.df$position <- as.numeric(pas.df$mode.pos)
# pas.df$start <- ifelse(pas.df$strand == "+", pas.df$position, pas.df$position - 1)
# pas.df$end <- ifelse(pas.df$strand == "+", pas.df$position + 1, pas.df$position)
pas.bed <- pas.df[, c("chr", "start", "end", "polyA_ID","score", "strand")]
pas.bed.output <- paste0(workdir, "/pas.bed")
pas.bed <- pas.bed[order(pas.bed$chr, pas.bed$start), ]
},
  error = function(e) {
    print(e)
    assign('pas.bed', 
           data.frame(chr = character(),
                      start = integer(),
                      end = integer(),
                      PeakID = character(),
                      score = integer(),
                      strand = character()),
           envir = .GlobalEnv)
  })

write.table(pas.bed, file = pas.bed.output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)