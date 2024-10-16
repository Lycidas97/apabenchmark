import sys
import pysam
import os 
sys.path.append("../apasim/")
from find_cs import *
import pandas as pd
from argparse import ArgumentParser
import os
import pickle
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-b", "--bam", dest="bam", help="path to config file", type=str)
parser.add_argument("-g", "--gtf", dest="gtf", help="gtf", type=str)
parser.add_argument("-p", "--prefix", dest="prefix", help="prefix", type=str,)
parser.add_argument("-j", "--n_jobs", dest="n_jobs", help="number of jobs", type=int, default=8)


args = parser.parse_args()
bam_path = args.bam
gtf_path = args.gtf
prefix = args.prefix
n_jobs = args.n_jobs

max_gap = 10
extension_length = 1000
allowed_mismatches_at_start = 2
match_length = 8
min_required_matches = 7

bam = pysam.AlignmentFile(bam_path, 'rb')

print('Reading GTF file...')
utr_bed = BedTool("/root/nfsdata/apabenchmark/data/utr_extend/mm10.3utr.sorted.bed")
# genecoord_bed = get_gene_coordinates_bed(gtf_path, extension_length)
for chr_name in bam.references:
    if chr_name.isdigit():  # 检查染色体名字是否为数字
        chrom_list = [x for x in bam.references if x.isdigit()]
        chrom_list.append('X')
        chrom_list.append('Y')
        chrom_starts_with_chr = False
    elif chr_name.startswith('chr'):
        chrom_list = [x for x in bam.references if x.startswith('chr')]
        chrom_starts_with_chr = True
        break

cs_bed_df, onsite_bed_df = find_cs(
    bam_path, 
    chrom_list,
    chrom_starts_with_chr, 
    utr_bed, 
    max_gap, 
    n_jobs, 
    allowed_mismatches_at_start, 
    match_length, 
    min_required_matches
    )

cs_bed_df.to_csv(prefix + '_cs.bed', sep='\t', index=False, header=False)
onsite_bed_df.to_csv(prefix + '_onsite.bed', sep='\t', index=False, header=False)