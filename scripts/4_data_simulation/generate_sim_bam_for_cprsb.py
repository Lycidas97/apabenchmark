import sys

sys.path.append("./apasim_dev/")
from apasim_dev import generate_sim_bam
import pandas as pd
from argparse import ArgumentParser
import os
import pickle

parser = ArgumentParser()
parser.add_argument(
    "-g", "--genome", dest="genome", help="path to genome fasta file", type=str
)
parser.add_argument(
    "-p", "--peak", dest="peak", help="path to peak pickle file", type=str
)
parser.add_argument(
    "-a",
    "--annotation",
    dest="annotation",
    help="path to pas annotation file",
    type=str,
)
parser.add_argument("-b", "--bam", dest="bam", help="path to output BAM file", type=str)
parser.add_argument("--barcode_num", dest="barcode_num", help="number of barcode", type=int)

if __name__ == "__main__":
    args = parser.parse_args()
    genome = args.genome
    peak = args.peak
    pas_annotation = args.annotation
    bam = args.bam
    barcode_num = args.barcode_num
    if not os.path.exists(peak):
        raise FileNotFoundError(f"Peak file {peak} not found")
    if not os.path.exists(genome):
        raise FileNotFoundError(f"Genome file {genome} not found")
    if not os.path.exists(pas_annotation):
        raise FileNotFoundError(f"Annotation file {pas_annotation} not found")
    # create bam dir if not exist
    if not os.path.exists(os.path.dirname(bam)):
        print(f"Creating directory {os.path.dirname(bam)}")
        os.makedirs(os.path.dirname(bam))

    with open(peak, "rb") as f:
        peaks = pickle.load(f)

    pas_df = pd.read_csv(pas_annotation, sep="\t")
    generate_sim_bam(
        genome,
        pas_df,
        peaks,
        bam,
        cell_num_per_type=int(barcode_num / 2),
    )
"""
python generate_sim_bam_for_cprsb.py -g /root/REFERENCE/GENOME/MOUSE/vm25/GRCm38.p6.genome.fa \
    -p /root/apabenchmark/simulation/sim_peak/simulated_peak_ps1000_rl40.pickle \
    -a /root/apabenchmark/simulation/pas_performance/mm10_sim_pas_gn2000_rep1.bed \
    -b /root/apabenchmark/simulation/mm10_sim_rep1.bam \
    --barcode_num 10000
"""
