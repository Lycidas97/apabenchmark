import sys

sys.path.append("../apasim/")
from simulation import generate_sim_bam
import pandas as pd
from argparse import ArgumentParser
import os
import pyarrow.feather as feather

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


if __name__ == "__main__":
    args = parser.parse_args()
    genome = args.genome
    peak = args.peak
    pas_annotation = args.annotation
    bam = args.bam

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

    # with open(peak, "rb") as f:
    #     peak_list = pickle.load(f)

    peaks_df = feather.read_feather(peak)
    peak_list = [group for _, group in peaks_df.groupby("pas")]
    pas_df = pd.read_csv(pas_annotation, sep="\t")

    # peak_list_filtered = [x for x in peak_list if (len(x) >= 20)]
    # peaks = []
    # drop_peaks = []
    # for x in peak_list_filtered:
    #     if x.strand.iloc[0] == "+":
    #         potential_tail = x.cigar_string.str.contains(r"[0-9]{2}S$")
    #     else:
    #         potential_tail = x.cigar_string.str.match(r"^[0-9]{2}S")
    #     potential_pas = x[potential_tail].groupby("distance_to_pas").size()
    #     potential_pas = potential_pas[potential_pas >= 3]
    #     if len(potential_pas) > 0:
    #         pas = potential_pas.idxmax()
    #         if abs(pas) <= 20:
    #             x["distance_to_pas"] = x["distance_to_pas"] - pas
    #             x = x[x["distance_to_pas"] <= 20]
    #             peaks.append(x)
    #         else:
    #             drop_peaks.append(x)
    #     else:
    #         if (x.distance_to_pas > 20).sum() <= 1:
    #             x = x[x["distance_to_pas"] <= 20]
    #             peaks.append(x)
    #         else:
    #             drop_peaks.append(x)

    peaks = [x[x["distance_to_pas"] <= 20] for x in peak_list if (len(x[x["distance_to_pas"] <= 20]) >= 50)]

    generate_sim_bam(
        genome,
        pas_df,
        peaks,
        bam,
    )
"""
python generate_sim_bam.py -g /root/REFERENCE/GENOME/MOUSE/vm25/GRCm38.p6.genome.fa \
    -p /root/apabenchmark/data/peaks/Visium_mouse_brain_V19L01-041-C1_peak.pickle \
    -a /root/apabenchmark/simulation/pas/mm10_sim_pas_rep1.bed \
    -b /root/apabenchmark/simulation/mm10_sim_rep1.bam
"""
