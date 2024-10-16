import pyarrow.feather as feather
import pandas as pd 
import numpy as np 
import sys
from pybedtools import BedTool
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest="input", help="path to input pickle file", type=str)
parser.add_argument("-o", "--output_dir", dest="output", help="directory to output", type=str)
parser.add_argument("-p", "--prefix", dest="prefix", help="prefix", type=str)
parser.add_argument("-a", "--annotation", dest="annotation", help="annotation", type=str)


utr_bed = BedTool("../../data/int_data/annotations/mm10.3utr.sorted.bed")


MIN_READS_PER_PEAK = 20
WINDOW_SIZE=20

args = parser.parse_args()
input_path = args.input
output_dir = args.output
prefix = args.prefix
annotation = args.annotation


peaks_df = feather.read_feather(input_path)
# with open(input_path, "rb") as f:
#     peaks_df = pickle.load(f)

annotation_bed = BedTool(annotation)

# peaks_df = pd.concat(peaks)
cs_counts = peaks_df[(peaks_df["distance_to_pas"] <= 20)].groupby("pas")["cigar_string"].count()
cs_counts = cs_counts[cs_counts >= MIN_READS_PER_PEAK]
filtered_cs = cs_counts.index
filtered_peaks_df = peaks_df[peaks_df["pas"].isin(filtered_cs)]
filtered_peaks_df["forward_unmatch"] = filtered_peaks_df["cigar_string"].str.endswith("S")
filtered_peaks_df["reverse_unmatch"] = filtered_peaks_df["cigar_string"].str.contains("^(\d+)S")
filtered_peaks_df["end_unmatch"] = (filtered_peaks_df["forward_unmatch"] & (filtered_peaks_df["strand"] == "+")) | (filtered_peaks_df["reverse_unmatch"] & (filtered_peaks_df["strand"] == "-"))
end_unmatch_counts = filtered_peaks_df[filtered_peaks_df["end_unmatch"] & (filtered_peaks_df["distance_to_pas"] >= -20) & (filtered_peaks_df["distance_to_pas"] <= 20)].groupby("pas")["distance_to_pas"].count()
valid_cs = end_unmatch_counts[end_unmatch_counts >= 2].index
cs_counts = cs_counts[cs_counts.index.isin(valid_cs)]
filtered_cs_df = pd.DataFrame(cs_counts).reset_index()
filtered_cs_df["chr"] = filtered_cs_df["pas"].str.split(":", expand=True)[0]
filtered_cs_df["start"] = filtered_cs_df["pas"].str.split(":", expand=True)[1].astype(int)
filtered_cs_df["strand"] = filtered_cs_df["pas"].str.split(":", expand=True)[3]
filtered_cs_df["end"] = filtered_cs_df["start"]+1
filtered_cs_df["name"] = filtered_cs_df["pas"]
filtered_cs_df["score"] = filtered_cs_df["cigar_string"]
filtered_cs_df["strand"] = filtered_cs_df["pas"].str.split(":", expand=True)[3]

if filtered_cs_df["chr"].str.startswith("chr").all():
    pass
else:
    filtered_cs_df["chr"] = "chr"+filtered_cs_df["chr"]
    
filtered_cs_df = filtered_cs_df[["chr", "start", "end", "name", "score", "strand"]].sort_values(["chr", "start"])
filtered_cs_bed = BedTool.from_dataframe(filtered_cs_df)

match_bed = filtered_cs_bed.window(annotation_bed, w=WINDOW_SIZE, sw=True, u=True)
unmatch_bed = filtered_cs_bed.window(annotation_bed, w=WINDOW_SIZE, sw=True, v=True)

match_utr_bed = match_bed.intersect(utr_bed, u=True,s=True)
match_nonutr_bed = match_bed.intersect(utr_bed, v=True,s=True)

unmatch_utr_bed = unmatch_bed.intersect(utr_bed, u=True,s=True)
unmatch_nonutr_bed = unmatch_bed.intersect(utr_bed, v=True,s=True)

filtered_cs_bed.saveas(f"{output_dir}/{prefix}_cs.bed")
match_utr_bed.saveas(f"{output_dir}/{prefix}_match_utr.bed")
match_nonutr_bed.saveas(f"{output_dir}/{prefix}_match_nonutr.bed")
unmatch_utr_bed.saveas(f"{output_dir}/{prefix}_unmatch_utr.bed")
unmatch_nonutr_bed.saveas(f"{output_dir}/{prefix}_unmatch_nonutr.bed")

