import pickle
import os
import pandas as pd 
import numpy as np 
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest="input", help="path to input fasta file", type=str)
parser.add_argument("-o", "--output_dir", dest="output", help="directory to output", type=str)


SEQUENCE_LENGTH=1001
pas_major_motif = motifs.create([Seq("AATAAA")])
pas_other_seq = ["ATTAAA", "AGTAAA", "TATAAA", "CATAAA", "TTTAAA","AATACA","AATATA", "GATAAA", "AAGAAA", "AATGAA", "AATAGA","ACTAAA"]
pas_other_motif = motifs.create([Seq(x) for x in pas_other_seq])
cfi_motif = motifs.create([Seq("TGTA")])
cfii_seq = [ "TGTGTG","TGTGTT","TGTTGG","TGTTGT","TTGGTG","TTGGTT","TTGTGG","TTGTGT"]
cfii_motif =  motifs.create([Seq(x) for x in cfii_seq ])

pas_search_region = np.array([-100, 0])
cfi_search_region = np.array([-100, 0])
cfii_search_region = np.array([0, 100])

pas_search_region = pas_search_region + SEQUENCE_LENGTH // 2
cfi_search_region = cfi_search_region + SEQUENCE_LENGTH // 2
cfii_search_region = cfii_search_region + SEQUENCE_LENGTH // 2

args = parser.parse_args()
input_path = args.input
output_dir = args.output

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

sample_id = "_".join(os.path.basename(input_path).split("_")[0:4])
match_type = "_".join(os.path.basename(input_path).replace("_extended","").split("_")[4:])
prefix = os.path.basename(input_path).split(".")[0].replace("_extended","")

with open(input_path) as f:
    sequences = [record.seq for record in SeqIO.parse(f, "fasta") if len(record.seq) == SEQUENCE_LENGTH]

def search_motif(seq, motif):
    motif_search_list = [pos for pos in seq.search(motif.alignment.sequences)]
    if motif_search_list:
        motif_search_result = np.array([x[0] for x in motif_search_list])
    else:
        motif_search_result = np.array([np.nan])
    return motif_search_result

motif_dict = {
    'pas_major': pas_major_motif,
    'pas_other': pas_other_motif,
    'cfi': cfi_motif,
    'cfii': cfii_motif
}

search_regions = {
    'pas_major': pas_search_region,
    'pas_other': pas_search_region,
    'cfi': cfi_search_region,
    'cfii': cfii_search_region
}

motif_search_results = {}
motif_found_vectors = {}
motif_positions = {}

for motif_name, motif in motif_dict.items():
    search_region = search_regions[motif_name]
    motif_search_results[motif_name] = [search_motif(x, motif) for x in sequences]
    motif_found_vectors[motif_name] = np.array([((search_region[0] <= x) & (x <= search_region[1])).any() for x in motif_search_results[motif_name]])
    motif_positions[motif_name] = np.concatenate(motif_search_results[motif_name])
    motif_positions[motif_name] = motif_positions[motif_name][~np.isnan(motif_positions[motif_name])]

categories = [
    "PAS Major + CFI/II",
    "PAS Major only",
    "PAS Other + CFI/II",
    "PAS Other only",
    "CFI/II only",
    "None"
]

has_pas_major = motif_found_vectors['pas_major']
has_pas_other = motif_found_vectors['pas_other']
has_cfi = motif_found_vectors['cfi']
has_cfii = motif_found_vectors['cfii']

site_categories = np.full(len(sequences), "None", dtype=object)

site_categories[(has_pas_major) & (has_cfi | has_cfii)] = categories[0]
site_categories[(has_pas_major) & ~(has_cfi | has_cfii)] = categories[1]
site_categories[has_pas_other & ~has_pas_major & (has_cfi | has_cfii)] = categories[2]
site_categories[(has_pas_other & ~has_pas_major) & ~(has_cfi | has_cfii)] = categories[3]
site_categories[(has_cfi | has_cfii) & ~(has_pas_major | has_pas_other)] = categories[4]
site_categories[~(has_pas_major | has_pas_other | has_cfi | has_cfii)] = categories[5]

category_counts = {}
for category in categories:
    count = np.sum(site_categories == category)
    category_counts[category] = count

category_df = pd.DataFrame(category_counts.items(), columns=["category", "count"])
category_df["prop"] = category_df["count"] / category_df["count"].sum()
category_df["sample_id"] = sample_id
category_df["type"] = match_type
category_df["protocol"] = sample_id.split("_")[0]

density_data = []
site_range = (0, SEQUENCE_LENGTH)
for motif_name, positions in motif_positions.items():
    counts, bins = np.histogram(positions, bins=range(site_range[0], site_range[1] + 1))
    density = counts / counts.sum()
    # 将密度、位点和 motif 名称添加到列表中
    for pos, dens in zip(bins[:-1], density):
        density_data.append([dens, pos, motif_name])

# 将密度直方图数据转换为 DataFrame
density_df = pd.DataFrame(density_data, columns=['density', 'position', 'motif'])
density_df["sample_id"] = sample_id
density_df["type"] = match_type
density_df["protocol"]= sample_id.split("_")[0]

with open(f"{output_dir}/{prefix}_motif_density.pkl", "wb") as f:
    pickle.dump(density_df, f)

with open(f"{output_dir}/{prefix}_motif_category.pkl", "wb") as f:
    pickle.dump(category_df, f)