import sys
import pandas as pd
from argparse import ArgumentParser
import os
import pickle
import numpy as np

def generate_junction_cigar(
    row,
    center,
    read_length
):
    if row["distance_to_pas"] <= center:
        return f"{read_length}M"
    elif row["distance_to_pas"] - read_length >= center:
        return f"{read_length}M"
    else:
        return f"{center - row['distance_to_pas'] + read_length}M1000N{row['distance_to_pas']-center}M"

def generate_sim_peak(
    peak_size,
    read_length,
    center=-50,
    std_dev=10,
    junction_prob=0.2,
    junction_fraction=0.1,
):
    distance_list = np.random.normal(center, std_dev, peak_size)
    distance_list = np.round(distance_list).astype(int)
    distance_list = np.clip(distance_list, -10000, 0)
    peak_df = pd.DataFrame({"distance_to_pas": distance_list})
    peak_df["cigar_string"] = f"{read_length}M"
    if np.random.rand() < junction_prob:
        peak_df["cigar_string"] = peak_df.apply(
            generate_junction_cigar, axis=1, center=center, read_length=read_length
        )

    peak_df["reference_length"] = read_length
    peak_df["strand"] = "+"
    peak_df["is_tail"] = False
    peak_df["is_gap"] = False
    return peak_df


peak_num = 500
peak_sizes = [250, 500, 1000, 2000, 4000]
read_lengths = [40, 70, 100, 130, 160]
center = -50
std_dev = 10

for peak_size in peak_sizes:
    for read_length in read_lengths:
        peak_list = [generate_sim_peak(peak_size, read_length) for i in range(peak_num)]
        with open(f"../../data/sim_data/sim_peak/simulated_peak_ps{peak_size}_rl{read_length}.pickle", "wb") as f:
            pickle.dump(peak_list, f)