import numpy as np
import pandas as pd
import os
import pickle
from pybedtools import BedTool
import pybedtools
from argparse import ArgumentParser
import warnings

def read_bed_like_file(file_path):
    """
    Read a BED-like file and return a pandas DataFrame.

    Parameters:
    file_path (str): The path to the input file.

    Returns:
    pandas.DataFrame: The DataFrame containing the data from the input file.

    Raises:
    ValueError: If the input file has less than 6 columns or if 'start' is not less than 'end'.
    """

    df = pd.read_csv(file_path, sep="\t", header=None, dtype=str)

    if len(df.columns) < 6:
        raise ValueError("Input file must have at least 6 columns.")

    if pd.to_numeric(df.iloc[0, 1], errors='coerce') is not None and pd.to_numeric(df.iloc[0, 2], errors='coerce') is not None:
        df = df.iloc[1:]

    if not all(df.iloc[:, 5].isin(['+', '-'])):
        df = df.iloc[1:]

    df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand'] + [str(x) for x in df.columns[6:]]
    df['start'] = pd.to_numeric(df['start']).astype(int)
    df['end'] = pd.to_numeric(df['end']).astype(int)

    if not all(df['start'] < df['end']):
        raise ValueError("'start' must be less than 'end'.")

    return df.sort_values(["chr","start","end"])

def rename_matched_columns(match_df, gt_col_num, res_col_num):
    """
    Renames the columns of the given DataFrame based on the number of columns in the ground truth (gt) and result (res).

    Args:
        match_df (pandas.DataFrame): The DataFrame containing the matched data.
        gt_col_num (int): The number of columns in the ground truth.
        res_col_num (int): The number of columns in the result.

    Returns:
        pandas.DataFrame: The DataFrame with renamed columns.

    """
    column_names = [
        'gt_chr', 'gt_start', 'gt_end', 'gt_name', 'gt_score', 'gt_strand',
        'res_chr', 'res_start', 'res_end', 'res_name', 'res_score', 'res_strand'
    ]

    if gt_col_num > 6:
        for i in range(6, gt_col_num):
            column_names.insert(i, f'gt_{i}')

    if res_col_num > 6:
        for i in range(gt_col_num + 6, gt_col_num + res_col_num):
            column_names.append(f'res_{i}')
            
    match_df.columns = column_names

    return match_df

def deduplicate_matches(match_df):
    """
    Deduplicates matches in the given DataFrame based on the minimum distance between 'gt_start' and 'res_start' columns.

    Args:
        match_df (DataFrame): The DataFrame containing the matches.

    Returns:
        DataFrame: The deduplicated DataFrame with an additional 'best_match' column indicating the best match for each unique combination of 'gt_chr', 'gt_start', 'gt_end', and 'gt_strand'.
    """
    match_df['distance'] = np.abs(match_df['gt_start'] - match_df['res_start'])
    drop_duplicated_df = match_df.loc[match_df.groupby("res_name")['distance'].idxmin()]
    gt_unique_cols = ['gt_chr', 'gt_start', 'gt_end', 'gt_strand']
    drop_duplicated_df.loc[drop_duplicated_df.groupby(gt_unique_cols)['distance'].idxmin()]

    gt_unique_cols = ['gt_chr', 'gt_start', 'gt_end', 'gt_strand']
    drop_duplicated_df["best_match"] = False
    drop_duplicated_df.loc[drop_duplicated_df.groupby(gt_unique_cols)['distance'].idxmin(), "best_match"] = True

    return drop_duplicated_df.sort_values(["gt_chr", "gt_start", "gt_end", "gt_strand"]).drop(columns=["distance"])

def generate_pas_pairs(group):
    sorted_group = group.sort_values('start')
    pas_pairs = []
    for i in range(len(sorted_group) - 1):
        pas1 = sorted_group.iloc[i]['gt_pas']
        pas2 = sorted_group.iloc[i+1]['gt_pas']
        gap = abs(sorted_group.iloc[i+1]['start'] - sorted_group.iloc[i]['start'])
        min_counts = min(sorted_group.iloc[i]['counts'], sorted_group.iloc[i+1]['counts'])
        pas_pairs.append((pas1, pas2, gap, min_counts))
    return pd.DataFrame(pas_pairs, columns=['pas1', 'pas2', 'gap', 'min_counts'])


parser = ArgumentParser()
parser.add_argument("--tool", dest="tool", help="tool name", required=True)
parser.add_argument("--sample", dest="sample", help="sample id", required=True)
parser.add_argument("--gt_pas", dest="gt_pas", help="ground truth PAS file", required=True)
parser.add_argument("--pd_pas", dest="pd_pas", help="predicted PAS file", required=True)
parser.add_argument("--gt_mtx", dest="gt_mtx", help="ground truth matrix file", required=True)
parser.add_argument("--output_dir", dest="output_dir", help="output directory", required=True)
parser.add_argument("--genome_size", dest="genome_size", help="genome size file", required=True)

args = parser.parse_args()
sample = args.sample
tool = args.tool
gt_pas_path = args.gt_pas
pd_pas_path = args.pd_pas
gt_mtx_path = args.gt_mtx
output_dir = args.output_dir
genome_size_path = args.genome_size

# read ground truth
gt_pas_df = read_bed_like_file(gt_pas_path)
# set te name and te type
gt_pas_df = gt_pas_df.rename(columns={"7": "te_name", "6": "te_type"})
gt_pas_df["gt_pas"] = gt_pas_df.apply(lambda row: f"{row['te_name']}|{row['chr']}:{row['start']}:{row['end']}:{row['strand']}", axis=1)
gt_pas_df["pas_num"] = gt_pas_df.groupby("te_name")["chr"].transform("count")
# grep te
gt_te_df = gt_pas_df.copy()
gt_te_df.drop(columns=["gt_pas"], inplace=True)
gt_te_df["pas_num"] = gt_te_df.groupby("te_name")["chr"].transform("count")
gt_te_df = gt_te_df.drop_duplicates(["te_name"])
gt_te_df["chr"] = gt_te_df["te_name"].str.split(":", expand=True)[2]
gt_te_df["start"] = gt_te_df["te_name"].str.split(":", expand=True)[3].astype(int)
gt_te_df["end"] = gt_te_df["te_name"].str.split(":", expand=True)[4].astype(int)
gt_te_df["strand"] = gt_te_df["te_name"].str.split(":", expand=True)[5]


gt_pas_bed = BedTool.from_dataframe(gt_pas_df.sort_values(["chr","start","end"]))

# annotate pa site count
gt_counts_df = pd.read_csv(gt_mtx_path, sep="\t", index_col=0)
gt_counts_df = gt_counts_df.sum(axis=0).drop("cell_type").rename("gt_counts")
gt_counts_df = pd.DataFrame(gt_counts_df).reset_index(names=["pas_name"])
gt_pas_df["counts"] = gt_pas_df["gt_pas"].map(gt_counts_df.set_index("pas_name")["gt_counts"])
gt_pas_df = gt_pas_df[gt_pas_df["pas_num"] > 1]

pas_pairs = gt_pas_df.groupby('te_name').apply(generate_pas_pairs)
result_te_df = gt_te_df[gt_te_df["pas_num"] > 1]
pas_pairs.reset_index().drop(columns=["level_1"]).groupby("te_name").apply(lambda x: max(x["gap"]))
result_te_df["max_gap"] = result_te_df["te_name"].map(pas_pairs.reset_index().drop(columns=["level_1"]).groupby("te_name").apply(lambda x: max(x["gap"])))
result_te_df["max_gap_min_counts"] = result_te_df["te_name"].map(pas_pairs.reset_index().drop(columns=["level_1"]).groupby("te_name").apply(lambda x: x[x["gap"] == max(x["gap"])]["min_counts"].values[0]))


# read result
pd_pas_df = read_bed_like_file(pd_pas_path)
pd_pas_df = pd_pas_df.drop_duplicates()
pd_pas_df["rank"] = pd_pas_df.groupby("name")["strand"].transform(
        lambda x: pd.Series(range(len(x)), index=x.index) / (len(x) - 1) if x.iloc[0] == "-" else pd.Series(
            range(len(x) - 1, -1, -1), index=x.index) / (len(x) - 1))
pd_pas_df = pd_pas_df.fillna(0)
pd_pas_df = pd_pas_df[pd_pas_df["rank"] == 0].drop(columns=["rank"])

predict_pas = BedTool.from_dataframe(pd_pas_df).sort()



# match at te level
gt_te_pd_map_df = BedTool.from_dataframe(gt_te_df).sort().intersect(predict_pas, wa=True, wb=True,s=True).to_dataframe(header=None)
gt_te_pd_map_df.columns = ['te_chr', 'te_start', 'te_end', 'gene_name', 'te_score', 'te_strand', 'te_type', 'te_name','pas_num','pas_chr', 'pas_start', 'pas_end', 'pas_name', 'pas_score', 'pas_strand']
gt_te_pd_map_df["pd_pas_num"] = gt_te_pd_map_df.groupby("te_name")["te_chr"].transform("count")
gt_te_pd_map_df["pd_pas"] = gt_te_pd_map_df["te_name"] + "|" + gt_te_pd_map_df["pas_name"]
gt_te_pd_map_df = gt_te_pd_map_df.sort_values(["te_chr", "te_start", "te_end", "te_strand", "pas_start", "pas_end", "pas_strand"])
gt_te_pd_map_df = gt_te_pd_map_df.drop_duplicates(["te_name"])
gt_te_pd_map_df = gt_te_pd_map_df.drop_duplicates(["te_name"])
gt_te_pd_map_df = gt_te_pd_map_df[gt_te_pd_map_df["te_type"] != "single_pas"]

result_te_df = pd.merge(gt_te_pd_map_df[["te_name", "pd_pas_num"]], result_te_df, on="te_name", how="right")

for slop in [50,100,150,200,250,300,350,400]:
    gt_te_pd_map_df = BedTool.from_dataframe(gt_te_df).sort().slop(b=slop, g=genome_size_path).intersect(predict_pas, wa=True, wb=True,s=True).to_dataframe(header=None)
    gt_te_pd_map_df.columns = ['te_chr', 'te_start', 'te_end', 'gene_name', 'te_score', 'te_strand', 'te_type', 'te_name','pas_num','pas_chr', 'pas_start', 'pas_end', 'pas_name', 'pas_score', 'pas_strand']
    gt_te_pd_map_df[f"pd_pas_num_slop_{slop}"] = gt_te_pd_map_df.groupby("te_name")["te_chr"].transform("count")
    gt_te_pd_map_df = gt_te_pd_map_df.drop_duplicates(["te_name"])
    gt_te_pd_map_df = gt_te_pd_map_df[gt_te_pd_map_df["te_type"] != "single_pas"]
    result_te_df = pd.merge(gt_te_pd_map_df[["te_name", f"pd_pas_num_slop_{slop}"]], result_te_df, on="te_name", how="right")

result_te_df = result_te_df.fillna(0)
result_te_df["tool"] = tool
result_te_df["sample"] = sample
result_te_df["protocol"] = sample.split("_")[0]

result_te_df.to_csv(f"{output_dir}/{tool}/{sample}_te_gap.tsv", sep="\t", index=False)