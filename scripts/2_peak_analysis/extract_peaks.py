import sys

import pandas as pd
import pyarrow.feather as feather
from argparse import ArgumentParser

sys.path.append("../apasim/")
from extraction import extract_peaks

# multiprocessing.set_start_method('forkserver', force=True)
import os
# os.environ['OPENBLAS_NUM_THREADS'] = '1'

parser = ArgumentParser()
parser.add_argument("-b", "--bam", dest="bam", help="path to config file", type=str)
parser.add_argument("-p", "--pas", dest="pas", help="pas", type=str)
parser.add_argument("-o", "--output", dest="output", help="output", type=str,)
parser.add_argument("-j", "--n_jobs", dest="n_jobs", help="number of jobs", type=int, default=8)
parser.add_argument("-u", "--upstream", dest="upstream", help="upstream", type=int, default=400)
parser.add_argument("-d", "--downstream", dest="downstream", help="downstream", type=int, default=100)

args = parser.parse_args()
bam_path = args.bam
pas_path = args.pas
output_path = args.output
n_jobs = args.n_jobs
upstream = args.upstream
downstream = args.downstream

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
    df['start'] = pd.to_numeric(df['start'])
    df['end'] = pd.to_numeric(df['end'])

    if not all(df['start'] < df['end']):
        raise ValueError("'start' must be less than 'end'.")

    return df.sort_values(["chr","start","end"])


if __name__ == "__main__":
    # pas_annotaion = pd.read_csv(pas_path, sep="\t", header=None)
    # pas_annotaion.columns = ["chr", "start", "end", "name", "score", "strand","pas_type", "gene_id", "gene_symbol", "TE_id"]
    pas_annotaion = read_bed_like_file(pas_path)
    pas_annotaion = pas_annotaion[pas_annotaion["chr"].str.match("chr[0-9XY]+$")]

    peak_list = extract_peaks(
        bam_path, pas_annotaion, upstream, downstream, num_processors=n_jobs
    )
    peak_list = [x for x in peak_list if len(x) > 0]
    peak_df = pd.concat(peak_list)
    feather.write_feather(peak_df, output_path)
    # with open(output_path, "wb") as f:
    #     pickle.dump(peak_df, f)