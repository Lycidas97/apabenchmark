import numpy as np
import pandas as pd
import pickle
import shortuuid
import random
import pandas as pd
import pysam
import logging
import subprocess
import os
import re
import hashlib

from scipy.stats import chi2_contingency
from typing import List
from array import array
from tqdm import tqdm
from itertools import chain
from Bio import SeqIO


class PAS:
    def __init__(self, chr, start, end, strand, read_df):
        self.chr = chr
        self.start = start
        self.end = end
        self.strand = strand
        self.read_df = read_df
        self.read_num = read_df.shape[0]

    def __repr__(self):
        return f"{self.chr}:{self.start}-{self.end}({self.strand})"


class TE:
    def __init__(self, te_id, pas_list, apa_type=None):
        self.te_id = te_id
        self.pas_list = pas_list
        self.pas_num = len(pas_list)
        valid_apa_types = ["diff", "non_diff", "single_pas"]
        if len(pas_list) == 1:
            self.apa_type = "single_pas"
        elif apa_type in valid_apa_types:
            self.apa_type = apa_type
        else:
            raise ValueError(
                f"apa_type must be one of {valid_apa_types}, but got '{apa_type}'"
            )

    def __repr__(self):
        return f"TE {self.te_id} with {self.pas_num} PAS"


def extract_seq_by_cigar(cigar, start_site, seq, cigar_ops, soft_clipping_fill="A"):
    sequence_parts = []
    index = start_site
    append_sequence = sequence_parts.append

    cigar_operations = cigar_ops.findall(cigar)
    for num, op in cigar_operations:
        num = int(num)
        if op == "M":
            append_sequence(seq[index : index + num])
            index += num
        elif op == "I":
            append_sequence("N" * num)
        elif op == "D" or op == "N":
            index += num
        elif op == "S":
            append_sequence(soft_clipping_fill * num)
    return "".join(sequence_parts)


def calculate_peak_range(read_df):
    seq_range_numbers = pd.to_numeric(
        read_df.cigar_string.str.extractall(r"(\d+)[MIDN]")[0]
    )
    seq_operations = read_df.cigar_string.str.extractall(r"(\d+)([MIDN])")[1]
    seq_range_numbers.loc[seq_operations == "I"] *= 0
    seq_range = seq_range_numbers.groupby(level=0).sum()
    start_range = read_df["distance_to_pas"] - seq_range
    end_range = read_df["distance_to_pas"]
    return start_range.min(), end_range.max()


def adjust_counts_to_match_total(c: np.ndarray, c_t: np.ndarray) -> np.ndarray:
    """
    Adjusts read counts among cells/spots to match the total counts.

    Parameters:
    - c (np.ndarray): Original counts for the PAS(s).
    - c_t (np.ndarray): Simulated counts for the PAS(s) after allocation.

    Returns:
    - np.ndarray: Adjusted simulated counts with errors corrected.
    """

    t_error = c - c_t.sum(axis=0)
    for i, error in enumerate(t_error):
        if error > 0:
            error_num = int(error)
            offset_base = error_num // c_t.shape[0]
            offset = np.array(
                [offset_base] * (c_t.shape[0] - error_num % c_t.shape[0])
                + [offset_base + 1] * (error_num % c_t.shape[0]),
                dtype=int,
            )
            np.random.shuffle(offset)
            c_t[:, i] += offset
        elif error < 0:
            error_num = abs(int(error))
            while error_num > 0:
                non_zero_indices = np.where(c_t[:, i] > 0)[0]
                subtract_num = min(len(non_zero_indices), error_num)
                if non_zero_indices.size > 0:
                    shuffle_indices = np.random.choice(
                        non_zero_indices, subtract_num, replace=False
                    )
                    c_t[shuffle_indices, i] -= 1
                error_num -= subtract_num
    # Validate the final adjusted counts
    if not np.isclose(c_t.sum(), c.sum()):
        raise ValueError("Adjusted counts do not match the original total counts.")
    return c_t


def allocate_reads_single_pas(count, cell_num=500, dispersion=0.1):
    """
    Allocate reads across a single PAS.

    Parameters:
    - count (int): The count data for the PAS.
    - cell_num (int, optional): The number of cells to simulate. Defaults to 500.
    - dispersion (float, optional): The dispersion parameter for the negative binomial distribution. Defaults to 0.1.

    Returns:
    - np.array: The allocated reads for the PAS.
    """
    phi = np.random.uniform(0.2, 0.8)
    e_t1 = count * phi
    mu_t1 = e_t1 / cell_num
    r_t1 = mu_t1 / dispersion
    p_t1 = r_t1 / (mu_t1 + r_t1)
    e_t1_vec = np.random.negative_binomial(n=r_t1, p=p_t1, size=cell_num)
    c_t1_vec = e_t1_vec.reshape(-1, 1)
    e_t2 = count - e_t1
    mu_t2 = e_t2 / cell_num
    r_t2 = mu_t2 / dispersion
    p_t2 = r_t2 / (mu_t2 + r_t2)
    e_t2_vec = np.random.negative_binomial(n=r_t2, p=p_t2, size=cell_num)
    c_t2_vec = e_t2_vec.reshape(-1, 1)
    c_t1 = np.array(int((count * phi).item()), dtype=int)
    c_t2 = np.array(count - c_t1, dtype=int)
    e_t1_vec = adjust_counts_to_match_total(c_t1, c_t1_vec)
    e_t2_vec = adjust_counts_to_match_total(c_t2, c_t2_vec)
    return e_t1_vec, e_t2_vec


def allocate_reads_multi_pas(
    counts,
    alpha,
    comparison_type,
    cell_num=500,
    dispersion=0.1,
    chi2_p_value=0.01,
    percentage_diff_threshold=0.2,
    swap_prob=0.1,
    pas_num=None,
):
    """
    Allocate reads across multiple PAS based on specified criteria.

    Parameters:
    - counts (np.array): The count data for PAS.
    - alpha (float): Parameter alpha for the beta distribution.
    - comparison_type (str): The type of comparison, 'diff' for differential, 'non_diff' for non-differential.
    - cell_num (int, optional): The number of cells to simulate. Defaults to 500.
    - dispersion (float, optional): The dispersion parameter for the negative binomial distribution. Defaults to 0.1.
    - chi2_p_value (float, optional): P-value threshold for chi-square test. Defaults to 0.05.
    - percentage_diff_threshold (float, optional): Threshold for percentage difference to consider. Defaults to 0.1.
    - pas_num (int, optional): The number of PAS. If not specified, inferred from counts.

    Returns:
    - (np.array, np.array): Tuple of two arrays representing allocated reads for two conditions.
    """
    total_reads = counts.sum()
    phi = np.random.uniform(0.2, 0.8)
    diff_flag = True
    if pas_num is None:
        pas_num = len(counts)
    times = 0
    while diff_flag:
        # Generate delta values using the beta distribution
        delta = np.sort(np.random.beta(alpha, alpha * phi, pas_num))
        # Randomly swap adjacent elements with 10% probability
        for i in range(len(delta) - 1):
            if np.random.rand() < swap_prob:
                delta[i], delta[i + 1] = delta[i + 1], delta[i]

        # Calculate theta for both conditions
        theta_t1 = (counts * delta) / (counts * delta).sum()
        theta_t2 = (counts * (1 - delta)) / (counts * (1 - delta)).sum()

        # Simulate reads allocation for condition 1
        e_t1 = (counts * delta).sum()
        # m_t1_nonzero = int(min(cell_num, max(int(e_t1), 1)))
        mu_t1 = e_t1 / cell_num
        r_t1 = mu_t1 / dispersion
        p_t1 = r_t1 / (mu_t1 + r_t1)
        e_t1_vec = np.random.negative_binomial(n=r_t1, p=p_t1, size=cell_num)
        c_t1_vec = np.array([np.random.multinomial(n, theta_t1) for n in e_t1_vec])

        # Simulate reads allocation for condition 2
        e_t2 = total_reads - e_t1
        # m_t2_nonzero = int(min(cell_num, max(int(e_t2), 1)))
        mu_t2 = e_t2 / cell_num
        r_t2 = mu_t2 / dispersion
        p_t2 = r_t2 / (mu_t2 + r_t2)
        e_t2_vec = np.random.negative_binomial(n=r_t2, p=p_t2, size=cell_num)
        c_t2_vec = np.array([np.random.multinomial(n, theta_t2) for n in e_t2_vec])

        # Adjust for errors in allocation
        c_t1 = (counts * delta).astype(int)
        c_t2 = counts - c_t1
        c_t1_vec = adjust_counts_to_match_total(c_t1, c_t1_vec)
        c_t2_vec = adjust_counts_to_match_total(c_t2, c_t2_vec)
        # Perform chi-square test and evaluate conditions
        chi2_data = [c_t1_vec.sum(axis=0), c_t2_vec.sum(axis=0)]
        chi2, p, dof, expected = chi2_contingency(chi2_data)
        percentage_diff = max(abs(theta_t1 - theta_t2))

        if comparison_type == "diff":
            if p < chi2_p_value and percentage_diff > percentage_diff_threshold:
                diff_flag = False
        elif comparison_type == "non_diff":
            if p > chi2_p_value and percentage_diff < percentage_diff_threshold:
                diff_flag = False
        else:
            raise ValueError("comparison_type must be 'diff' or 'non_diff'")
        times += 1

    return c_t1_vec, c_t2_vec


"""
# test
c = np.array([20, 200, 25, 30, 40])  # np.random.randint(20, 100000, 5)
alpha = 10
cell_num = 100
dispersion = 0.1
c_t1_vec, c_t2_vec = allocate_reads_multi_pas(c, alpha, "diff", cell_num, dispersion)
print(c_t1_vec.sum(axis=0))
print(c_t2_vec.sum(axis=0))
a, b = allocate_reads_single_pas(np.array([100]), 100)
"""


def reverse_cigar(cigar_str):
    elements = re.findall(r"(\d+)([MIDNSHP=X])", cigar_str)
    return "".join(num + op for num, op in elements[::-1])


def extract_seq_by_cigar(cigar, start_site, seq, cigar_ops, soft_clipping_fill="A"):
    sequence_parts = []
    index = start_site
    append_sequence = sequence_parts.append

    cigar_operations = cigar_ops.findall(cigar)
    for num, op in cigar_operations:
        num = int(num)
        if op == "M":
            append_sequence(seq[index : index + num])
            index += num
        elif op == "I":
            append_sequence("N" * num)
        elif op == "D" or op == "N":
            index += num
        elif op == "S":
            append_sequence(soft_clipping_fill * num)
    return "".join(sequence_parts)


def generate_read_list(
    te_list, genome_dict, header, cigar_ops, onsite_flag, logger, console_out=False
):
    read_list = []
    if console_out:
        te_list = tqdm(te_list)
    for te in te_list:
        for pas in te.pas_list:
            current_pas = pas
            read_df = current_pas.read_df.copy()
            sam_flag = 0 if current_pas.strand == "+" else 16
            # reverse cigar
            if current_pas.strand != read_df.strand.iloc[0]:
                read_df.loc[:, "cigar_string"] = read_df.cigar_string.apply(
                    reverse_cigar
                )

            # get read range
            seq_range_numbers = pd.to_numeric(
                read_df.cigar_string.str.extractall(r"(\d+)[MIDN]")[0]
            )
            seq_operations = read_df.cigar_string.str.extractall(r"(\d+)([MIDN])")[1]
            seq_range_numbers.loc[seq_operations == "I"] *= 0
            seq_range = seq_range_numbers.groupby(level=0).sum()

            read_df["start"] = np.where(
                current_pas.strand == "+",
                current_pas.start + read_df["distance_to_pas"] - seq_range,
                current_pas.end - read_df["distance_to_pas"],
            ).astype(int)

            read_df["end"] = np.where(
                current_pas.strand == "+",
                current_pas.start + read_df["distance_to_pas"],
                current_pas.end - read_df["distance_to_pas"] + seq_range,
            ).astype(int)

            current_pas.range_start = read_df["start"].min()
            current_pas.range_end = read_df["end"].max() + 1
            if (
                current_pas.range_start
                < 0 | current_pas.range_end
                > len(genome_dict[current_pas.chr].seq)
            ):
                logger.error(f"start or end out of range: {current_pas}")
                continue
            current_pas.seq = str(
                genome_dict[current_pas.chr].seq[
                    current_pas.range_start : current_pas.range_end
                ]
            )
            read_df.loc[:, "start_site"] = read_df["start"] - current_pas.range_start
            read_df.loc[:, "end_site"] = read_df["end"] - current_pas.range_start

            # extract_seq
            cell_barcode_list = list(current_pas.cell_barcode_generator)
            random.shuffle(cell_barcode_list)
            chr_id = header.references.index(current_pas.chr)
            for i, barcode in enumerate(cell_barcode_list):
                row = read_df.iloc[i]
                if row.is_tail & onsite_flag:
                    soft_clipping_fill = "A" if current_pas.strand == "+" else "T"
                else:
                    soft_clipping_fill = "N"
                read = pysam.AlignedSegment(header=header)
                read_seq = extract_seq_by_cigar(
                    row.cigar_string,
                    row.start_site,
                    current_pas.seq,
                    cigar_ops,
                    soft_clipping_fill,
                )
                read.query_sequence = read_seq
                read.query_name = ":".join(
                    [
                        current_pas.chr,
                        str(current_pas.start),
                        current_pas.strand,
                        str(i),
                    ]
                )
                read.flag = sam_flag
                read.reference_id = chr_id
                read.reference_start = row.start
                read.mapping_quality = 255
                read.cigarstring = row.cigar_string
                read.query_qualities = array("B", [30] * len(read_seq))
                read.set_tag("CB", barcode)
                ub_string = hashlib.md5(
                    "_".join(
                        [
                            barcode,
                            current_pas.chr,
                            str(current_pas.start),
                            current_pas.strand,
                            str(i),
                        ]
                    ).encode()
                ).hexdigest()

                read.set_tag(
                    "UB",
                    ub_string,
                )
                read_list.append(read.to_string())
    return read_list


def generate_sim_bam(
    genome_fasta_path: str,
    pas_df: pd.DataFrame,
    peaks: List[pd.DataFrame],
    output_bam_path: str,
    min_reads_per_peak: int = 20,
    cell_num_per_type: int = 100,
    diff_alpha: float = 10,
    non_diff_alpha: float = 100,
    dispersion: float = 0.1,
    barcode_length: int = 6,
    progress_flag: bool = False,
):
    """
    Generate a simulated BAM file with reads assigned to different PAS (poly(A) sites) and cell barcodes.

    Parameters:
    - genome_fasta_path (str): Path to the genome FASTA file.
    - pas_sample (pd.DataFrame): DataFrame containing PAS information. It should include the following columns:
      - 'chr' (str): The chromosome where the PAS is located.
      - 'start' (int): The start position of the PAS.
      - 'end' (int): The end position of the PAS.
      - 'strand' (str): The strand of the PAS, either "+" or "-".
      - 'exon_id' (str): The exon ID to which the PAS belongs. This is used to group PAS into corresponding transcripts or genes.
      - 'pas_type' (str): The type of the PAS, which can be one of the following values:
        - "diff_apa": Indicates that the PAS belongs to a differential APA event.
        - "non_diff_apa": Indicates that the PAS belongs to a non-differential APA event.
        - "single_pas": Indicates that the PAS is part of a single PAS gene.
    - peaks (List[pd.DataFrame]): A list of DataFrames containing retained peaks extracted from real data.
    - output_bam_path (str): Path to save the generated BAM file.
    - min_reads_per_peak (int, optional): Minimum number of reads per peak. Default is 20.
    - cell_num_per_type (int, optional): Number of cells per cell type. Default is 100.
    - diff_alpha (float, optional): Alpha parameter for differential TEs. Default is 10.
    - non_diff_alpha (float, optional): Alpha parameter for non-differential TEs. Default is 100.
    - dispersion (float, optional): Dispersion parameter for read allocation. Default is 0.1.
    - barcode_length (int, optional): Length of cell barcodes. Default is 6.
    - progress_flag (bool, optional): Flag to display progress. Default is False.

    Returns:
    None

    Description:
    The function generates a simulated BAM file with reads assigned to different PAS and cell barcodes.
    """
    genome_dict = SeqIO.to_dict(list(SeqIO.parse(genome_fasta_path, "fasta")))
    pas_df.loc[:, "chr_end"] = pas_df.apply(
        lambda x: len(genome_dict[x.chr].seq), axis=1
    )
    pas_df["end_gap"] = pas_df["chr_end"] - pas_df["end"]

    pas_df.loc[:, "chr_end"] = pas_df.apply(
        lambda x: len(genome_dict[x.chr].seq), axis=1
    )
    # set logger
    logger = logging.getLogger(__name__ + ".generate_sim_bam")
    logger.setLevel(logging.INFO)
    logger.basicConfig = False
    log_file = f"{output_bam_path}.log"

    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # process pas_df
    diff_te_list = pas_df[pas_df["pas_type"] == "diff_apa"]["exon_id"].unique()
    non_diff_te_list = pas_df[pas_df["pas_type"] == "non_diff_apa"]["exon_id"].unique()
    single_pas_gene_list = pas_df[pas_df["pas_type"] == "single_pas"][
        "exon_id"
    ].unique()
    pas_num = pas_df.shape[0]

    # select peak
    # keep more peaks than needed to provide enougt choices for TE, prevent lack of eligible peaks
    retained_peaks = [x for x in peaks if len(x) >= min_reads_per_peak]
    if len(retained_peaks) >= int(pas_num * 1.1):
        retained_peaks = random.sample(retained_peaks, int(pas_num * 1.1))
        random.shuffle(retained_peaks)
    else:
        # if not enough peaks, duplicate some peaks
        retained_peaks_additional = [
            random.choice(retained_peaks)
            for _ in range(int(pas_num * 1.1) - len(retained_peaks))
        ]
        retained_peaks += retained_peaks_additional
        random.shuffle(retained_peaks)

    logger.info(f"Number of PAS: {pas_num}")
    logger.info(f"Number of retained peaks: {len(retained_peaks)}")

    t1_cell_barcode_list = [
        shortuuid.ShortUUID().random(length=barcode_length)
        for _ in range(cell_num_per_type)
    ]
    t2_cell_barcode_list = [
        shortuuid.ShortUUID().random(length=barcode_length)
        for _ in range(cell_num_per_type)
    ]

    peak_read_num = [len(x) for x in retained_peaks]
    peak_range = [calculate_peak_range(x) for x in retained_peaks]
    remaining_set = set(zip(peak_read_num, range(len(peak_read_num)), peak_range))
    original_remaining_set = remaining_set.copy()

    logger.info(f"Number of remaining peaks: {len(remaining_set)}")
    te_group = []

    # set alpha=10 to generate differential TEs
    # alpha controls the dispersion of the beta distribution
    # smaller alpha leads to larger dispersion
    alpha = 10
    for te in diff_te_list:
        pas_group = pas_df[pas_df["exon_id"] == te]
        pas_num = len(pas_group)
        pas_start_gap = pas_group["start"]
        pas_end_gap = pas_group["end_gap"]
        pas_strand = pas_group["strand"].iloc[0]
        if pas_strand == "+":
            pass
        else:
            pas_start_gap, pas_end_gap = pas_end_gap, pas_start_gap
        while True:
            current_group = set(random.sample(list(remaining_set), pas_num))
            # check if the peaks are within the range of the chromsome
            start_positions = pas_start_gap + np.array([x[2][0] for x in current_group])
            end_positions = pas_end_gap - np.array([x[2][1] for x in current_group])
            if (start_positions < 0).any() or (end_positions < 0).any():
                continue

            values_sorted = sorted(current_group, key=lambda x: x[0], reverse=True)
            max_value_in_group = values_sorted[0][0]
            second_max_value_in_group = values_sorted[1][0]
            if max_value_in_group <= 2 * second_max_value_in_group:
                remaining_set -= current_group
                break

            eligible_second_peaks = [
                x
                for x in remaining_set
                if max_value_in_group * 0.5 <= x[0] <= max_value_in_group * 2
            ]
            eligible_second_peaks = [
                x for x in eligible_second_peaks if x not in current_group
            ]

            if not eligible_second_peaks:
                eligible_second_peaks = [
                    x
                    for x in original_remaining_set
                    if max_value_in_group * 0.5 <= x[0] <= max_value_in_group * 2
                ]
                eligible_second_peaks = [
                    x for x in eligible_second_peaks if x not in current_group
                ]

            if eligible_second_peaks:
                new_second_peak = random.sample(eligible_second_peaks, 1)
                old_second_max_value_peak = values_sorted[1]
                current_group.discard(old_second_max_value_peak)
                current_group |= set(new_second_peak)

                # Check if the updated peaks are within the range of the chromosome
                start_positions = pas_start_gap + np.array(
                    [x[2][0] for x in current_group]
                )
                end_positions = pas_end_gap - np.array([x[2][1] for x in current_group])
                if (start_positions < 0).any() or (end_positions < 0).any():
                    continue
                remaining_set -= current_group
                break
            else:
                remaining_set |= current_group
                remaining_set.remove(max(current_group, key=lambda x: x[0]))

        current_peak_list = [retained_peaks[x] for _, x, i in current_group]
        pas_list = [
            PAS(x.chr, x.start, x.end, x.strand, current_peak_list[_])
            for _, x in pas_group.sort_values(["chr", "start"]).reset_index().iterrows()
        ]
        pas_counts = np.array([x.read_num for x in pas_list])
        # generate read distribution
        c_t1_vec, c_t2_vec = allocate_reads_multi_pas(
            pas_counts, alpha, "diff", cell_num_per_type, dispersion
        )

        # allocate reads to cell barcode
        for i, pas in enumerate(pas_list):
            pas.cell_barcode_generator = (
                barcode
                for barcode in chain(
                    np.repeat(t1_cell_barcode_list, c_t1_vec[:, i]),
                    np.repeat(t2_cell_barcode_list, c_t2_vec[:, i]),
                )
            )
        te_group.append(TE(te, pas_list, "diff"))
        te_group[-1].c_t1_vec = c_t1_vec
        te_group[-1].c_t2_vec = c_t2_vec

    logger.info(f"Number of remaining peaks: {len(remaining_set)}")

    # set alpha=100 to generate non-differential TEs
    alpha = 100
    for te in non_diff_te_list:
        pas_group = pas_df[pas_df["exon_id"] == te]
        pas_num = len(pas_group)
        pas_start_gap = pas_group["start"]
        pas_end_gap = pas_group["end_gap"]
        pas_strand = pas_group["strand"].iloc[0]
        if pas_strand == "+":
            pass
        else:
            pas_start_gap, pas_end_gap = pas_end_gap, pas_start_gap
        while True:
            current_group = set(random.sample(list(remaining_set), pas_num))

            # check if the peaks are within the range of the chromsome
            start_positions = pas_start_gap + np.array([x[2][0] for x in current_group])
            end_positions = pas_end_gap - np.array([x[2][1] for x in current_group])
            if (start_positions < 0).any() or (end_positions < 0).any():
                continue

            values_sorted = sorted(current_group, key=lambda x: x[0], reverse=True)
            max_value_in_group = values_sorted[0][0]
            second_max_value_in_group = values_sorted[1][0]
            if max_value_in_group <= 2 * second_max_value_in_group:
                remaining_set -= current_group
                break

            eligible_second_peaks = [
                x
                for x in remaining_set
                if max_value_in_group * 0.5 <= x[0] <= max_value_in_group * 2
            ]
            eligible_second_peaks = [
                x for x in eligible_second_peaks if x not in current_group
            ]

            if not eligible_second_peaks:
                eligible_second_peaks = [
                    x
                    for x in original_remaining_set
                    if max_value_in_group * 0.5 <= x[0] <= max_value_in_group * 2
                ]
                eligible_second_peaks = [
                    x for x in eligible_second_peaks if x not in current_group
                ]

            if eligible_second_peaks:
                new_second_peak = random.sample(eligible_second_peaks, 1)
                old_second_max_value_peak = values_sorted[1]
                current_group.discard(old_second_max_value_peak)
                current_group |= set(new_second_peak)

                # Check if the updated peaks are within the range of the chromosome
                start_positions = pas_start_gap + np.array(
                    [x[2][0] for x in current_group]
                )
                end_positions = pas_end_gap - np.array([x[2][1] for x in current_group])
                if (start_positions < 0).any() or (end_positions < 0).any():
                    continue

                remaining_set -= current_group
                break
            else:
                remaining_set |= current_group
                remaining_set.remove(max(current_group, key=lambda x: x[0]))

        current_peak_list = [retained_peaks[x] for _, x, i in current_group]
        pas_list = [
            PAS(x.chr, x.start, x.end, x.strand, current_peak_list[_])
            for _, x in pas_group.sort_values(["chr", "start"]).reset_index().iterrows()
        ]
        pas_counts = np.array([x.read_num for x in pas_list])

        # generate read distribution
        c_t1_vec, c_t2_vec = allocate_reads_multi_pas(
            pas_counts, alpha, "non_diff", cell_num_per_type, dispersion
        )

        # allocate reads to cell barcode
        for i, pas in enumerate(pas_list):
            pas.cell_barcode_generator = (
                barcode
                for barcode in chain(
                    np.repeat(t1_cell_barcode_list, c_t1_vec[:, i]),
                    np.repeat(t2_cell_barcode_list, c_t2_vec[:, i]),
                )
            )
        te_group.append(TE(te, pas_list, "non_diff"))
        te_group[-1].c_t1_vec = c_t1_vec
        te_group[-1].c_t2_vec = c_t2_vec

    logger.info(f"Number of remaining peaks: {len(remaining_set)}")

    for te in single_pas_gene_list:
        pas_group = pas_df[pas_df["exon_id"] == te]
        pas_num = len(pas_group)

        if pas_num != 1:
            raise ValueError("pas_num must be 1 for single_pas_gene_list")

        pas_start_gap = pas_group["start"]
        pas_end_gap = pas_group["end_gap"]
        pas_strand = pas_group["strand"].iloc[0]
        if pas_strand == "+":
            pass
        else:
            pas_start_gap, pas_end_gap = pas_end_gap, pas_start_gap

        while True:
            current_group = set(random.sample(list(remaining_set), 1))
            # check if the peaks are within the range of the chromsome
            start_positions = pas_start_gap + np.array([x[2][0] for x in current_group])
            end_positions = pas_end_gap - np.array([x[2][1] for x in current_group])
            if (start_positions < 0).any() or (end_positions < 0).any():
                continue
            remaining_set -= current_group
            break

        current_peak_list = [retained_peaks[x] for _, x, i in current_group]
        pas_list = [
            PAS(x.chr, x.start, x.end, x.strand, current_peak_list[_])
            for _, x in pas_group.sort_values(["chr", "start"]).reset_index().iterrows()
        ]
        pas_counts = np.array([x.read_num for x in pas_list])

        # generate read distribution
        c_t1_vec, c_t2_vec = allocate_reads_single_pas(
            pas_counts, cell_num_per_type, dispersion
        )

        # allocate reads to cell barcode
        for i, pas in enumerate(pas_list):
            pas.cell_barcode_generator = (
                barcode
                for barcode in chain(
                    np.repeat(t1_cell_barcode_list, c_t1_vec[:, i]),
                    np.repeat(t2_cell_barcode_list, c_t2_vec[:, i]),
                )
            )

        te_group.append(TE(te, pas_list, "single_pas"))
        te_group[-1].c_t1_vec = c_t1_vec
        te_group[-1].c_t2_vec = c_t2_vec

    logger.info(f"Number of remaining peaks: {len(remaining_set)}")
    # generate read list
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.0"},
            "SQ": [{"LN": len(x.seq), "SN": x.id} for x in genome_dict.values()],
        }
    )
    cigar_ops = re.compile(r"(\d+)([MIDNSHP=X])")
    onsite_flag = True
    read_list = generate_read_list(
        te_group, genome_dict, header, cigar_ops, onsite_flag, logger, progress_flag
    )

    logger.info(f"Number of reads generated: {len(read_list)}")

    # write to bam
    bam = pysam.AlignmentFile(output_bam_path, "wb", header=header)
    for read in read_list:
        bam.write(pysam.AlignedSegment.fromstring(read, header=bam.header))
    bam.close()

    logger.info(f"Simulated bam file written to {output_bam_path}")

    result = subprocess.run(
        ["samtools", "sort", "-o", f"{output_bam_path}.sorted", output_bam_path]
    )
    if result.returncode != 0:
        logger.error(
            f"Error occurred while running 'samtools sort'. Return code: {result.returncode}"
        )
        raise RuntimeError("Error occurred while running 'samtools sort'")

    os.remove(output_bam_path)
    os.rename(f"{output_bam_path}.sorted", output_bam_path)
    result = subprocess.run(["samtools", "index", output_bam_path])
    if result.returncode != 0:
        logger.error(
            f"Error occurred while running 'samtools index'. Return code: {result.returncode}"
        )
        raise RuntimeError("Error occurred while running 'samtools index'")

    # save expr mtx
    pas_name_list = list(
        chain.from_iterable(
            [
                [
                    "|".join(
                        [x.te_id, ":".join([y.chr, str(y.start), str(y.end), y.strand])]
                    )
                    for y in x.pas_list
                ]
                for x in te_group
            ]
        )
    )
    cell_list = t1_cell_barcode_list + t2_cell_barcode_list
    pas_expr_mtx = np.hstack([np.vstack([x.c_t1_vec, x.c_t2_vec]) for x in te_group])
    pas_expr_df = pd.DataFrame(pas_expr_mtx, index=cell_list, columns=pas_name_list)
    cell_type_df = pd.DataFrame(
        np.array(["T1"] * cell_num_per_type + ["T2"] * cell_num_per_type),
        index=cell_list,
        columns=["cell_type"],
    )
    pas_expr_df = pd.concat([cell_type_df, pas_expr_df], axis=1)
    pas_expr_df.to_csv(f"{output_bam_path}.expr.tsv", sep="\t", index=True)
    logger.info(f"Expression matrix written to {output_bam_path}.expr.tsv")

    file_handler.close()
    logger.removeHandler(file_handler)
