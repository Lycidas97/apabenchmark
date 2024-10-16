import re
import pandas as pd
import numpy as np
from typing import List, Union
import pysam
from dask import delayed, compute
from dask.diagnostics import ProgressBar
from itertools import chain


def extract_peaks(
    bam_path: str,
    pas_annotaion: pd.DataFrame,
    upstream: int = 1000,
    downstream: int = 1000,
    num_processors: int = 4,
    pA_allowed_mismatches_at_start: int=2,
    pA_match_length: int=10,
    pA_min_required_matches: int=8,
    progress_bar: bool = False,
) -> pd.DataFrame:
    """
    Extract reads surrounding PolyA sites from a BAM file and compile them into a DataFrame.

    Parameters:
        bam_path (str): Path to the BAM file.
        pas_annotation (pd.DataFrame): DataFrame containing PolyA site annotations.
        upstream (int): Number of bases to consider upstream of the PolyA site.
        downstream (int): Number of bases to consider downstream of the PolyA site.
        num_processors (int): Number of processors to use for parallel processing.

    Returns:
        pd.DataFrame: A DataFrame with the compiled reads data.
    """

    peak_chunks = []
    for chr in pas_annotaion["chr"].unique():
        chr_pas_annotaion = pas_annotaion[pas_annotaion["chr"] == chr]
        chunk = delayed(extract_surrounding_reads_in_chunk)(
            bam_path,
            chr_pas_annotaion,
            upstream,
            downstream,
            pA_allowed_mismatches_at_start,
            pA_match_length,
            pA_min_required_matches
        )
        peak_chunks.append(chunk)
    if progress_bar:
        with ProgressBar():
            peak_list = compute(
                *peak_chunks, num_workers=num_processors, scheduler="processes"
            )
    else:
        peak_list = compute(
            *peak_chunks, num_workers=num_processors, scheduler="processes"
        )
    peak_list = list(chain.from_iterable(peak_list))
    return peak_list


def extract_surrounding_reads_in_chunk(
    bam_path: str,
    pas_chunk: pd.DataFrame,
    upstream: int = 1000,
    downstream: int = 1000,
    pA_allowed_mismatches_at_start: int=2,
    pA_match_length: int=10,
    pA_min_required_matches: int=8,
) -> pd.DataFrame:
    """
    Extract reads surrounding PolyA sites from a specific chromosome chunk.

    This function processes chunks of PAS annotations for a given chromosome to
    extract surrounding reads based on specified upstream and downstream distances.

    Parameters are the same as in the `extract_peaks` function, targeting a specific chunk.

    Returns:
        pd.DataFrame: A DataFrame of reads surrounding PolyA sites within the chunk.
    """
    peak_list = []
    for i, row in pas_chunk.iterrows():
        peak = extract_surrounding_reads_in_pas(
            bam_path, row, upstream, downstream, pA_allowed_mismatches_at_start, pA_match_length, pA_min_required_matches
        )
        peak_list.append(peak)
    return peak_list


def extract_surrounding_reads_in_pas(
    bam_path: str,
    row: Union[pd.Series, pd.DataFrame],
    upstream: int = 1000,
    downstream: int = 1000,
    pA_allowed_mismatches_at_start: int=2,
    pA_match_length: int=10,
    pA_min_required_matches: int=8,
) -> pd.DataFrame:
    """
    Extract reads for a single PolyA site based on its annotation row.

    This function fetches reads from the BAM file that are within the specified
    upstream and downstream range of a single PolyA site. It also filters reads
    based on the provided regex patterns for forward and reverse strands.

    Parameters are similar to those in the `extract_peaks` function, but targeted
    at a single PAS annotation row.

    Returns:
        pd.DataFrame: A DataFrame containing information on the reads surrounding the specific PolyA site.
    """

    read = row
    bam = pysam.AlignmentFile(bam_path, "rb")
    chrom = read["chr"]
    strand = read["strand"]
    if strand == "+":
        pas = read.start
        start = read.start - upstream
        end = read.start + downstream
    elif strand == "-":
        pas = read.end
        start = read.end - downstream
        end = read.end + upstream
    else:
        raise ValueError('strand must be either "+" or "-"')

    start = max(0, start)
    end = max(0, end)

    for c in bam.references:
        if c.startswith("chr"):
            break
        else:
            chrom = chrom.replace("chr", "")
            break
    reads = bam.fetch(chrom, start, end)
    if reads is None:
        return pd.DataFrame()

    data = []

    # check tail reads
    cigar_pattern_1 = re.compile(r"[0-9]*M[0-9]{2,}S")
    cigar_pattern_2 = re.compile(r"[0-9]{2,}S[0-9]*M")

    unmatch_pattern_1 = re.compile(r"(\d+)S$")
    unmatch_pattern_2 = re.compile(r"^(\d+)S")

    seqpattern_1 = re.compile(fr'^[^A]{{0,{pA_allowed_mismatches_at_start}}}A[A-Z]{{{pA_match_length-1}}}')
    seqpattern_2 = re.compile(fr'^[^T]{{0,{pA_allowed_mismatches_at_start}}}T[A-Z]{{{pA_match_length-1}}}')

    for read in reads:
        read_strand = "+" if not read.is_reverse else "-"
        if read_strand != strand:
            continue

        if strand == "+":
            distance_to_pas = read.reference_end - pas
        else:
            distance_to_pas = pas - read.reference_start

        read_is_tail = False
        cigar_string = read.cigarstring
        if read.is_reverse:
            if cigar_pattern_2.match(cigar_string):
                unmatch = int(unmatch_pattern_2.search(read.cigarstring)[1])
                unmatch_seq = read.query_sequence[:unmatch][::-1]
                m = seqpattern_2.search(unmatch_seq)
                read_is_tail = m.group(0).count("T") >= pA_min_required_matches if m is not None else False
        else:
            if cigar_pattern_1.match(cigar_string):
                unmatch = int(unmatch_pattern_1.search(read.cigarstring)[1])
                unmatch_seq = read.query_sequence[-unmatch:]
                m = seqpattern_1.search(unmatch_seq)
                read_is_tail = m.group(0).count("A") >= pA_min_required_matches if m is not None else False

        reference_length = read.query_alignment_length
        data.append(
            {
                "distance_to_pas": distance_to_pas,
                "cigar_string": cigar_string,
                "reference_length": reference_length,
                "strand": read_strand,
                "is_tail": read_is_tail,
            }
        )

    if len(data) == 0:
        return pd.DataFrame()
    data = pd.DataFrame(data)
    data = data[
        (data["distance_to_pas"] >= -upstream) & (data["distance_to_pas"] <= downstream)
    ]
    data["is_gap"] = data["cigar_string"].str.contains(r"\d{2}N")
    data["pas"] = ":".join([chrom, str(row["start"]),str(row["end"]), strand])
    return data
