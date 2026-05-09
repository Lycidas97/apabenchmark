import pandas as pd
from pybedtools import BedTool, create_interval_from_list
from Bio import SeqIO
import re
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import pysam




def match_seq(pattern, mrm, letter, unmatch_seq):
    m = pattern.search(unmatch_seq)
    flag = m.group(0).count(letter) >= mrm if m is not None else False
    return flag

def find_cs_by_chrom(
        bam_path: str,
        chrom: str,
        allowed_mismatches_at_start: int=2,
        match_length: int=10,
        min_required_matches: int=8,
    ) -> pd.DataFrame:
    """
    Identify and count cs in sequences for a specific chromosome in a BAM file.

    This function identifies cs in sequences for a specific chromosome based on 
    specific CIGAR string patterns and sequence patterns. It returns a DataFrame 
    where each row represents a unique combination of chromosome, strand, and 
    coordinate, with an additional column for the count of such combinations in 
    the BAM file.

    Parameters
    ----------
    bam_path : str
        Path to the BAM file.
    chrom : str
        The specific chromosome to process.

    Returns
    -------
    pd.DataFrame
        DataFrame where each row represents a unique combination of chromosome, 
        strand, and coordinate, with an additional column for the count of such 
        combinations.
    """
    bam = pysam.AlignmentFile(bam_path, 'rb')
    cigar_pattern_1 = re.compile(r'[0-9]*M[0-9]{2,}S')
    cigar_pattern_2 = re.compile(r'[0-9]{2,}S[0-9]*M')

    unmatch_pattern_1 = re.compile(r"(\d+)S$")
    unmatch_pattern_2 = re.compile(r"^(\d+)S")

    seqpattern_1 = re.compile(fr'^[^A]{{0,{allowed_mismatches_at_start}}}A[A-Z]{{{match_length-1}}}')
    seqpattern_2 = re.compile(fr'^[^T]{{0,{allowed_mismatches_at_start}}}T[A-Z]{{{match_length-1}}}')
    # seqpattern_1 = re.compile(r'(A{3,}[^A]{0,2})*A{6,}([^A]{0,2}A{3,})*.{0,2}?$')
    # seqpattern_2 = re.compile(r'(T{3,}[^T]{0,2})*T{6,}([^T]{0,2}T{3,})*.{0,2}?$')

    sequences = [read for read in bam.fetch(contig=chrom)]
    strand = '+'
    s_1 = ([read.reference_name, '+', read.reference_end, read.cigarstring, read.query_sequence[-int(unmatch_pattern_1.search(read.cigarstring)[1]):]] for read in sequences
            if cigar_pattern_1.match(read.cigarstring) and read.is_reverse == False)
    cs_pos = [s for s in s_1 if match_seq(seqpattern_1, min_required_matches, "A", s[4])]

    strand = '-'
    s_2 = ([read.reference_name, '-', read.reference_start, read.cigarstring, read.query_sequence[:int(unmatch_pattern_2.search(read.cigarstring)[1])][::-1]] for read in sequences
            if cigar_pattern_2.match(read.cigarstring) and read.is_reverse == True)
    cs_neg = [s for s in s_2 if match_seq(seqpattern_2, min_required_matches, "T", s[4])]


    final_cs_pos = pd.DataFrame(cs_pos, columns=['chrom', 'strand', 'coord', 'cigar', 'unmatch_seq']).drop(columns=['cigar', 'unmatch_seq'])
    final_cs_neg = pd.DataFrame(cs_neg, columns=['chrom', 'strand', 'coord', 'cigar', 'unmatch_seq']).drop(columns=['cigar', 'unmatch_seq'])


    cs = pd.concat([final_cs_pos, final_cs_neg])
    cs_grouped = cs.groupby(['chrom', 'strand', 'coord']).size().reset_index(name='count')

    return cs_grouped



def find_cs(
        bam_path,
        chrom_list,
        chrom_starts_with_chr,
        genecoord_bed,
        max_gap=10, 
        num_processors=8,
        allowed_mismatches_at_start: int=2,
        match_length: int=10,
        min_required_matches: int=8,
        ):

    with ProcessPoolExecutor(max_workers=num_processors) as executor:
        results = list(
            tqdm(executor.map(
                find_cs_by_chrom,
                [bam_path]*len(chrom_list),
                chrom_list,
                [allowed_mismatches_at_start]*len(chrom_list),
                [match_length]*len(chrom_list),
                [min_required_matches]*len(chrom_list)),
                total=len(chrom_list)
                ))
    
    cs_df = pd.concat(results)
    cs_df['start'] = cs_df.apply(lambda x: x['coord'] - 1 if x['strand'] == '+' else x['coord'], axis=1)
    cs_df['end'] = cs_df.apply(lambda x: x['coord'] if x['strand'] == '+' else x['coord'] + 1, axis=1)
    cs_df['name'] = '.'
    cs_df['score'] = cs_df['count']
    
    if chrom_starts_with_chr:
        pass
    else:
        cs_df['chrom'] = 'chr' + cs_df['chrom']

    cs_bed = BedTool.from_dataframe(cs_df[['chrom', 'start', 'end', 'name', 'score', 'strand', 'count']])
    def merge_bed(bed, mg):
            bed_pos = bed.filter(lambda b: b.strand == '+').sort().merge(d=mg, s=True, c='4,5,6', o='last,sum,last')
            bed_neg = bed.filter(lambda b: b.strand == '-').sort().merge(d=mg, s=True, c='4,5,6', o='first,sum,first')
            return bed_pos.cat(bed_neg, postmerge=False).sort()
    cs_bed = merge_bed(cs_bed, max_gap)
    cs_bed_df = cs_bed.to_dataframe()
    cs_bed_df["start"] = cs_bed_df.apply(lambda x: x.start if x.strand == '-' else x.end-1, axis=1)
    cs_bed_df["end"] = cs_bed_df.apply(lambda x: x.end if x.strand == '+' else x.start+1, axis=1)

    onsite_bed = cs_bed.intersect(genecoord_bed, wa=True, s=True, u=True)
    onsite_bed_df = onsite_bed.to_dataframe()
    onsite_bed_df["start"] = onsite_bed_df.apply(lambda x: x.start if x.strand == '-' else x.end-1, axis=1)
    onsite_bed_df["end"] = onsite_bed_df.apply(lambda x: x.end if x.strand == '+' else x.start+1, axis=1)
    return cs_bed_df, onsite_bed_df


def get_gene_coordinates_bed(gtf_file: str, extension_length: int) -> BedTool:
    """
    Get gene coordinates from a GTF file, extend them by a specified length, and
    convert them into a BedTool object.

    This function reads a GTF file into a DataFrame, selects rows where the feature 
    is 'gene', extends the coordinates of these genes by a specified length at downstream,
    converts the coordinates into a BedTool object, and returns this object.

    Parameters
    ----------
    gtf_file : str
        Path to the GTF file.
    extension_length : int
        Length by which to extend the gene coordinates.

    Returns
    -------
    BedTool
        BedTool object with the extended gene coordinates.
    """

    df = pd.read_csv(gtf_file, comment='#', sep='\t', header=None, 
                     names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

    df_genes = df[df['feature'] == 'gene'].copy()

    df_genes.loc[df_genes['strand'] == '+', 'end'] = df_genes.loc[df_genes['strand'] == '+', 'end'] + extension_length
    df_genes.loc[df_genes['strand'] == '-', 'start'] = df_genes.loc[df_genes['strand'] == '-', 'start'] - extension_length

    df_genes['start'] = df_genes['start'].clip(lower=1)

    bed_intervals = []
    for row in df_genes.itertuples():
        interval = create_interval_from_list([row.seqname, str(row.start), str(row.end), '.', '0', row.strand])
        bed_intervals.append(interval)
    
    result_bed = BedTool(bed_intervals)
    return result_bed