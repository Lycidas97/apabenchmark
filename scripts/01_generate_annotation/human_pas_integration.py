#!/usr/bin/env python3
"""Script export of human_pas_integration.ipynb."""


# %%
import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pybedtools import BedTool

# %%
# download the data from the following link

# %%
# set source data path
GENCODEPOLYA_PATH="../../data/raw_data/annotations/gencode.v40.polyAs.gtf"
POLYASITES_PATH="../../data/raw_data/annotations/atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz"
POLYADB_PATH="../../data/raw_data/annotations/HumanPas/hg38.PAS.main.tsv"
GENCODE_PATH="../../data/raw_data/annotations/gencode.v40.annotation.bed"

# set interval data path
GENCODEPOLYA_PROCESSED_PATH="../../data/raw_data/annotations/gencode.v40.polya.processed.bed"
POLYASITES_PROCESSED_PATH="../../data/raw_data/annotations/polyasites.human.processed.bed"
POLYADB_PROCESSED_PATH="../../data/raw_data/annotations/humanpas.hg38.processed.bed"

# set output data path
INTEGRATED_PAS_OUTPUT_PATH="../../data/int_data/annotations/human_integrated_pas.bed"
GENCODE_PROCESSED_PATH="../../data/raw_data/annotations/gencode.v40.annotation.bed"

# Venn overlap stats are expensive on full PolyASite3.0 inputs; enable only when needed.
RUN_VENN=False

os.makedirs("../../data/raw_data/annotations", exist_ok=True)
os.makedirs("../../data/int_data/annotations", exist_ok=True)

# %%
# process genome annotation
gencode_df = pd.read_csv(
    GENCODE_PATH, sep="\t", header=None
)
gencode_df = gencode_df.rename(
    columns={
        0: "chr",
        1: "start",
        2: "stop",
        4: "score",
        5: "strand",
        6: "source",
        7: "type",
        9: "annotation",
    }
)

# %%
gencode_exon_df = gencode_df[gencode_df.type == "exon"]
gencode_exon_df[['gene_id', 'transcript_id', 'gene_name', 'exon_number']] = gencode_exon_df['annotation'].str.extract('gene_id "([^"]*)";.*transcript_id "([^"]*)";.*gene_name "([^"]*)";.*exon_number ([0-9]*);')
gencode_exon_df['exon_number'] = gencode_exon_df['exon_number'].astype(int)

# classify exon
gencode_exon_df['exon_type'] = 'Internal exon'
first_exons = gencode_exon_df.groupby('transcript_id')['exon_number'].idxmin()
last_exons = gencode_exon_df.groupby('transcript_id')['exon_number'].idxmax()
gencode_exon_df.loc[first_exons, 'exon_type'] = '5\' most exon'
gencode_exon_df.loc[last_exons, 'exon_type'] = '3\' most exon'
exon_counts = gencode_exon_df.groupby('transcript_id').size()
single_exon_transcripts = exon_counts[exon_counts == 1].index
gencode_exon_df.loc[gencode_exon_df['transcript_id'].isin(single_exon_transcripts), 'exon_type'] = 'Single exon'

# %%
gencode_genebody_df = gencode_df[gencode_df.type == "gene"]
gencode_genebody_df[['gene_id', 'gene_name',]] = gencode_genebody_df['annotation'].str.extract('gene_id "([^"]*)";.*gene_name "([^"]*)";')
gencode_genebody_df = gencode_genebody_df.drop_duplicates(subset=['gene_id'])

# %%
# process gencode polya annotation
gencode_polya = pd.read_csv(
    GENCODEPOLYA_PATH, sep="\t", header=None, skiprows=5
)
gencode_polya = gencode_polya.rename(
    columns={
        0: "chr",
        1: "source",
        2: "tag",
        3: "start",
        4: "stop",
        5: "score",
        6: "strand",
        7: "phase",
        8: "annotation",
    }
)
gencode_polya["source"] = "Gencode"
gencode_polya["score"] = "1"
gencode_polya_to_save = gencode_polya[
    gencode_polya.tag == "polyA_site"
][["chr", "start", "stop", "source", "score", "strand"]]
gencode_polya_to_save["stop"] = gencode_polya_to_save.apply(lambda x: x["start"] + 1 if x["strand"] == "+" else x["stop"], axis=1)
gencode_polya_to_save["start"] = gencode_polya_to_save["stop"] - 1
gencode_polya_to_save.to_csv(
    GENCODEPOLYA_PROCESSED_PATH, index=False, header=None, sep="\t"
)

# %%
# process polyasite annotation
polyasite = pd.read_csv(
    POLYASITES_PATH, sep="	", header=None, compression="infer"
)
polyasite = polyasite.rename(
    columns={
        0: "chr",
        1: "start",
        2: "stop",
        3: "id",
        4: "score",
        5: "strand",
        9: "annotation",
    }
)
polyasite["chr"] = polyasite["chr"].astype(str)
missing_chr_prefix = ~polyasite["chr"].str.startswith("chr")
polyasite.loc[missing_chr_prefix, "chr"] = "chr" + polyasite.loc[missing_chr_prefix, "chr"]
polyasite["position"] = polyasite.id.astype(str).str.split(":", expand=True)[1].astype(int)
polyasite["start"] = np.where(
    polyasite.strand == "+", polyasite.position, polyasite.position - 1
)
polyasite["stop"] = np.where(
    polyasite.strand == "+", polyasite.position + 1, polyasite.position
)
polyasite["source"] = "PolyASite3.0"
polyasite["score"] = "1"
polyasite["annotation"] = polyasite["annotation"].fillna("unknown")
polyasite_to_save = polyasite[["chr", "start", "stop", "annotation", "score", "strand"]].copy()
polyasite_to_save["stop"] = np.where(
    polyasite_to_save["strand"] == "+", polyasite_to_save["start"] + 1, polyasite_to_save["stop"]
)
polyasite_to_save["start"] = polyasite_to_save["stop"] - 1
polyasite_to_save.to_csv(
    POLYASITES_PROCESSED_PATH, index=False, header=None, sep="	"
)

# %%
# process HumanPas annotation; coordinates are already hg38, so no liftover is needed
polyadb = pd.read_csv(
    POLYADB_PATH,
    sep="	",
    usecols=["PAS_ID", "GeneSymbol", "PAS_type_RefSeq_label"],
)
polyadb[["Chromosome", "Strand", "Position"]] = polyadb["PAS_ID"].str.extract(
    r"^(chr[^:]+):([+-]):(\d+)$"
)
polyadb = polyadb.dropna(subset=["Chromosome", "Strand", "Position"]).copy()
polyadb["Position"] = polyadb["Position"].astype(int)
polyadb["start"] = np.where(
    polyadb.Strand == "+", polyadb.Position, polyadb.Position - 1
)
polyadb["stop"] = np.where(
    polyadb.Strand == "+", polyadb.Position + 1, polyadb.Position
)
polyadb["GeneSymbol"] = polyadb["GeneSymbol"].fillna("unknown").astype(str)
polyadb["PAS_type_RefSeq_label"] = polyadb["PAS_type_RefSeq_label"].fillna("unknown").astype(str)
polyadb["annotation"] = (
    polyadb["PAS_type_RefSeq_label"].str.replace(" ", "_", regex=False)
    + ":"
    + polyadb["GeneSymbol"]
)
polyadb["score"] = "1"
polyadb_to_save = polyadb[
    ["Chromosome", "start", "stop", "annotation", "score", "Strand"]
].copy()
polyadb_to_save["stop"] = np.where(
    polyadb_to_save["Strand"] == "+", polyadb_to_save["start"] + 1, polyadb_to_save["stop"]
)
polyadb_to_save["start"] = polyadb_to_save["stop"] - 1
polyadb_to_save.to_csv(
    POLYADB_PROCESSED_PATH, index=False, header=None, sep="	"
)

# %%
gencode_exon_bed = BedTool.from_dataframe(gencode_exon_df)
gencode_genebody_bed = BedTool.from_dataframe(gencode_genebody_df)
polyasites_bed = BedTool.from_dataframe(polyasite_to_save)
polyadb_bed = BedTool(POLYADB_PROCESSED_PATH)
gencode_polya_bed = BedTool.from_dataframe(gencode_polya_to_save)

# %%
polyadb_bed

# %%
# define function to merge bed
def merge_bed(bed):
    bed_pos = bed.filter(lambda b: b.strand == '+').sort().merge(d=10, s=True, c='4,5,6', o='last')
    bed_neg = bed.filter(lambda b: b.strand == '-').sort().merge(d=10, s=True, c='4,5,6', o='first')
    bed_merged = bed_pos.cat(bed_neg, postmerge=False).sort().to_dataframe()
    bed_merged['start'] = bed_merged.apply(lambda x: x['end']-1 if x['strand'] == '+' else x['start'], axis=1)
    bed_merged['end'] = bed_merged['start'] + 1
    return BedTool.from_dataframe(bed_merged)


# 合并并去重各个数据集
gencode_polya_bed_merged = merge_bed(gencode_polya_bed)
polyadb_bed_merged = merge_bed(polyadb_bed)
merged_polyasites = merge_bed(polyasites_bed)

print(f"gencode_polya_bed_merged: {len(gencode_polya_bed_merged)}")
print(f"polyadb_bed_merged: {len(polyadb_bed_merged)}")
print(f"merged_polyasites: {len(merged_polyasites)}")

# 获取Gencode M25的PASs与HumanPas的PASs在±10nt范围内有重叠的PAS
overlapping_gencode_polyadb = gencode_polya_bed_merged.window(polyadb_bed_merged, w=20, sm=True, u=True)
print(f"overlapping_gencode_polyadb: {len(overlapping_gencode_polyadb)}")
polyadb_subtracted = polyadb_bed_merged.subtract(overlapping_gencode_polyadb, s=True, A=True)
print(f"polyadb_subtracted: {len(polyadb_subtracted)}")

# 合并gencode和polyadb_subtracted
merged_gencode_polyadb = gencode_polya_bed_merged.cat(polyadb_subtracted, postmerge=False).sort()

# 获取合并后的gencode_polyadb与PolyASite的PASs在±10nt范围内有重叠的PAS
overlapping_merged_polyasites = merged_gencode_polyadb.window(merged_polyasites, w=20, sm=True, u=True)
print(f"overlapping_merged_polyasites: {len(overlapping_merged_polyasites)}")
merged_polyasites_subtracted = merged_polyasites.subtract(overlapping_merged_polyasites, s=True, A=True)
print(f"merged_polyasites_subtracted: {len(merged_polyasites_subtracted)}")

# 最终合并所有数据集
merged_polyasites_addall = merged_gencode_polyadb.cat(merged_polyasites_subtracted, postmerge=False).sort()

print(f"Final merged PASs: {len(merged_polyasites_addall)}")

# %%
# # 定义函数用于合并数据集
# def merge_bed(bed):
#     bed_pos = bed.filter(lambda b: b.strand == '+').sort().merge(d=10, s=True, c='4,5,6', o='last')
#     bed_neg = bed.filter(lambda b: b.strand == '-').sort().merge(d=10, s=True, c='4,5,6', o='first')
#     bed_merged = bed_pos.cat(bed_neg, postmerge=False).sort().to_dataframe()
#     bed_merged['start'] = bed_merged.apply(lambda x: x['end']-1 if x['strand'] == '+' else x['start'], axis=1)
#     bed_merged['end'] = bed_merged['start'] + 1
#     return BedTool.from_dataframe(bed_merged)

# # 定义一个辅助函数用于获取重叠的PASs
# # def get_overlapping_pases(bed1, bed2):
# #     return bed1.window(bed1, w=20, sm=True, u=True)

# # 合并并去重各个数据集
# merged_polyasites = merge_bed(polyasites_bed)
# polyadb_bed_merged = merge_bed(polyadb_bed)
# gencode_polya_bed_merged = merge_bed(gencode_polya_bed)

# print(f"merged_polyasites: {len(merged_polyasites)}")

# # 获取HumanPas的PASs与当前PAS集合中在±10nt范围内有重叠的PAS
# # overlapping_polyadb = get_overlapping_pases(polyadb_bed_merged, merged_polyasites)
# overlapping_polyadb = polyadb_bed_merged.window(merged_polyasites, w=20, sm=True, u=True)
# print(f"overlapping_polyadb: {len(overlapping_polyadb)}")
# polyadb_subtracted = polyadb_bed_merged.subtract(overlapping_polyadb, s=True, A=True)
# print(f"polyadb_subtracted: {len(polyadb_subtracted)}")
# merged_polyasites_addpd = merged_polyasites.cat(polyadb_subtracted, postmerge=False).sort()

# # 同样的方法处理Gencode M25的PASs
# # overlapping_gencode = get_overlapping_pases(gencode_polya_bed_merged, merged_polyasites_addpd)
# overlapping_gencode = gencode_polya_bed_merged.window(merged_polyasites_addpd, w=20, sm=True, u=True)
# print(f"overlapping_gencode: {len(overlapping_gencode)}")
# gencode_polya_subtracted = gencode_polya_bed_merged.subtract(overlapping_gencode, s=True, A=True)
# print(f"gencode_polya_subtracted: {len(gencode_polya_subtracted)}")
# merged_polyasites_addall = merged_polyasites_addpd.cat(gencode_polya_subtracted, postmerge=False).sort()

# %%
if RUN_VENN:
    from matplotlib_venn import venn3
    import matplotlib.pyplot as plt

    from pybedtools import BedTool

    # polyasites_bed_merged = merge_bed(polyasites_bed)
    # polyadb_bed_merged = merge_bed(polyadb_bed)
    # gencode_polya_bed_merged = merge_bed(gencode_polya_bed)

    # 计算各个集合的大小
    polyasites_size = len(merged_polyasites_addall.window(merged_polyasites, w=20, sm=True, u=True))
    polyadb_size = len(merged_polyasites_addall.window(polyadb_bed_merged, w=20, sm=True, u=True))
    gencode_size = len(merged_polyasites_addall.window(gencode_polya_bed_merged, w=20, sm=True, u=True))

    # 计算交集的大小
    polyasites_polyadb_overlap = len(merged_polyasites_addall.window(merged_polyasites, w=20, sm=True, u=True).window(polyadb_bed_merged, w=20, sm=True, u=True))
    polyasites_gencode_overlap = len(merged_polyasites_addall.window(merged_polyasites, w=20, sm=True, u=True).window(gencode_polya_bed_merged, w=20, sm=True, u=True))
    polyadb_gencode_overlap = len(merged_polyasites_addall.window(gencode_polya_bed_merged, w=20, sm=True, u=True).window(polyadb_bed_merged, w=20, sm=True, u=True))

    # # 对window的结果进行处理，只保留前六列
    # def process_window_result(result):
    #     processed_result = []
    #     for line in str(result).split('\n'):
    #         columns = line.split('\t')
    #         if len(columns) >= 6:
    #             processed_result.append('\t'.join(columns[:6]))
    #     return BedTool('\n'.join(processed_result), from_string=True)

    # 计算三个数据集的交集的大小
    polyasites_polyadb_gencode_overlap = len(merged_polyasites_addall.window(merged_polyasites, w=20, sm=True, u=True).window(polyadb_bed_merged, w=20, sm=True, u=True).window(gencode_polya_bed_merged, w=20, sm=True, u=True))
    # 计算只在两个数据集中的元素的数量
    polyasites_polyadb_only = polyasites_polyadb_overlap - polyasites_polyadb_gencode_overlap
    polyasites_gencode_only = polyasites_gencode_overlap - polyasites_polyadb_gencode_overlap
    polyadb_gencode_only = polyadb_gencode_overlap - polyasites_polyadb_gencode_overlap

    # 计算只在一个数据集中的元素的数量
    polyasites_only = polyasites_size - polyasites_polyadb_overlap - polyasites_gencode_overlap + polyasites_polyadb_gencode_overlap
    polyadb_only = polyadb_size - polyasites_polyadb_overlap - polyadb_gencode_overlap + polyasites_polyadb_gencode_overlap
    gencode_only = gencode_size - polyasites_gencode_overlap - polyadb_gencode_overlap + polyasites_polyadb_gencode_overlap
else:
    print("Skipping Venn overlap calculation; set RUN_VENN=True to enable.")

# %%
if RUN_VENN:
    import pickle
    subsets=(polyasites_only, polyadb_only, polyasites_polyadb_only, gencode_only, polyasites_gencode_only, polyadb_gencode_only, polyasites_polyadb_gencode_overlap)
    with open("../../data/raw_data/annotations/hg38_venn.pkl", "wb") as f:
        pickle.dump(subsets, f)
    

# %%
if RUN_VENN:
    # 生成韦恩图
    with open("../../data/raw_data/annotations/hg38_venn.pkl", "rb") as f:
        subsets = pickle.load(f)
    venn3(subsets=subsets,
          set_labels=('PolyASite3.0', 'HumanPas', 'Gencode'))
    plt.show()
    plt.savefig("../../data/raw_data/annotations/hg38_venn3.svg", dpi=300)

# %%
# Filter 3' most exon and single exon
gencode_exon_df_filtered = gencode_exon_df[gencode_exon_df['exon_type'].isin(['3\' most exon', 'Single exon'])]
downstream_len=2000

# %%
import bisect
def find_next_start(current, sorted_list, strand):
    # 使用二分查找找到比current_stop大的最小的start
    if strand == '+':
        idx = bisect.bisect_right(sorted_list, current)
    elif strand == '-':
        idx = bisect.bisect_right(sorted_list, current) - 1
    else:
        raise ValueError('Strand must be + or -')

    if idx != len(sorted_list) and idx != -1:
        return sorted_list[idx]
    else:
        raise ValueError('No next start found')

chrs = gencode_genebody_df['chr'].unique()
gencode_downstream_df = pd.DataFrame()
for chr in chrs:
    chr_genbody_df = gencode_genebody_df[gencode_genebody_df['chr'] == chr]
    chr_exon_df = gencode_exon_df_filtered[gencode_exon_df_filtered['chr'] == chr]
    #chr_df = pd.merge(chr_genbody_df, chr_exon_df[['chr', 'start', 'stop', 'gene_id', 'next_start']], left_on='gene_id', right_on='gene_id', how='left')
    
    chr_genbody_df_pos = chr_genbody_df[chr_genbody_df['strand'] == '+'].sort_values(['start'])
    chr_genbody_df_neg = chr_genbody_df[chr_genbody_df['strand'] == '-'].sort_values(['stop'], ascending=False)
    chr_genbody_df_pos['next_start'] = chr_genbody_df_pos['stop'] + downstream_len
    chr_genbody_df_neg['next_start'] = chr_genbody_df_neg['start'] - downstream_len
    chr_genbody_df_pos.loc[chr_genbody_df_pos.index[:-1], 'next_start'] = chr_genbody_df_pos['start'].shift(-1)[:-1]
    chr_genbody_df_neg.loc[chr_genbody_df_neg.index[:-1], 'next_start'] = chr_genbody_df_neg['stop'].shift(-1)[:-1]
   
    exceptional_gene_pos = chr_genbody_df_pos[chr_genbody_df_pos['next_start'] < chr_genbody_df_pos['stop']]
    exceptional_gene_neg = chr_genbody_df_neg[chr_genbody_df_neg['next_start'] > chr_genbody_df_neg['start']]


    sorted_starts = sorted(chr_genbody_df_pos['start'].unique())
    sorted_stops = sorted(chr_genbody_df_neg['stop'].unique())

    if len(exceptional_gene_pos) > 0:
        for index, gene in exceptional_gene_pos.iterrows():
            try:
                chr_genbody_df_pos.loc[index, 'next_start'] = find_next_start(gene['stop'], sorted_starts, '+')
            except ValueError:
                chr_genbody_df_pos.loc[index, 'next_start'] = gene['stop'] + downstream_len
            except Exception as e:
                raise e
    if len(exceptional_gene_neg) > 0:
        for index, gene in exceptional_gene_neg.iterrows():
            try:
                chr_genbody_df_neg.loc[index, 'next_start'] = find_next_start(gene['start'], sorted_stops, '-')
            except ValueError:
                chr_genbody_df_neg.loc[index, 'next_start'] = gene['start'] - downstream_len
            except Exception as e:
                raise e
    chr_exon_df_pos = chr_exon_df[chr_exon_df['strand'] == '+'].sort_values(['start'])
    chr_exon_df_neg = chr_exon_df[chr_exon_df['strand'] == '-'].sort_values(['stop'], ascending=False)

    chr_exon_df_pos = pd.merge(chr_exon_df_pos, chr_genbody_df_pos[["gene_id","next_start"]],  left_on='gene_id', right_on='gene_id', how='left')
    chr_exon_df_neg = pd.merge(chr_exon_df_neg, chr_genbody_df_neg[["gene_id","next_start"]],  left_on='gene_id', right_on='gene_id', how='left')

    chr_exon_df_pos['downstream_start'] = chr_exon_df_pos["stop"]
    chr_exon_df_pos['downstream_stop'] = np.where(chr_exon_df_pos['downstream_start'] + downstream_len > chr_exon_df_pos['next_start'],
                                        chr_exon_df_pos['next_start'], chr_exon_df_pos['downstream_start'] + downstream_len)
    chr_exon_df_neg['downstream_stop'] = chr_exon_df_neg["start"]
    chr_exon_df_neg['downstream_start'] = np.where(chr_exon_df_neg['downstream_stop'] - downstream_len < chr_exon_df_neg['next_start'],
                                        chr_exon_df_neg['next_start'], chr_exon_df_neg['downstream_stop'] - downstream_len)
    
    chr_downstream_df = pd.concat([chr_exon_df_neg, chr_exon_df_pos])[['chr', 'downstream_start', 'downstream_stop', 'gene_id', 'score', 'strand', 'gene_name']]
    chr_downstream_df["exon_type"] = "downstream"

    gencode_downstream_df = pd.concat([gencode_downstream_df, chr_downstream_df])
gencode_downstream_df = gencode_downstream_df.rename(columns={"downstream_start": "start", "downstream_stop": "stop"})
gencode_downstream_df = gencode_downstream_df.drop_duplicates(['chr', 'start', 'stop', 'gene_id', 'score', 'strand', 'gene_name'])

# %%

gencode_exon_df_final = pd.concat([gencode_exon_df, gencode_downstream_df])
gencode_exon_df_final = gencode_exon_df_final.rename(columns={"stop": "end", "chr":"chrom"})
gencode_exon_df_final = gencode_exon_df_final.drop_duplicates(["chrom", "start", "end", "gene_id", "exon_type", "strand"], keep="first")

gencode_exon_df_final["name"] = gencode_exon_df_final["gene_id"].str.cat(gencode_exon_df_final[["gene_name"]], sep=":")
exon_type_dict = {
    "5' most exon": "FE",
    "Internal exon": "IE",
    "3' most exon": "TE",
    "Single exon": "SE",
    "downstream": "DS",
    "Intron": "IT",
}

gencode_exon_df_final["exon_type"] = gencode_exon_df_final["exon_type"].map(exon_type_dict)

gencode_te_to_collapse = gencode_exon_df_final[gencode_exon_df_final['exon_type'].isin(["TE"])]
gencode_te_to_collapse_bed = BedTool.from_dataframe(gencode_te_to_collapse[["chrom", "start", "end", "name", "score", "strand"]])
gencode_te = gencode_te_to_collapse_bed.merge(s=True, c="4,5,6", o="distinct").to_dataframe()
gencode_te["name"] = gencode_te["name"].str.split(",")
gencode_te = gencode_te.explode("name")


gencode_te["exon_type"] = "TE"
gencode_te["exon_id"] = gencode_te.apply(lambda x: f"{x['name']}:{x.chrom}:{x.start}:{x.end}:{x.strand}:TE", axis=1)

gencode_se_to_collapse = gencode_exon_df_final[gencode_exon_df_final['exon_type'].isin(["SE"])]
gencode_se_to_collapse_bed = BedTool.from_dataframe(gencode_se_to_collapse[["chrom", "start", "end", "name", "score", "strand"]])
gencode_se = gencode_se_to_collapse_bed.merge(s=True, c="4,5,6", o="distinct").to_dataframe()
gencode_se["name"] = gencode_se["name"].str.split(",")
gencode_se = gencode_se.explode("name")


gencode_se["exon_type"] = "SE"
gencode_se["exon_id"] = gencode_se.apply(lambda x: f"{x['name']}:{x.chrom}:{x.start}:{x.end}:{x.strand}:SE", axis=1)

gencode_exon_df_final = gencode_exon_df_final[~gencode_exon_df_final['exon_type'].isin(["TE", "SE"])]
gencode_exon_df_final["exon_id"] = gencode_exon_df_final.apply(lambda x: f"{x['name']}:{x.chrom}:{x.start}:{x.end}:{x.strand}:{x.exon_type}", axis=1)
gencode_exon_df_final = pd.concat([gencode_exon_df_final, gencode_te, gencode_se])
gencode_exon_df_final = gencode_exon_df_final.loc[:,["chrom", "start", "end", "exon_id", "score", "strand"]]

# %%
gencode_genebody_df_final = pd.concat([gencode_genebody_df, gencode_downstream_df])
gencode_genebody_df_final = gencode_genebody_df_final.rename(columns={"stop": "end", "chr": "chrom"})
gencode_genebody_df_final["name"] = gencode_genebody_df_final["gene_id"].str.cat(gencode_genebody_df_final[["gene_name"]], sep=":")
gencode_genebody_df_final = gencode_genebody_df_final.drop_duplicates(["chrom", "start", "end", "gene_id", "strand", "gene_name"], keep="first")
gencode_genebody_df_final = gencode_genebody_df_final.loc[:,["chrom", "start", "end", "name", "score", "strand"]]

gencode_exon_bed = BedTool.from_dataframe(gencode_exon_df_final[["chrom", "start", "end", "exon_id", "score", "strand"]])
gencode_genebody_bed = BedTool.from_dataframe(gencode_genebody_df_final[["chrom", "start", "end", "name", "score", "strand"]])

gencode_genebody_antisense_df = gencode_genebody_df_final.copy()
gencode_genebody_antisense_df["strand"] = gencode_genebody_antisense_df["strand"].map({"+": "-", "-": "+"})
gencode_genebody_antisense_bed = BedTool.from_dataframe(gencode_genebody_antisense_df[["chrom", "start", "end", "name", "score", "strand"]])
gencode_genebody_antisense_bed = gencode_genebody_antisense_bed.subtract(gencode_genebody_bed, s=True, A=True,).sort()
gencode_genebody_antisense_df = gencode_genebody_antisense_bed.to_dataframe()

# %%
gencode_exon_antisense_df = gencode_exon_df_final.copy()
gencode_exon_antisense_df["strand"] = gencode_exon_antisense_df["strand"].map({"+": "-", "-": "+"})
gencode_exon_antisense_bed = BedTool.from_dataframe(gencode_exon_antisense_df[["chrom", "start", "end", "exon_id", "score", "strand"]])
gencode_exon_antisense_bed = gencode_exon_antisense_bed.subtract(gencode_genebody_bed, s=True, A=True,).sort()
gencode_exon_antisense_df = gencode_exon_antisense_bed.to_dataframe()
gencode_exon_antisense_df = gencode_exon_antisense_df.rename(columns={"name": "exon_id"})

# %%
integrated_pas = merged_polyasites_addall.to_dataframe()
integrated_pas.columns = ["chr", "start", "stop", "name", "score", "strand"]
integrated_pas["name"] = integrated_pas["chr"].str.cat(integrated_pas[["start", "stop", "strand",]].astype(str), sep=":")

# %%
# gencode_exon_bed = BedTool.from_dataframe(gencode_exon_df_merged[["chrom", "start", "end", "exon_id", "score", "strand"]])
gencode_exon_bed = BedTool.from_dataframe(gencode_exon_df_final[["chrom", "start", "end", "exon_id", "score", "strand"]])
gencode_exon_antisense_bed = BedTool.from_dataframe(gencode_exon_antisense_df[["chrom", "start", "end", "exon_id", "score", "strand"]])
gencode_genebody_bed = BedTool.from_dataframe(gencode_genebody_df_final[["chrom", "start", "end", "name", "score", "strand"]])
gencode_genebody_antisense_bed = BedTool.from_dataframe(gencode_genebody_antisense_df[["chrom", "start", "end", "name", "score", "strand"]])

integrated_pas_bed = BedTool.from_dataframe(integrated_pas)
intersected_pas_genebody = integrated_pas_bed.intersect(gencode_genebody_bed, wa=True, wb=True, s=True)
intersected_pas_exon = integrated_pas_bed.intersect(gencode_exon_bed, wa=True, wb=True, s=True)
intersected_pas_genebody_antisense = integrated_pas_bed.intersect(gencode_genebody_antisense_bed, wa=True, wb=True, s=True)
intersected_pas_exon_antisense = integrated_pas_bed.intersect(gencode_exon_antisense_bed, wa=True, wb=True, s=True)

# %%
def _empty_annotation_frame():
    return pd.DataFrame(columns=["name", "gene_id", "gene_name", "pas_type", "exon_id"])


def _read_intersection_columns(bedtool, usecols, names):
    try:
        df = pd.read_csv(
            bedtool.fn,
            sep="\t",
            header=None,
            usecols=usecols,
            dtype="string",
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=names)
    df.columns = names
    return df


def _parse_genebody_intersection(bedtool):
    df = _read_intersection_columns(bedtool, [3, 9], ["name", "gene_info"])
    if df.empty:
        return _empty_annotation_frame()
    gene_parts = df["gene_info"].str.split(":", n=1, expand=True)
    df["gene_id"] = gene_parts[0]
    df["gene_name"] = gene_parts[1] if gene_parts.shape[1] > 1 else "unknown"
    df["pas_type"] = pd.NA
    df["exon_id"] = "unknown"
    return df[["name", "gene_id", "gene_name", "pas_type", "exon_id"]]


def _parse_exon_intersection(bedtool, antisense=False):
    df = _read_intersection_columns(bedtool, [3, 9], ["name", "exon_id"])
    if df.empty:
        return _empty_annotation_frame()
    exon_parts = df["exon_id"].str.split(":", expand=True)
    df["gene_id"] = exon_parts[0]
    df["gene_name"] = exon_parts[1] if exon_parts.shape[1] > 1 else "unknown"
    df["pas_type"] = exon_parts[6] if exon_parts.shape[1] > 6 else "unknown"
    if antisense:
        df["pas_type"] = df["pas_type"].replace(
            {"TE": "AE", "IE": "AE", "FE": "AE", "SE": "AE", "DS": "AU"}
        )
    return df[["name", "gene_id", "gene_name", "pas_type", "exon_id"]]


intersected_pas_genebody_df = _parse_genebody_intersection(intersected_pas_genebody)
intersected_pas_exon_df = _parse_exon_intersection(intersected_pas_exon)
intersected_pas_genebody_antisense_df = _parse_genebody_intersection(intersected_pas_genebody_antisense)
intersected_pas_exon_antisense_df = _parse_exon_intersection(
    intersected_pas_exon_antisense, antisense=True
)

exon_names = pd.Index(intersected_pas_exon_df["name"].drop_duplicates())
antisense_exon_names = pd.Index(
    intersected_pas_exon_antisense_df.loc[
        intersected_pas_exon_antisense_df["pas_type"] != "AU", "name"
    ].drop_duplicates()
)
antisense_upstream_names = pd.Index(
    intersected_pas_exon_antisense_df.loc[
        intersected_pas_exon_antisense_df["pas_type"] == "AU", "name"
    ].drop_duplicates()
)
antisense_non_intron_names = antisense_exon_names.union(antisense_upstream_names)

intersected_pas_intron_df = intersected_pas_genebody_df.loc[
    ~intersected_pas_genebody_df["name"].isin(exon_names)
].copy()
intersected_pas_intron_df["pas_type"] = "IT"

intersected_pas_antisense_intron_df = intersected_pas_genebody_antisense_df.loc[
    ~intersected_pas_genebody_antisense_df["name"].isin(antisense_non_intron_names)
].copy()
intersected_pas_antisense_intron_df["pas_type"] = "AI"

intersected_pas_antisense_exon_df = intersected_pas_exon_antisense_df.loc[
    intersected_pas_exon_antisense_df["name"].isin(antisense_exon_names)
].copy()
intersected_pas_antisense_exon_df["pas_type"] = "AE"

intersected_pas_antisense_upstream_df = intersected_pas_exon_antisense_df.loc[
    intersected_pas_exon_antisense_df["name"].isin(antisense_upstream_names)
].copy()
intersected_pas_antisense_upstream_df["pas_type"] = "AU"

annotation_cols = ["name", "gene_id", "gene_name", "pas_type", "exon_id"]
intersected_pas_df = pd.concat(
    [
        intersected_pas_exon_df,
        intersected_pas_intron_df,
        intersected_pas_antisense_exon_df,
        intersected_pas_antisense_intron_df,
        intersected_pas_antisense_upstream_df,
    ],
    ignore_index=True,
)[annotation_cols]

priority_dict = {
    "TE": 5,
    "SE": 4,
    "IE": 3,
    "FE": 2,
    "IT": 1,
    "DS": 0,
    "AE": -1,
    "AI": -2,
    "AU": -3,
}
intersected_pas_df["priority"] = intersected_pas_df["pas_type"].map(priority_dict).astype("Int8")

te_se_mask = intersected_pas_df["pas_type"].isin(["TE", "SE"])
te_se_df = intersected_pas_df.loc[te_se_mask, annotation_cols].copy()
te_se_names = pd.Index(te_se_df["name"].drop_duplicates())
other_df = intersected_pas_df.loc[
    ~te_se_mask & ~intersected_pas_df["name"].isin(te_se_names),
    annotation_cols + ["priority"],
].copy()

if other_df.empty:
    other_df_drop_duplicates = pd.DataFrame(columns=annotation_cols)
else:
    max_priority = other_df.groupby("name")["priority"].transform("max")
    other_df_multiple_max = other_df.loc[other_df["priority"] == max_priority]
    ambiguous_names = other_df_multiple_max.drop_duplicates(["name", "gene_id"]).groupby("name").size()
    ambiguous_names = pd.Index(ambiguous_names[ambiguous_names > 1].index)

    best_idx = other_df.groupby("name")["priority"].idxmax()
    other_df_drop_duplicates = other_df.loc[best_idx, annotation_cols].copy()
    other_df_drop_duplicates.loc[
        other_df_drop_duplicates["name"].isin(ambiguous_names), ["gene_id", "gene_name"]
    ] = "unknown"

result_df = pd.concat([te_se_df, other_df_drop_duplicates], ignore_index=True)
print(
    "Annotation rows: "
    f"exon={len(intersected_pas_exon_df)}, "
    f"intron={len(intersected_pas_intron_df)}, "
    f"antisense={len(intersected_pas_antisense_exon_df) + len(intersected_pas_antisense_intron_df) + len(intersected_pas_antisense_upstream_df)}, "
    f"selected={len(result_df)}"
)

# %%
integrated_pas_final = integrated_pas.merge(result_df, on="name", how="left", sort=False)
integrated_pas_final = integrated_pas_final.rename(columns={"chr": "chrom"})
integrated_pas_final["pas_type"] = integrated_pas_final["pas_type"].fillna("IG")
integrated_pas_final.loc[~integrated_pas_final["pas_type"].isin(["TE", "SE"]), "exon_id"] = "unknown"
integrated_pas_final = integrated_pas_final.fillna("unknown")

os.makedirs(os.path.dirname(INTEGRATED_PAS_OUTPUT_PATH), exist_ok=True)
integrated_pas_final[
    ["chrom", "start", "stop", "name", "score", "strand", "gene_id", "gene_name", "pas_type", "exon_id"]
].to_csv(INTEGRATED_PAS_OUTPUT_PATH, sep="\t", index=False, header=False)
print(f"Wrote {len(integrated_pas_final)} rows to {INTEGRATED_PAS_OUTPUT_PATH}")

# %%
# Replaced by the optimized annotation merge above.

# %%
# Replaced by the optimized annotation merge above.

# %%
# Replaced by the optimized annotation merge above.

# %%
# Replaced by the optimized annotation merge above.

# %%
# Replaced by the optimized annotation merge above.

# %%
# Replaced by the optimized annotation merge above.

# %%
# Replaced by the optimized annotation merge above.

# %%
# gencode_genebody_df_final = pd.concat([gencode_genebody_df, gencode_downstream_df])
# gencode_exon_df_final = pd.concat([gencode_exon_df, gencode_downstream_df])
# gencode_exon_df_final = gencode_exon_df_final.drop_duplicates(["chr", "start", "stop", "gene_id", "exon_type", "strand"], keep="first")

# gencode_exon_df_final["name"] = gencode_exon_df_final["gene_id"].str.cat(gencode_exon_df_final[["gene_name", "exon_type"]], sep=":")
# gencode_genebody_df_final["name"] = gencode_genebody_df_final["gene_id"].str.cat(gencode_genebody_df_final[["gene_name"]], sep=":")

# gencode_exon_bed = BedTool.from_dataframe(gencode_exon_df_final[["chr", "start", "stop", "name", "score", "strand"]])
# gencode_genebody_bed = BedTool.from_dataframe(gencode_genebody_df_final[["chr", "start", "stop", "name", "score", "strand"]])

# %%
# integrated_pas = merged_polyasites_addall.to_dataframe()
# integrated_pas.columns = ["chr", "start", "stop", "name", "score", "strand"]
# integrated_pas["name"] = integrated_pas["chr"].str.cat(integrated_pas[["start", "stop", "strand",]].astype(str), sep=":")

# %%
# integrated_pas

# %%
# integrated_pas_bed = BedTool.from_dataframe(integrated_pas)
# intersected_pas_genebody = integrated_pas_bed.intersect(gencode_genebody_bed, wa=True, wb=True, s=True)
# intersected_pas_exon = integrated_pas_bed.intersect(gencode_exon_bed, wa=True, wb=True, s=True)

# %%
# intersected_pas_genebody_df = intersected_pas_genebody.to_dataframe()
# intersected_pas_exon_df = intersected_pas_exon.to_dataframe()

# intersected_pas_genebody_df[['gene_id', 'gene_name']] = intersected_pas_genebody_df.iloc[:,9].str.split(':', expand=True)
# intersected_pas_exon_df[['gene_id', 'gene_name', 'exon_type']] = intersected_pas_exon_df.iloc[:,9].str.split(':', expand=True)

# %%
# intronic_pas = intersected_pas_genebody_df[~intersected_pas_genebody_df["name"].isin(intersected_pas_exon_df["name"])]["name"].copy().drop_duplicates().tolist()
# intersected_pas_intron_df = intersected_pas_genebody_df[intersected_pas_genebody_df["name"].isin(intronic_pas)].copy()
# intersected_pas_intron_df["exon_type"] = "intron"
# intersected_pas_df = pd.concat([intersected_pas_exon_df, intersected_pas_intron_df])

# %%
# priority_dict = {
#     "3' most exon": 5,
#     "Single exon": 4,
#     "Internal exon": 3,
#     "5' most exon": 2,
#     "downstream": 0,
#     "intron": 1
# }
# intersected_pas_df['priority'] = intersected_pas_df['exon_type'].map(priority_dict)

# %%
# intersected_pas_df_max_priority = intersected_pas_df.groupby('name')['priority'].max()
# intersected_pas_df_with_max_priority = intersected_pas_df.merge(intersected_pas_df_max_priority, on='name', suffixes=('', '_max'))
# intersected_pas_df_multiple_max = intersected_pas_df_with_max_priority[intersected_pas_df_with_max_priority['priority'] == intersected_pas_df_with_max_priority['priority_max']]
# grouped_df = intersected_pas_df_multiple_max.drop_duplicates(["name", "gene_id"]).groupby("name").size()
# unknown_pas = grouped_df[grouped_df > 1].index.tolist()
# intersected_pas_df_drop_duplicates = intersected_pas_df.sample(frac=1).sort_values(['name', 'priority'], ascending=[True, False]).groupby('name').first().reset_index()
# intersected_pas_df_drop_duplicates.loc[intersected_pas_df_drop_duplicates["name"].isin(unknown_pas), ["gene_id", "gene_name"]] = "unknown"

# %%
# intronic_pas = intersected_pas_df_drop_duplicates[intersected_pas_df_drop_duplicates["exon_type"] == "intron"]["name"].copy().drop_duplicates().tolist()
# most3Exon = intersected_pas_df_drop_duplicates[intersected_pas_df_drop_duplicates["exon_type"] == "3' most exon"]["name"].copy().drop_duplicates().tolist()
# most5Exon = intersected_pas_df_drop_duplicates[intersected_pas_df_drop_duplicates["exon_type"] == "5' most exon"]["name"].copy().drop_duplicates().tolist()
# singleExon = intersected_pas_df_drop_duplicates[intersected_pas_df_drop_duplicates["exon_type"] == "Single exon"]["name"].copy().drop_duplicates().tolist()
# internalExon = intersected_pas_df_drop_duplicates[intersected_pas_df_drop_duplicates["exon_type"] == "Internal exon"]["name"].copy().drop_duplicates().tolist()
# downstream = intersected_pas_df_drop_duplicates[intersected_pas_df_drop_duplicates["exon_type"] == "downstream"]["name"].copy().drop_duplicates().tolist()

# %%
# integrated_pas["pas_type"] = "intergenic"
# integrated_pas.loc[integrated_pas["name"].isin(intronic_pas), "pas_type"] = "Intron"
# integrated_pas.loc[integrated_pas["name"].isin(most5Exon), "pas_type"] = "5' most exon"
# integrated_pas.loc[integrated_pas["name"].isin(most3Exon), "pas_type"] = "3' most exon"
# integrated_pas.loc[integrated_pas["name"].isin(downstream), "pas_type"] = "downstream"
# integrated_pas.loc[integrated_pas["name"].isin(singleExon), "pas_type"] = "Single exon"
# integrated_pas.loc[integrated_pas["name"].isin(internalExon), "pas_type"] = "Internal exon"

# %%
# integrated_pas_final = pd.merge(left=integrated_pas, right=intersected_pas_df_drop_duplicates, on="name", how="left", suffixes=('', '_y'), ).fillna("unknown").loc[:, ["chr", "start", "stop", "name", "score", "strand", "pas_type", "gene_id", "gene_name"]]
# integrated_pas_final_bed_df = integrated_pas_final.copy()
# integrated_pas_final_bed = BedTool.from_dataframe(integrated_pas_final_bed_df)
# # add polyadb exon information for pas
# # intersect_polyadb = get_overlapping_pases(integrated_pas_final_bed, polyadb_bed_merged).to_dataframe(header=None)
# intersect_polyadb = integrated_pas_final_bed.window(polyadb_bed_merged, w=20, sm=True).to_dataframe(header=None)
# intersect_polyadb[["polyadb_type","polyadb_gene"]] = intersect_polyadb.iloc[:,12].str.split(':', expand=True)
# intersect_polyadb = intersect_polyadb[~intersect_polyadb["polyadb_gene"].isin(["unknown", "na", "nan"])]
# intersect_polyadb["polyadb_type"] = intersect_polyadb["polyadb_type"].str.replace("_", " ")
# intersect_polyadb = intersect_polyadb.rename(columns={3: "name", 6: "pas_type"})
# intersect_polyadb = intersect_polyadb[["name", "polyadb_type"]]
# # add polyasites TE information for pas
# # intersect_polyasites = get_overlapping_pases(integrated_pas_final_bed, polyasites_bed_merged).to_dataframe(header=None)
# intersect_polyasites = integrated_pas_final_bed.window(merged_polyasites, w=20, sm=True).to_dataframe(header=None)
# intersect_polyasites = intersect_polyasites.rename(columns={3: "name", 12:"polyasites_type"})
# intersect_polyasites = intersect_polyasites[intersect_polyasites["polyasites_type"] == "TE"][["name", "polyasites_type"]]
# intersect_polyasites["polyasites_type"] = "3' most exon"

# integrated_pas_final_add_annotation = pd.merge(integrated_pas_final, intersect_polyadb, how="left")
# integrated_pas_final_add_annotation = pd.merge(integrated_pas_final_add_annotation, intersect_polyasites, how="left")

# integrated_pas_final_add_annotation["polyadb_type"].combine_first(integrated_pas_final_add_annotation["polyasites_type"]).combine_first(integrated_pas_final_add_annotation["pas_type"])
# integrated_pas_final_add_annotation["integrated_pas_type"] = integrated_pas_final_add_annotation["polyadb_type"].combine_first(integrated_pas_final_add_annotation["polyasites_type"]).combine_first(integrated_pas_final_add_annotation["pas_type"])

# %%
# integrated_pas_final = integrated_pas_final_add_annotation.loc[:, ["chr", "start", "stop", "name", "score", "strand", "integrated_pas_type", "gene_id", "gene_name"]]

# %%
# integrated_pas_final.to_csv("../../data/int_data/annotations/human_integrated_pas.bed", sep="\t", index=False, header=False)
