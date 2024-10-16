import numpy as np
import pandas as pd
from pybedtools import BedTool, create_interval_from_list
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
# SAMPLE_GENE_NUM = 4000
# INTERGENIC_NUM = 1000
# MIN_GAP_LEN = 50
# PA_RATIO = 0.05
# MAX_PAS_PER_GENE = 10
# PAS_PATH = "./data/pas/mouse_integrated_pas.bed"

GENE_NUMBERS = [1000, 2000, 4000, 8000, 16000]
# DIFF_APA_TE_NUM = 2000
# NON_DIFF_APA_TE_NUM = 2000

MAX_PAS_PER_TE = 5
MIN_GAP_LEN = 100
PAS_PATH = "../../data/int_data/annotations/mouse_pas.bed"


# filter pas by chr
pas_df = pd.read_csv(PAS_PATH, sep='\t', header=None)
pas_df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'gene_id', 'gene_name', "pas_type", "exon_id"]

# filter pas by chr
pas_df = pas_df[pas_df["chr"].str.contains(r"^chr[0-9XY]")]
duplicated_pas = pas_df.groupby(["name"]).filter(lambda x: len(x) > 1)["name"].unique()
pas_df = pas_df.loc[~pas_df["name"].isin(duplicated_pas)]

# grep single pas gene
single_pas_gene_list = pas_df.groupby("gene_id").filter(lambda x: set(x['pas_type']).issubset({'SE', 'TE'})).groupby("gene_id").filter(lambda x: len(x) == 1)["gene_id"].drop_duplicates()

pas_df_gene = pas_df[(pas_df['pas_type'].isin(["SE", "TE"]))]

pas_df_gene_forward = pas_df_gene[pas_df_gene["strand"] == "+"].sort_values(["chr", "start"])
pas_df_gene_reverse = pas_df_gene[pas_df_gene["strand"] == "-"].sort_values(["chr", "start"])

pas_df_gene_forward["gap"] = pas_df_gene_forward.groupby("gene_id")["start"].shift(-1) - pas_df_gene_forward["end"]
pas_df_gene_reverse["gap"] = pas_df_gene_reverse.groupby("gene_id")["start"].shift(-1) - pas_df_gene_reverse["end"]

pas_df_gene = pd.concat([pas_df_gene_forward, pas_df_gene_reverse])
pas_df_gene = pas_df_gene[(pas_df_gene["gap"] >= MIN_GAP_LEN) | (pas_df_gene["gap"] < 0) | (pas_df_gene["gap"].isna())]

multi_pas_te_list = pas_df_gene.groupby("exon_id").filter(lambda x: len(x) > 1)["exon_id"].drop_duplicates()


for gene_number in GENE_NUMBERS:
    for rep in range(1, 3):
        sample_bed = f"../../data/sim_data/sim_pas_for_cprsb/mm10_sim_pas_gn{gene_number}_rep{rep}.bed"
        
        # calculate the number of each gene type based on the ratio
        single_pas_gene_num = int(gene_number * 0.1)
        diff_apa_te_num = int(gene_number * 0.45)
        non_diff_apa_te_num = gene_number - single_pas_gene_num - diff_apa_te_num
        
        # sample genes
        single_pas_genes = single_pas_gene_list.sample(single_pas_gene_num, replace=False)
        diff_apa_tes = multi_pas_te_list.sample(diff_apa_te_num, replace=False)
        non_diff_apa_tes = multi_pas_te_list.drop(diff_apa_tes.index).sample(non_diff_apa_te_num, replace=False)
        
        # sample pas
        single_pas_gene_pas = pas_df_gene[pas_df_gene['gene_id'].isin(single_pas_genes)]
        diff_apa_te_pas = pas_df_gene[pas_df_gene['exon_id'].isin(diff_apa_tes)]
        diff_apa_te_pas = (diff_apa_te_pas.groupby('exon_id').
            apply(lambda x: x if len(x) <= MAX_PAS_PER_TE else x.sample(np.random.randint(2, MAX_PAS_PER_TE+1)).reset_index(drop=True))
            )
        non_diff_apa_te_pas = pas_df_gene[pas_df_gene['exon_id'].isin(non_diff_apa_tes)]
        non_diff_apa_te_pas = (non_diff_apa_te_pas.groupby('exon_id').
            apply(lambda x: x if len(x) <= MAX_PAS_PER_TE else x.sample(np.random.randint(2, MAX_PAS_PER_TE+1)).reset_index(drop=True))
            )
        
        single_pas_gene_pas.loc[:, "pas_type"] = "single_pas"
        diff_apa_te_pas.loc[:, "pas_type"] = "diff_apa"
        non_diff_apa_te_pas.loc[:, "pas_type"] = "non_diff_apa"
        
        pas_sample = pd.concat([single_pas_gene_pas, diff_apa_te_pas, non_diff_apa_te_pas]).reset_index(drop=True)
        pas_sample = pas_sample.sort_values(by=["chr", "start"])
        pas_sample = pas_sample.loc[:, ["chr", "start", "end", "gene_id", "score", "strand", "pas_type", "exon_id"]]
        # break
        pas_sample.to_csv(sample_bed, header=True, index=False, sep="\t")