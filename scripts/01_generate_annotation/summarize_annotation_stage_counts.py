from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
ANNOTATION_DIR = ROOT / "data" / "int_data" / "annotations"
SIM_PAS_DIR = ROOT / "data" / "sim_data" / "sim_pas"
OUTPUT_PATH = ANNOTATION_DIR / "supplementary_annotation_stage_counts.tsv"

SOURCE_COUNTS = {
    "Human": {
        "PolyA_DB v4.1 PAS after within-source collapse": 281_091,
        "PolyASite v3.0 PAS after within-source collapse": 18_432_135,
        "Gencode transcript-end PAS after collapse": 42_278,
    },
    "Mouse": {
        "PolyA_DB v4.1 PAS after within-source collapse": 251_177,
        "PolyASite v3.0 PAS after within-source collapse": 1_750_661,
        "Gencode transcript-end PAS after collapse": 27_180,
    },
}

STAGES = [
    "PolyA_DB v4.1 PAS after within-source collapse",
    "PolyASite v3.0 PAS after within-source collapse",
    "Gencode transcript-end PAS after collapse",
    "integrated raw unique PAS catalog",
    "assigned to Gencode terminal exons",
    "assigned unambiguously to one gene / TE",
    "candidate single-PAS genes",
    "candidate multi-PAS terminal exons",
    "sampled annotation sets",
]

COLUMNS = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "gene_id",
    "gene_name",
    "pas_type",
    "exon_id",
]


def read_integrated_pas(path):
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=COLUMNS,
        usecols=["chrom", "start", "end", "name", "strand", "gene_id", "pas_type", "exon_id"],
        dtype={
            "chrom": "string",
            "start": "int64",
            "end": "int64",
            "name": "string",
            "strand": "string",
            "gene_id": "string",
            "pas_type": "string",
            "exon_id": "string",
        },
    )


def count_candidate_sets(df):
    df = df[df["chrom"].str.contains(r"^chr[0-9XY]", regex=True, na=False)].copy()
    duplicated_names = df.loc[df.duplicated("name", keep=False), "name"].unique()
    df = df.loc[~df["name"].isin(duplicated_names)].copy()

    terminal_df = df[df["pas_type"].isin(["SE", "TE"])].copy()
    single_pas_genes = terminal_df.groupby("gene_id").filter(lambda x: len(x) == 1)[
        "gene_id"
    ].nunique()

    terminal_forward = terminal_df[terminal_df["strand"] == "+"].sort_values(["chrom", "start"])
    terminal_reverse = terminal_df[terminal_df["strand"] == "-"].sort_values(["chrom", "start"])
    terminal_forward["gap"] = (
        terminal_forward.groupby("gene_id")["start"].shift(-1) - terminal_forward["end"]
    )
    terminal_reverse["gap"] = (
        terminal_reverse.groupby("gene_id")["start"].shift(-1) - terminal_reverse["end"]
    )
    terminal_df = pd.concat([terminal_forward, terminal_reverse], ignore_index=True)
    terminal_df = terminal_df[
        (terminal_df["gap"] >= 200) | (terminal_df["gap"] < 0) | terminal_df["gap"].isna()
    ]
    multi_pas_terminal_exons = terminal_df.groupby("exon_id").filter(lambda x: len(x) > 1)[
        "exon_id"
    ].nunique()

    return single_pas_genes, multi_pas_terminal_exons


def summarize_species(label, path, sample_prefix):
    df = read_integrated_pas(path)

    counts = dict(SOURCE_COUNTS[label])
    counts["integrated raw unique PAS catalog"] = df["name"].nunique()

    terminal_df = df[df["pas_type"].isin(["SE", "TE"])].copy()
    counts["assigned to Gencode terminal exons"] = terminal_df["name"].nunique()

    assignment_counts = terminal_df.groupby("name").agg(
        gene_count=("gene_id", "nunique"),
        exon_count=("exon_id", "nunique"),
    )
    counts["assigned unambiguously to one gene / TE"] = len(
        assignment_counts[
            (assignment_counts["gene_count"] == 1) & (assignment_counts["exon_count"] == 1)
        ]
    )

    single_pas_genes, multi_pas_terminal_exons = count_candidate_sets(df)
    counts["candidate single-PAS genes"] = single_pas_genes
    counts["candidate multi-PAS terminal exons"] = multi_pas_terminal_exons

    counts["sampled annotation sets"] = len(
        list(SIM_PAS_DIR.glob(f"{sample_prefix}_sim_pas_*.bed"))
    )
    return counts


def main():
    human_counts = summarize_species(
        "Human", ANNOTATION_DIR / "human_integrated_pas.bed", "hg38"
    )
    mouse_counts = summarize_species(
        "Mouse", ANNOTATION_DIR / "mouse_integrated_pas.bed", "mm10"
    )

    table = pd.DataFrame(
        {
            "Stage": STAGES,
            "Human": [human_counts[stage] for stage in STAGES],
            "Mouse": [mouse_counts[stage] for stage in STAGES],
        }
    )
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUTPUT_PATH, sep="\t", index=False)
    print(table.to_string(index=False))
    print(f"Wrote {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
