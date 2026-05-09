#!/usr/bin/env python3
"""Prepare DE-APA performance joined with peak-shape descriptors."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys
from typing import Any

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_SCRIPT_DIR = SCRIPT_DIR.parent
PLOT_ROOT = SIM_SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))
if str(SIM_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_SCRIPT_DIR))

from _shared.io import ensure_dir, read_table, write_json, write_parquet, write_tsv  # noqa: E402
from _shared.paths import topic_intermediate_dir  # noqa: E402
from _plot_helpers import _pick_apa_filter_set, _prepare_filter_comb  # noqa: E402

SIM_TOPIC = "sim_data_performance"
PEAK_TOPIC = "peak_model_params"
APA_RESIDUAL_INPUT = "sim_data_apa_residuals.parquet"
SINGLE_FILTER_INPUT = "sim_data_differential_apa_single_filter_performance.parquet"
PEAK_PARAMS_INPUT = "peak_model_params_prepared.tsv"

OUTPUT_DIR = SCRIPT_DIR / "output"
OUTPUT_TABLE = "de_apa_descriptor_performance.parquet"
OUTPUT_METADATA = "de_apa_descriptor_performance_metadata.json"
OUTPUT_UNMATCHED = "unmatched_sim_samples.tsv"
OUTPUT_DESCRIPTOR_SUMMARY = "descriptor_summary.tsv"
SINGLE_FILTER_OUTPUT_TABLE = "de_apa_single_filter_descriptor_performance.parquet"
SINGLE_FILTER_OUTPUT_METADATA = "de_apa_single_filter_descriptor_performance_metadata.json"
SINGLE_FILTER_OUTPUT_UNMATCHED = "single_filter_unmatched_sim_samples.tsv"
SINGLE_FILTER_OUTPUT_DESCRIPTOR_SUMMARY = "single_filter_descriptor_summary.tsv"

PERFORMANCE_METRICS = ["precision", "recall", "f1"]
SINGLE_FILTER_TYPE_1 = "scmapa_adPval_0.05"
SINGLE_FILTER_TYPE_2 = "scmapa_or_0.1"
EXPECTED_SINGLE_FILTER_TOOLS = ["scAPAtrap", "SCAPE", "Infernape", "SCAPTURE", "scAPA", "Sierra", "DaPars2", "scMAPA"]
SINGLE_FILTER_MATCH_TYPE_HARMONIZATION = {
    "Sample_1_vs_Sample_2": "T1_vs_T2",
}
DESCRIPTOR_COLUMNS = [
    "shape_per_peak_upstream_asymmetry_log2_mean",
    "shape_per_peak_effective_width_q10_q90_mean",
    "shape_per_peak_average_distance_mean",
]
DESCRIPTOR_LABELS = {
    "shape_per_peak_upstream_asymmetry_log2_mean": "Upstream asymmetry (log2)",
    "shape_per_peak_effective_width_q10_q90_mean": "Effective width (q10-q90)",
    "shape_per_peak_average_distance_mean": "Average distance",
}
SIM_SAMPLE_RE = re.compile(
    r"^(?P<empirical_sample_id>.+)_(?P<genome>hg38|mm10)_(?P<pas_scenario>pas[^_]+)_gn(?P<gene_count>\d+)_rep(?P<sim_rep>\d+)$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--apa-residuals",
        type=Path,
        default=topic_intermediate_dir(SIM_TOPIC) / APA_RESIDUAL_INPUT,
        help="DE-APA residual/performance table produced by 02_prepare_apa_residuals.py.",
    )
    parser.add_argument(
        "--single-filter-performance",
        type=Path,
        default=topic_intermediate_dir(SIM_TOPIC) / SINGLE_FILTER_INPUT,
        help=(
            "DaPars2/scMAPA-compatible single-filter table produced by "
            "03_prepare_differential_apa_single_filter_performance.py."
        ),
    )
    parser.add_argument(
        "--peak-params",
        type=Path,
        default=topic_intermediate_dir(PEAK_TOPIC) / PEAK_PARAMS_INPUT,
        help="Peak descriptor table produced by peak_model_params/00_prepare_data.py.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Directory for joined table and QC outputs.",
    )
    parser.add_argument(
        "--allow-unmatched",
        action="store_true",
        help="Write output even when some sim empirical samples cannot be matched to peak descriptors.",
    )
    return parser.parse_args()


def parse_sim_sample(sample: str) -> dict[str, Any]:
    match = SIM_SAMPLE_RE.match(str(sample))
    if not match:
        return {
            "empirical_sample_id": str(sample),
            "genome": pd.NA,
            "pas_scenario": pd.NA,
            "gene_count": pd.NA,
            "sim_rep": pd.NA,
            "sample_parse_ok": False,
        }
    data: dict[str, Any] = match.groupdict()
    data["gene_count"] = int(data["gene_count"])
    data["sim_rep"] = int(data["sim_rep"])
    data["sample_parse_ok"] = True
    return data


def load_apa_residuals(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing APA residual table: {path}. Run 02_prepare_apa_residuals.py first.")
    usecols = [
        "sample",
        "tool",
        "protocol",
        "match_type",
        "filter_type_1",
        "filter_type_2",
        "filter_type_comb",
        "precision",
        "recall",
        "f1",
        "gt_recall",
        "f1_residuals",
    ]
    df = pd.read_parquet(path, columns=usecols)
    return _prepare_filter_comb(df)


def select_main_criteria(df: pd.DataFrame) -> list[str]:
    selected = _pick_apa_filter_set(df)
    if not selected:
        raise ValueError("Could not select main-text APA filter combinations from residual table.")
    return selected


def collapse_main_criteria(df: pd.DataFrame, selected_criteria: list[str]) -> pd.DataFrame:
    work = df[df["filter_type_comb"].astype(str).isin(selected_criteria)].copy()
    if work.empty:
        raise ValueError("No APA residual rows remain after selected criteria filtering.")

    work = add_sim_sample_metadata(work)

    group_cols = [
        "tool",
        "sample",
        "empirical_sample_id",
        "protocol",
        "match_type",
        "genome",
        "pas_scenario",
        "gene_count",
        "sim_rep",
        "sample_parse_ok",
    ]
    grouped = (
        work.groupby(group_cols, dropna=False, as_index=False)
        .agg(
            precision=("precision", "mean"),
            recall=("recall", "mean"),
            f1=("f1", "mean"),
            gt_recall=("gt_recall", "mean"),
            n_filter_combinations=("filter_type_comb", "nunique"),
            filter_type_comb=("filter_type_comb", lambda values: "|;|".join(sorted(set(map(str, values))))),
        )
        .reset_index(drop=True)
    )
    return grouped


def add_sim_sample_metadata(df: pd.DataFrame) -> pd.DataFrame:
    if "sample" not in df.columns:
        raise KeyError("Missing required column 'sample' for sim sample parsing.")
    out = df.copy()
    sample_meta = pd.DataFrame([parse_sim_sample(sample) for sample in out["sample"].astype(str)])
    sample_meta.index = out.index
    return pd.concat([out, sample_meta], axis=1)


def load_single_filter_performance(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Missing single-filter DE-APA performance table: {path}. "
            "Run 03_prepare_differential_apa_single_filter_performance.py first."
        )
    required = {
        "sample",
        "tool",
        "protocol",
        "match_type",
        "filter_type_1",
        "filter_type_2",
        "precision",
        "recall",
        "f1",
    }
    df = read_table(path)
    missing = sorted(required.difference(df.columns))
    if missing:
        raise KeyError(f"Missing required single-filter performance columns: {missing}")

    work = df.copy()
    work["filter_type_1"] = work["filter_type_1"].astype("string")
    work["filter_type_2"] = work["filter_type_2"].astype("string")
    invalid_filter = work[
        (work["filter_type_1"] != SINGLE_FILTER_TYPE_1)
        | (work["filter_type_2"] != SINGLE_FILTER_TYPE_2)
    ]
    if not invalid_filter.empty:
        observed = (
            invalid_filter[["filter_type_1", "filter_type_2"]]
            .drop_duplicates()
            .astype(str)
            .to_dict("records")
        )
        raise ValueError(
            "Single-filter performance table contains rows outside the expected "
            f"{SINGLE_FILTER_TYPE_1}|{SINGLE_FILTER_TYPE_2} criteria pair: {observed[:10]}"
        )

    for metric in PERFORMANCE_METRICS:
        work[metric] = pd.to_numeric(work[metric], errors="coerce")
    if "gt_count" in work.columns:
        work["gt_count"] = pd.to_numeric(work["gt_count"], errors="coerce")
    if "pd_count" in work.columns:
        work["pd_count"] = pd.to_numeric(work["pd_count"], errors="coerce")

    keep_cols = [
        col
        for col in [
            "tool",
            "sample",
            "protocol",
            "match_type",
            "filter_type_1",
            "filter_type_2",
            "precision",
            "recall",
            "f1",
            "gt_count",
            "pd_count",
            "source_kind",
            "source_file",
            "tool_raw",
        ]
        if col in work.columns
    ]
    out = work[keep_cols].copy()
    out = add_sim_sample_metadata(out)
    out["filter_type_comb"] = SINGLE_FILTER_TYPE_1 + "|" + SINGLE_FILTER_TYPE_2
    out["match_type_harmonized"] = (
        out["match_type"].astype(str).replace(SINGLE_FILTER_MATCH_TYPE_HARMONIZATION).astype("string")
    )
    return out


def load_peak_descriptors(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing peak descriptor table: {path}. Run peak_model_params/00_prepare_data.py first.")
    required = ["sample_id", "sample_full_id", "protocol", "species", "tissue", *DESCRIPTOR_COLUMNS]
    df = read_table(path)
    missing = sorted(set(required).difference(df.columns))
    if missing:
        raise KeyError(f"Missing required peak descriptor columns: {missing}")

    out = df[required].copy()
    out = out.rename(
        columns={
            "sample_id": "peak_sample_id",
            "sample_full_id": "peak_sample_full_id",
            "protocol": "peak_protocol",
        }
    )
    for col in DESCRIPTOR_COLUMNS:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    return out.drop_duplicates(subset=["peak_sample_id"], keep="first")


def add_descriptor_zscores(joined: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    out = joined.copy()
    matched_samples = out.loc[out["descriptor_matched"], ["empirical_sample_id", *DESCRIPTOR_COLUMNS]].drop_duplicates(
        subset=["empirical_sample_id"]
    )
    summary_rows: list[dict[str, Any]] = []
    for col in DESCRIPTOR_COLUMNS:
        values = matched_samples[col].dropna()
        mean = float(values.mean()) if not values.empty else np.nan
        sd = float(values.std(ddof=0)) if len(values) > 1 else np.nan
        z_col = f"{col}_z"
        if pd.notna(sd) and sd > 0:
            out[z_col] = (pd.to_numeric(out[col], errors="coerce") - mean) / sd
        else:
            out[z_col] = np.nan
        summary_rows.append(
            {
                "descriptor": col,
                "label": DESCRIPTOR_LABELS[col],
                "n_matched_samples": int(values.shape[0]),
                "missing_matched_samples": int(matched_samples[col].isna().sum()) if not matched_samples.empty else 0,
                "mean": mean,
                "sd": sd,
                "min": float(values.min()) if not values.empty else np.nan,
                "max": float(values.max()) if not values.empty else np.nan,
                "zscore_column": z_col,
                "zscore_scope": "unique matched empirical samples",
            }
        )
    return out, summary_rows


def join_peak_descriptors(performance: pd.DataFrame, peak: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, Any]]]:
    joined = performance.merge(
        peak,
        left_on="empirical_sample_id",
        right_on="peak_sample_id",
        how="left",
        validate="many_to_one",
    )
    joined["descriptor_matched"] = joined["peak_sample_id"].notna()
    joined["missing_descriptor_reason"] = np.where(
        joined["descriptor_matched"],
        "",
        np.where(joined["sample_parse_ok"], "no_peak_sample_id_match", "sample_parse_failed"),
    )
    joined, descriptor_summary = add_descriptor_zscores(joined)
    unmatched = build_unmatched_table(joined)
    return joined, unmatched, descriptor_summary


def build_unmatched_table(joined: pd.DataFrame) -> pd.DataFrame:
    unmatched = joined[~joined["descriptor_matched"]].copy()
    if unmatched.empty:
        return pd.DataFrame(
            columns=[
                "empirical_sample_id",
                "n_sim_samples",
                "n_rows",
                "protocols",
                "sample_parse_ok",
                "example_sample",
                "missing_descriptor_reason",
            ]
        )
    return (
        unmatched.groupby("empirical_sample_id", dropna=False)
        .agg(
            n_sim_samples=("sample", "nunique"),
            n_rows=("sample", "size"),
            protocols=("protocol", lambda values: "|".join(sorted(set(map(str, values))))),
            sample_parse_ok=("sample_parse_ok", "all"),
            example_sample=("sample", "first"),
            missing_descriptor_reason=("missing_descriptor_reason", "first"),
        )
        .reset_index()
        .sort_values(["protocols", "empirical_sample_id"], kind="stable")
    )


def tool_counts(df: pd.DataFrame) -> dict[str, int]:
    return {str(k): int(v) for k, v in df["tool"].value_counts().sort_index().to_dict().items()}


def write_descriptor_outputs(
    *,
    joined: pd.DataFrame,
    unmatched: pd.DataFrame,
    descriptor_summary: list[dict[str, Any]],
    out_dir: Path,
    output_table: str,
    output_metadata: str,
    output_unmatched: str,
    output_descriptor_summary: str,
    metadata: dict[str, Any],
    allow_unmatched: bool,
) -> None:
    unmatched_path = write_tsv(unmatched, out_dir / output_unmatched)
    if not allow_unmatched and not unmatched.empty:
        raise ValueError(
            f"Found {len(unmatched)} unmatched empirical samples. "
            f"Wrote details to {unmatched_path}. Rerun with --allow-unmatched to write the joined table anyway."
        )

    output_path = write_parquet(joined, out_dir / output_table)
    descriptor_summary_path = write_tsv(pd.DataFrame(descriptor_summary), out_dir / output_descriptor_summary)
    metadata = {
        **metadata,
        "output_table": str(output_path),
        "descriptor_columns": DESCRIPTOR_COLUMNS,
        "descriptor_labels": DESCRIPTOR_LABELS,
        "rows_output": int(len(joined)),
        "tool_counts": tool_counts(joined),
        "n_tools": int(joined["tool"].dropna().astype(str).nunique()),
        "n_match_types": int(joined["match_type"].dropna().astype(str).nunique()),
        "n_sim_samples": int(joined["sample"].dropna().astype(str).nunique()),
        "n_empirical_samples": int(joined["empirical_sample_id"].dropna().astype(str).nunique()),
        "n_unmatched_empirical_samples": int(len(unmatched)),
        "descriptor_match_rate_rows": float(joined["descriptor_matched"].mean()) if len(joined) else np.nan,
        "performance_metrics": PERFORMANCE_METRICS,
        "zscore_scope": "unique matched empirical samples",
        "unmatched_table": str(unmatched_path),
        "descriptor_summary": str(descriptor_summary_path),
    }
    write_json(metadata, out_dir / output_metadata)

    print(f"Wrote joined table: {output_path}")
    print(f"Wrote metadata: {out_dir / output_metadata}")
    print(f"Wrote unmatched samples: {unmatched_path}")
    print(f"Wrote descriptor summary: {descriptor_summary_path}")


def main() -> None:
    args = parse_args()
    out_dir = ensure_dir(args.output_dir)

    apa = load_apa_residuals(args.apa_residuals)
    selected_criteria = select_main_criteria(apa)
    collapsed = collapse_main_criteria(apa, selected_criteria)

    peak = load_peak_descriptors(args.peak_params)
    joined, unmatched, descriptor_summary = join_peak_descriptors(collapsed, peak)
    write_descriptor_outputs(
        joined=joined,
        unmatched=unmatched,
        descriptor_summary=descriptor_summary,
        out_dir=out_dir,
        output_table=OUTPUT_TABLE,
        output_metadata=OUTPUT_METADATA,
        output_unmatched=OUTPUT_UNMATCHED,
        output_descriptor_summary=OUTPUT_DESCRIPTOR_SUMMARY,
        allow_unmatched=args.allow_unmatched,
        metadata={
            "analysis_regime": "main_text_three_criteria",
            "apa_residuals": str(args.apa_residuals),
            "peak_params": str(args.peak_params),
            "selected_filter_combinations": selected_criteria,
            "criteria_merge_rule": "mean precision/recall/f1 within tool x sample x match_type across selected criteria",
            "rows_input_apa": int(len(apa)),
            "rows_collapsed_performance": int(len(collapsed)),
            "n_filter_combinations_per_row": {
                str(key): int(value)
                for key, value in joined["n_filter_combinations"].value_counts(dropna=False).sort_index().items()
            },
            "n_peak_samples": int(peak["peak_sample_id"].dropna().astype(str).nunique()),
        },
    )

    single = load_single_filter_performance(args.single_filter_performance)
    single_joined, single_unmatched, single_descriptor_summary = join_peak_descriptors(single, peak)
    present_tools = set(single_joined["tool"].dropna().astype(str))
    missing_expected_tools = [tool for tool in EXPECTED_SINGLE_FILTER_TOOLS if tool not in present_tools]
    source_kind_counts = (
        single_joined.groupby(["tool", "source_kind"], dropna=False)
        .size()
        .reset_index(name="rows")
        .to_dict("records")
        if "source_kind" in single_joined.columns
        else []
    )
    write_descriptor_outputs(
        joined=single_joined,
        unmatched=single_unmatched,
        descriptor_summary=single_descriptor_summary,
        out_dir=out_dir,
        output_table=SINGLE_FILTER_OUTPUT_TABLE,
        output_metadata=SINGLE_FILTER_OUTPUT_METADATA,
        output_unmatched=SINGLE_FILTER_OUTPUT_UNMATCHED,
        output_descriptor_summary=SINGLE_FILTER_OUTPUT_DESCRIPTOR_SUMMARY,
        allow_unmatched=args.allow_unmatched,
        metadata={
            "analysis_regime": "dapars2_scmapa_compatible_single_filter",
            "single_filter_performance": str(args.single_filter_performance),
            "peak_params": str(args.peak_params),
            "filter_type_1": SINGLE_FILTER_TYPE_1,
            "filter_type_2": SINGLE_FILTER_TYPE_2,
            "criteria_merge_rule": "no criteria averaging; table is already restricted to one common filter pair",
            "match_type_harmonization": SINGLE_FILTER_MATCH_TYPE_HARMONIZATION,
            "instance_match_type_column": "match_type_harmonized",
            "expected_tools": EXPECTED_SINGLE_FILTER_TOOLS,
            "missing_expected_tools": missing_expected_tools,
            "all_expected_tools_present": not missing_expected_tools,
            "source_selection_policy": (
                "non-dapars2/scmapa use *_de_apa_dapars_parallel_performance.tsv; "
                "dapars2/scmapa use *_de_apa_performance.tsv, as prepared by "
                "03_prepare_differential_apa_single_filter_performance.py"
            ),
            "source_kind_counts": source_kind_counts,
            "rows_input_single_filter": int(len(single)),
            "n_peak_samples": int(peak["peak_sample_id"].dropna().astype(str).nunique()),
        },
    )


if __name__ == "__main__":
    main()
