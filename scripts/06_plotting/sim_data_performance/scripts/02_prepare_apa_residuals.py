#!/usr/bin/env python3
"""Prepare APA residual intermediate tables from DE-APA performance outputs."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import PROTOCOL_MAP, TOOL_MAP
from _shared.io import read_table, write_json, write_parquet
from _shared.paths import resolve_project_path, topic_intermediate_dir

TOPIC = "sim_data_performance"
DEFAULT_MAIN_PROJECT_ROOT = resolve_project_path("")
DE_APA_SUFFIX = "_de_apa_performance.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate DE-APA performance tables and compute per-group F1 residuals."
    )
    parser.add_argument(
        "--data-root",
        default=str(DEFAULT_MAIN_PROJECT_ROOT),
        help=(
            "Path to benchmark project root (contains data/) or data directory. "
            "Default: repository checkout root."
        ),
    )
    parser.add_argument(
        "--filter-type-1-suffix",
        default="0.05",
        help=(
            "Keep rows where filter_type_1 ends with this suffix before residual "
            "regression. Use empty string to disable filtering."
        ),
    )
    parser.add_argument(
        "--gt-count-total",
        type=float,
        default=5000.0,
        help="Ground-truth total used for gt_recall = gt_count / gt_count_total.",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="Optional debug cap on number of source files (0 means all files).",
    )
    parser.add_argument(
        "--f1-lower-quantile",
        type=float,
        default=0.01,
        help="Lower quantile bound for F1 trimming before regression (default: 0.01).",
    )
    parser.add_argument(
        "--f1-upper-quantile",
        type=float,
        default=0.99,
        help="Upper quantile bound for F1 trimming before regression (default: 0.99).",
    )
    return parser.parse_args()


def resolve_data_paths(data_root_arg: str) -> tuple[Path, Path]:
    root = Path(data_root_arg).expanduser().resolve()
    as_project_data = root / "data"
    as_data_dir = root

    project_perf = as_project_data / "result" / "performance"
    data_perf = as_data_dir / "result" / "performance"

    if project_perf.exists():
        return root, project_perf
    if data_perf.exists():
        return as_data_dir.parent, data_perf

    raise FileNotFoundError(
        "Cannot resolve performance directory from --data-root. "
        f"Tried: {project_perf} and {data_perf}"
    )


def _safe_relative(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _fill_identity_columns(df: pd.DataFrame, *, tool: str, sample: str) -> pd.DataFrame:
    if "tool" in df.columns:
        df["tool"] = df["tool"].fillna(tool).replace("", tool)
    else:
        df["tool"] = tool

    if "sample" in df.columns:
        df["sample"] = df["sample"].fillna(sample).replace("", sample)
    else:
        df["sample"] = sample

    default_protocol = sample.split("_")[0] if sample else ""
    if "protocol" in df.columns:
        df["protocol"] = df["protocol"].fillna(default_protocol).replace("", default_protocol)
    else:
        df["protocol"] = default_protocol
    return df


def _fit_linear_with_residuals(x: pd.Series, y: pd.Series) -> tuple[float, float, float, pd.Series]:
    x_num = pd.to_numeric(x, errors="coerce")
    y_num = pd.to_numeric(y, errors="coerce")
    valid = x_num.notna() & y_num.notna()

    residuals = pd.Series(np.nan, index=x.index, dtype=float)
    if int(valid.sum()) < 2:
        return np.nan, np.nan, np.nan, residuals

    xv = x_num.loc[valid].to_numpy(dtype=float)
    yv = y_num.loc[valid].to_numpy(dtype=float)
    design = np.column_stack((np.ones(len(xv)), xv))
    coef = np.linalg.lstsq(design, yv, rcond=None)[0]
    y_hat = design @ coef
    residuals.loc[valid] = yv - y_hat

    ss_res = float(np.sum((yv - y_hat) ** 2))
    ss_tot = float(np.sum((yv - np.mean(yv)) ** 2))
    r2 = np.nan if ss_tot == 0 else 1.0 - (ss_res / ss_tot)
    intercept = float(coef[0])
    slope = float(coef[1])
    return intercept, slope, r2, residuals


def _trim_group_by_quantile(
    group: pd.DataFrame, column: str, *, lower_q: float, upper_q: float
) -> pd.DataFrame:
    if column not in group.columns:
        return group.iloc[0:0].copy()
    out = group.copy()
    out[column] = pd.to_numeric(out[column], errors="coerce")
    valid = out[column].notna()
    if not bool(valid.any()):
        return out.iloc[0:0].copy()
    low = float(out.loc[valid, column].quantile(lower_q))
    high = float(out.loc[valid, column].quantile(upper_q))
    if pd.isna(low) or pd.isna(high):
        return out.iloc[0:0].copy()
    if low > high:
        low, high = high, low
    keep = valid & out[column].between(low, high, inclusive="both")
    return out.loc[keep].copy()


def _compute_group_residuals(group: pd.DataFrame, *, lower_q: float, upper_q: float) -> pd.DataFrame:
    current = _trim_group_by_quantile(group, "f1", lower_q=lower_q, upper_q=upper_q)
    current["gt_recall"] = pd.to_numeric(current["gt_recall"], errors="coerce")

    f1_intercept, f1_coef, f1_r2, f1_residuals = _fit_linear_with_residuals(
        current["gt_recall"], current["f1"]
    )
    current["f1_intercept"] = f1_intercept
    current["f1_coefficient"] = f1_coef
    current["f1_r2"] = f1_r2
    current["f1_residuals"] = f1_residuals

    if "rank" in current.columns:
        rank_intercept, rank_coef, rank_r2, rank_residuals = _fit_linear_with_residuals(
            current["gt_recall"], current["rank"]
        )
        current["rank_intercept"] = rank_intercept
        current["rank_coefficient"] = rank_coef
        current["rank_r2"] = rank_r2
        current["rank_residuals"] = rank_residuals
    else:
        current["rank_intercept"] = np.nan
        current["rank_coefficient"] = np.nan
        current["rank_r2"] = np.nan
        current["rank_residuals"] = np.nan

    return current


def main() -> None:
    args = parse_args()
    if not (0.0 <= float(args.f1_lower_quantile) < float(args.f1_upper_quantile) <= 1.0):
        raise ValueError("Require 0 <= --f1-lower-quantile < --f1-upper-quantile <= 1.")
    data_project_root, performance_dir = resolve_data_paths(args.data_root)

    de_apa_files = sorted(performance_dir.glob(f"*/*{DE_APA_SUFFIX}"))
    if args.max_files and args.max_files > 0:
        de_apa_files = de_apa_files[: args.max_files]
    if not de_apa_files:
        raise FileNotFoundError(f"No DE-APA performance files found in: {performance_dir}")

    frames: list[pd.DataFrame] = []
    for idx, table_path in enumerate(de_apa_files, start=1):
        tool_name = table_path.parent.name
        sample_name = table_path.name.removesuffix(DE_APA_SUFFIX)
        df = read_table(table_path)
        df = _fill_identity_columns(df, tool=tool_name, sample=sample_name)
        df["source_file"] = _safe_relative(table_path, data_project_root)
        frames.append(df)
        if idx % 2000 == 0:
            print(f"Loaded {idx}/{len(de_apa_files)} DE-APA tables ...")

    merged = pd.concat(frames, ignore_index=True, sort=False)
    for col in ["f1", "gt_count", "precision", "recall", "pd_count"]:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")

    if "tool" in merged.columns:
        merged["tool"] = merged["tool"].map(TOOL_MAP).fillna(merged["tool"])
    if "protocol" in merged.columns:
        merged["protocol"] = merged["protocol"].map(PROTOCOL_MAP).fillna(merged["protocol"])

    if "gt_count" not in merged.columns:
        raise KeyError("Missing required column 'gt_count' for residual calculation.")
    if "f1" not in merged.columns:
        raise KeyError("Missing required column 'f1' for residual calculation.")

    merged["gt_recall"] = merged["gt_count"] / float(args.gt_count_total)
    if {"filter_type_1", "filter_type_2"}.issubset(merged.columns):
        merged["filter_type_1"] = (
            merged["filter_type_1"]
            .astype("string")
            .str.strip()
            .replace({"": pd.NA, "nan": pd.NA, "None": pd.NA, "<NA>": pd.NA})
        )
        merged["filter_type_2"] = (
            merged["filter_type_2"]
            .astype("string")
            .str.strip()
            .replace({"": pd.NA, "nan": pd.NA, "None": pd.NA, "<NA>": pd.NA})
        )
        valid_filter = merged["filter_type_1"].notna() & merged["filter_type_2"].notna()
        merged["filter_type_comb"] = pd.Series(pd.NA, index=merged.index, dtype="string")
        merged.loc[valid_filter, "filter_type_comb"] = (
            merged.loc[valid_filter, "filter_type_1"]
            + "|"
            + merged.loc[valid_filter, "filter_type_2"]
        )

    residual_input = merged.copy()
    if args.filter_type_1_suffix and "filter_type_1" in residual_input.columns:
        suffix = str(args.filter_type_1_suffix)
        residual_input = residual_input[
            residual_input["filter_type_1"].astype("string").str.endswith(suffix, na=False)
        ].copy()

    required_group_cols = ["sample", "tool", "match_type"]
    missing_group_cols = [name for name in required_group_cols if name not in residual_input.columns]
    if missing_group_cols:
        raise KeyError(f"Missing required grouping columns: {missing_group_cols}")

    processed_groups: list[pd.DataFrame] = []
    grouped = residual_input.groupby(required_group_cols, dropna=False, sort=False)
    for idx, (_, group_df) in enumerate(grouped, start=1):
        processed_groups.append(
            _compute_group_residuals(
                group_df,
                lower_q=float(args.f1_lower_quantile),
                upper_q=float(args.f1_upper_quantile),
            )
        )
        if idx % 2000 == 0:
            print(f"Computed residuals for {idx} groups ...")

    if processed_groups:
        final_df = pd.concat(processed_groups, ignore_index=True, sort=False)
    else:
        final_df = residual_input.copy()
        for col in [
            "f1_intercept",
            "f1_coefficient",
            "f1_r2",
            "f1_residuals",
            "rank_intercept",
            "rank_coefficient",
            "rank_r2",
            "rank_residuals",
        ]:
            final_df[col] = np.nan

    summary_cols = [name for name in ["f1_residuals", "rank_residuals", "gt_recall"] if name in final_df.columns]
    if "filter_type_comb" in final_df.columns and summary_cols:
        by_filter_input = final_df[final_df["filter_type_comb"].notna()].copy()
        residual_by_filter = (
            by_filter_input.groupby("filter_type_comb", dropna=False)[summary_cols]
            .mean(numeric_only=True)
            .reset_index()
        )
        split = residual_by_filter["filter_type_comb"].astype("string").str.split(
            "|",
            n=1,
            expand=True,
            regex=False,
        )
        if split.shape[1] == 2:
            residual_by_filter["filter_type_1"] = split.iloc[:, 0]
            residual_by_filter["filter_type_2"] = split.iloc[:, 1]
        residual_by_filter = residual_by_filter.dropna(subset=["filter_type_1", "filter_type_2"]).copy()
        residual_by_filter = residual_by_filter[
            (residual_by_filter["filter_type_1"].astype(str).str.strip() != "")
            & (residual_by_filter["filter_type_2"].astype(str).str.strip() != "")
        ].copy()
    else:
        residual_by_filter = pd.DataFrame()

    sort_cols = [name for name in ["tool", "sample", "match_type", "filter_type_1", "filter_type_2"] if name in final_df.columns]
    if sort_cols:
        final_df = final_df.sort_values(sort_cols).reset_index(drop=True)

    out_dir = topic_intermediate_dir(TOPIC)
    final_path = out_dir / "sim_data_apa_residuals.parquet"
    by_filter_path = out_dir / "sim_data_apa_residuals_by_filter.parquet"
    metadata_path = out_dir / "sim_data_apa_residuals_metadata.json"

    write_parquet(final_df, final_path)
    write_parquet(residual_by_filter, by_filter_path)
    write_json(
        {
            "topic": TOPIC,
            "data_root_argument": str(args.data_root),
            "resolved_data_project_root": str(data_project_root),
            "resolved_performance_dir": str(performance_dir),
            "source_file_count": int(len(de_apa_files)),
            "filter_type_1_suffix": args.filter_type_1_suffix,
            "gt_count_total": float(args.gt_count_total),
            "f1_lower_quantile": float(args.f1_lower_quantile),
            "f1_upper_quantile": float(args.f1_upper_quantile),
            "rows_merged": int(len(merged)),
            "rows_residual_input": int(len(residual_input)),
            "rows_final": int(len(final_df)),
            "rows_by_filter": int(len(residual_by_filter)),
            "columns_final": list(final_df.columns),
            "columns_by_filter": list(residual_by_filter.columns),
        },
        metadata_path,
    )
    print(f"Wrote residual table: {final_path}")
    print(f"Wrote residual-by-filter table: {by_filter_path}")


if __name__ == "__main__":
    main()
