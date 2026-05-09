#!/usr/bin/env python3
"""Prepare recall/gt-only table using configured filter pairs and gt threshold."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.io import read_table, write_json, write_parquet, write_tsv
from _shared.paths import topic_intermediate_dir, topic_root

TOPIC = "raw_data_performance"
DEFAULT_INPUT_NAME = "raw_data_performance_prepared.parquet"
DEFAULT_CONFIG_PATH = topic_root(TOPIC) / "config" / "recall_gt_filter.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter raw-data metrics to recall/gt with configured filter pairs and gt rules."
    )
    parser.add_argument(
        "--input",
        default=str(topic_intermediate_dir(TOPIC) / DEFAULT_INPUT_NAME),
        help="Input merged metrics parquet from 00_prepare_data.py",
    )
    parser.add_argument(
        "--config",
        default=str(DEFAULT_CONFIG_PATH),
        help="JSON config for filter pairs and gt rules.",
    )
    return parser.parse_args()


def _load_config(config_path: Path) -> dict:
    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    metric_domains = cfg.get("metric_domains", [])
    filter_pairs = cfg.get("filter_pairs", [])
    gt_rules = cfg.get("gt_rules", {})
    recall_rules = cfg.get("recall_rules", {})

    if not metric_domains:
        raise ValueError("config.metric_domains is required")
    if not filter_pairs:
        raise ValueError("config.filter_pairs is required")

    min_gt_count = int(gt_rules.get("min_gt_count", 1))
    drop_gt_zero = bool(gt_rules.get("drop_gt_zero", True))
    drop_all_tool_zero = bool(recall_rules.get("drop_all_tool_zero_recall_sample_pairs", False))

    if min_gt_count < 0:
        raise ValueError("gt_rules.min_gt_count must be >= 0")

    return {
        "metric_domains": metric_domains,
        "filter_pairs": filter_pairs,
        "gt_rules": {
            "drop_gt_zero": drop_gt_zero,
            "min_gt_count": min_gt_count,
        },
        "recall_rules": {
            "drop_all_tool_zero_recall_sample_pairs": drop_all_tool_zero,
        },
    }


def _pair_key(filter_type_1: str, filter_type_2: str) -> str:
    return f"{filter_type_1} | {filter_type_2}"


def _distribution_stats(df: pd.DataFrame, min_gt_count: int, *, stage: str) -> pd.DataFrame:
    group_cols = ["filter_type_1", "filter_type_2"]
    rows = []
    for (f1, f2), sub in df.groupby(group_cols, dropna=False):
        gt = pd.to_numeric(sub["gt_count"], errors="coerce")
        rows.append(
            {
                "stage": stage,
                "filter_type_1": f1,
                "filter_type_2": f2,
                "filter_pair": _pair_key(str(f1), str(f2)),
                "rows": int(len(sub)),
                "gt_non_na": int(gt.notna().sum()),
                "gt_zero": int((gt == 0).sum()),
                "gt_lt_min": int((gt < min_gt_count).sum()),
                "gt_min": float(gt.min()) if gt.notna().any() else None,
                "gt_q25": float(gt.quantile(0.25)) if gt.notna().any() else None,
                "gt_median": float(gt.median()) if gt.notna().any() else None,
                "gt_q75": float(gt.quantile(0.75)) if gt.notna().any() else None,
                "gt_max": float(gt.max()) if gt.notna().any() else None,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).expanduser().resolve()
    config_path = Path(args.config).expanduser().resolve()

    if not input_path.exists():
        raise FileNotFoundError(f"Missing input metrics file: {input_path}")
    if not config_path.exists():
        raise FileNotFoundError(f"Missing config file: {config_path}")

    cfg = _load_config(config_path)
    metric_domains = cfg["metric_domains"]
    filter_pairs = cfg["filter_pairs"]
    min_gt_count = int(cfg["gt_rules"]["min_gt_count"])
    drop_gt_zero = bool(cfg["gt_rules"]["drop_gt_zero"])
    drop_all_tool_zero = bool(cfg["recall_rules"]["drop_all_tool_zero_recall_sample_pairs"])

    df = read_table(input_path)
    required_cols = {
        "metric_domain",
        "filter_type_1",
        "filter_type_2",
        "recall",
        "gt_count",
    }
    missing = sorted(required_cols - set(df.columns))
    if missing:
        raise ValueError(f"Input table missing required columns: {missing}")

    pair_set = {
        (str(item["filter_type_1"]), str(item["filter_type_2"]))
        for item in filter_pairs
    }

    sub = df[df["metric_domain"].isin(metric_domains)].copy()
    sub = sub[
        sub.apply(
            lambda x: (str(x["filter_type_1"]), str(x["filter_type_2"])) in pair_set,
            axis=1,
        )
    ].copy()

    sub["recall"] = pd.to_numeric(sub["recall"], errors="coerce")
    sub["gt_count"] = pd.to_numeric(sub["gt_count"], errors="coerce")
    sub = sub[sub["recall"].notna() & sub["gt_count"].notna()].copy()

    before = sub.copy()

    if drop_gt_zero:
        sub = sub[sub["gt_count"] != 0].copy()
    sub = sub[sub["gt_count"] >= min_gt_count].copy()

    zero_pair_rows_removed = 0
    zero_pair_count_removed = 0
    if drop_all_tool_zero:
        key_candidates = ["filter_type_1", "filter_type_2", "sample_id", "pair_label"]
        group_cols = [col for col in key_candidates if col in sub.columns]
        if {"sample_id", "pair_label"}.issubset(set(group_cols)):
            per_pair = (
                sub.groupby(group_cols, dropna=False)["recall"]
                .max()
                .reset_index(name="max_recall")
            )
            zero_pairs = per_pair[per_pair["max_recall"] == 0].copy()
            if not zero_pairs.empty:
                zero_pair_count_removed = int(len(zero_pairs))
                zero_keys = set(
                    tuple(row[col] for col in group_cols)
                    for _, row in zero_pairs.iterrows()
                )
                keep_mask = sub.apply(
                    lambda row: tuple(row[col] for col in group_cols) not in zero_keys,
                    axis=1,
                )
                zero_pair_rows_removed = int((~keep_mask).sum())
                sub = sub[keep_mask].copy()

    keep_cols = [
        c
        for c in [
            "metric_domain",
            "tool",
            "protocol",
            "sample",
            "sample_id",
            "pair_label",
            "match_type",
            "filter_type_1",
            "filter_type_2",
            "recall",
            "gt_count",
            "raw_run_id",
            "source_file",
        ]
        if c in sub.columns
    ]
    filtered = sub[keep_cols].copy()
    filtered["filter_pair"] = (
        filtered["filter_type_1"].astype(str)
        + " | "
        + filtered["filter_type_2"].astype(str)
    )

    dist_before = _distribution_stats(before, min_gt_count, stage="before_gt_filter")
    dist_after = _distribution_stats(filtered, min_gt_count, stage="after_gt_filter")
    dist = pd.concat([dist_before, dist_after], ignore_index=True)

    out_dir = topic_intermediate_dir(TOPIC)
    filtered_path = out_dir / "raw_data_recall_gt_filtered.parquet"
    dist_path = out_dir / "raw_data_recall_gt_distribution.tsv"
    metadata_path = out_dir / "raw_data_recall_gt_metadata.json"

    write_parquet(filtered, filtered_path)
    write_tsv(dist, dist_path)
    write_json(
        {
            "topic": TOPIC,
            "input": str(input_path),
            "config": str(config_path),
            "metric_domains": metric_domains,
            "filter_pairs": [
                {
                    "filter_type_1": f1,
                    "filter_type_2": f2,
                    "filter_pair": _pair_key(f1, f2),
                }
                for (f1, f2) in sorted(pair_set)
            ],
            "gt_rules": {
                "drop_gt_zero": drop_gt_zero,
                "min_gt_count": min_gt_count,
            },
            "recall_rules": {
                "drop_all_tool_zero_recall_sample_pairs": drop_all_tool_zero,
                "zero_pair_count_removed": zero_pair_count_removed,
                "zero_pair_rows_removed": zero_pair_rows_removed,
            },
            "rows": {
                "before_gt_filter": int(len(before)),
                "after_gt_filter": int(len(filtered)),
            },
            "columns": list(filtered.columns),
            "outputs": {
                "filtered": str(filtered_path),
                "distribution": str(dist_path),
            },
        },
        metadata_path,
    )

    print(f"Wrote filtered table: {filtered_path}")
    print(f"Wrote gt distribution: {dist_path}")
    print(f"Rows before gt filter: {len(before)}")
    print(f"Rows after gt filter: {len(filtered)}")


if __name__ == "__main__":
    main()
