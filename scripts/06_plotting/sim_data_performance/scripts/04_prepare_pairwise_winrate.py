#!/usr/bin/env python3
"""Prepare pairwise tool win-rate table for sim_data_performance heatmaps."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Iterable

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import TOOL_ORDER
from _shared.io import write_json, write_parquet
from _shared.paths import topic_intermediate_dir

from _plot_helpers import (
    TOPIC,
    build_apa_plot_table,
    build_apa_plot_table_by_filter,
    load_metric_domain,
    load_single_filter_apa_table,
    trim_metric_to_quantile,
)

OUT_TABLE = "sim_data_pairwise_winrate.parquet"
OUT_META = "sim_data_pairwise_winrate_metadata.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build pairwise win-rate table from existing sim_data_performance inputs."
    )
    parser.add_argument(
        "--epsilon-default",
        type=float,
        default=0.01,
        help="Default epsilon for win/tie/loss thresholding.",
    )
    parser.add_argument(
        "--epsilon-zero-metrics",
        default="",
        help="Comma-separated metric names that should use epsilon=0.",
    )
    parser.add_argument(
        "--balanced-by",
        choices=["none", "match_type"],
        default="none",
        help="Optional balanced win-rate by averaging per stratum.",
    )
    return parser.parse_args()


def _slug(text: str) -> str:
    out = []
    for ch in text.lower():
        if ch.isalnum():
            out.append(ch)
        else:
            out.append("_")
    s = "".join(out)
    while "__" in s:
        s = s.replace("__", "_")
    return s.strip("_")


def _resolve_tool_order(frames: Iterable[pd.DataFrame]) -> list[str]:
    observed: set[str] = set()
    for frame in frames:
        if "tool" not in frame.columns:
            continue
        observed.update(frame["tool"].dropna().astype(str).unique())
    ordered = [name for name in TOOL_ORDER if name in observed]
    leftovers = sorted(name for name in observed if name not in ordered)
    return ordered + leftovers


def _panel_win_rate(
    df: pd.DataFrame,
    *,
    dataset: str,
    metric: str,
    analysis_type: str,
    analysis_value: str,
    analysis_protocol: str | None,
    analysis_filter_comb: str | None,
    eps: float,
    tool_order: list[str],
    balanced_by: str,
    unit_cols_candidates: list[str] | None = None,
) -> list[dict[str, object]]:
    needed = {"tool", metric}
    if not needed.issubset(df.columns):
        return []

    if unit_cols_candidates is None:
        unit_cols_candidates = ["sample", "match_type"]
    unit_cols = [col for col in unit_cols_candidates if col in df.columns]
    if not unit_cols:
        return []

    work = trim_metric_to_quantile(df, metric)
    if work.empty:
        return []

    score = pd.to_numeric(work[metric], errors="coerce")
    if metric.startswith("mape"):
        score = -score
    work = work.copy()
    work["_score"] = score
    work = work.dropna(subset=["_score", "tool"] + unit_cols)
    if work.empty:
        return []

    wide = work.pivot_table(index=unit_cols, columns="tool", values="_score", aggfunc="mean")
    if wide.empty:
        return []

    present_tools = [name for name in tool_order if name in set(wide.columns.astype(str))]
    rows: list[dict[str, object]] = []
    panel_id = f"{dataset}__{metric}__{analysis_type}__{_slug(analysis_value)}"

    for tool_a in present_tools:
        for tool_b in present_tools:
            if tool_a == tool_b:
                continue
            pair = wide[[tool_a, tool_b]].dropna()
            if pair.empty:
                continue

            diff = pair[tool_a].to_numpy(dtype=float) - pair[tool_b].to_numpy(dtype=float)
            win = np.where(diff > eps, 1.0, np.where(np.abs(diff) <= eps, 0.5, 0.0))
            win_rate = float(np.mean(win))

            win_rate_balanced = np.nan
            if balanced_by != "none" and balanced_by in unit_cols:
                pair_reset = pair.reset_index()
                pair_reset["_win"] = win
                by_stratum = pair_reset.groupby(balanced_by, dropna=False)["_win"].mean()
                if len(by_stratum) > 0:
                    win_rate_balanced = float(by_stratum.mean())

            rows.append(
                {
                    "panel_id": panel_id,
                    "dataset": dataset,
                    "metric": metric,
                    "analysis_type": analysis_type,
                    "analysis_value": analysis_value,
                    "analysis_protocol": analysis_protocol,
                    "analysis_filter_comb": analysis_filter_comb,
                    "tool_a": tool_a,
                    "tool_b": tool_b,
                    "win_rate": win_rate,
                    "win_rate_balanced": win_rate_balanced,
                    "n_units": int(len(pair)),
                    "epsilon": float(eps),
                }
            )
    return rows


def main() -> None:
    args = parse_args()
    eps_zero_metrics = {
        item.strip() for item in str(args.epsilon_zero_metrics).split(",") if item.strip()
    }

    pas_detect = load_metric_domain("match_performance")
    pas_quant = load_metric_domain("pas_quantify_performance")
    apa_detect = build_apa_plot_table()
    apa_detect_by_filter = build_apa_plot_table_by_filter()
    apa_single = load_single_filter_apa_table()

    frames = [pas_detect, pas_quant, apa_detect, apa_detect_by_filter, apa_single]
    common_order = _resolve_tool_order(frames)

    datasets: list[tuple[str, pd.DataFrame, list[str], list[str], bool]] = [
        ("pas_detect", pas_detect, ["precision", "recall", "f1"], ["sample", "match_type"], False),
        (
            "pas_quantification",
            pas_quant,
            ["cor_pas", "mape_pas", "mape_pas_ct"],
            ["sample", "match_type"],
            False,
        ),
        ("apa_detect", apa_detect, ["precision", "recall", "f1"], ["sample", "match_type"], False),
        (
            "apa_detect_filter_comb",
            apa_detect_by_filter,
            ["precision", "recall", "f1"],
            ["sample", "match_type"],
            True,
        ),
        (
            "apa_detect_single_filter",
            apa_single,
            ["precision", "recall", "f1"],
            # Use sample only so DaPars2/scMAPA rows align with sidecar tools.
            ["sample"],
            False,
        ),
    ]

    all_rows: list[dict[str, object]] = []

    for dataset, df, metrics, unit_cols_candidates, split_filter_comb in datasets:
        if df.empty:
            continue
        protocols = sorted(df["protocol"].dropna().astype(str).unique()) if "protocol" in df.columns else []
        filter_combs = (
            sorted(df["filter_type_comb"].dropna().astype(str).unique())
            if split_filter_comb and "filter_type_comb" in df.columns
            else []
        )
        for metric in metrics:
            if metric not in df.columns:
                continue
            eps = 0.0 if metric in eps_zero_metrics else float(args.epsilon_default)

            if filter_combs:
                for filter_comb in filter_combs:
                    subset_fc = df[df["filter_type_comb"].astype(str) == filter_comb].copy()
                    if subset_fc.empty:
                        continue
                    all_rows.extend(
                        _panel_win_rate(
                            subset_fc,
                            dataset=dataset,
                            metric=metric,
                            analysis_type="overall_filter_comb",
                            analysis_value=filter_comb,
                            analysis_protocol=None,
                            analysis_filter_comb=filter_comb,
                            eps=eps,
                            tool_order=common_order,
                            balanced_by=str(args.balanced_by),
                            unit_cols_candidates=unit_cols_candidates,
                        )
                    )
                    for protocol in protocols:
                        subset = subset_fc[subset_fc["protocol"].astype(str) == protocol].copy()
                        if subset.empty:
                            continue
                        all_rows.extend(
                            _panel_win_rate(
                                subset,
                                dataset=dataset,
                                metric=metric,
                                analysis_type="protocol_filter_comb",
                                analysis_value=f"{protocol}__{filter_comb}",
                                analysis_protocol=protocol,
                                analysis_filter_comb=filter_comb,
                                eps=eps,
                                tool_order=common_order,
                                balanced_by=str(args.balanced_by),
                                unit_cols_candidates=unit_cols_candidates,
                            )
                        )
                continue

            all_rows.extend(
                _panel_win_rate(
                    df,
                    dataset=dataset,
                    metric=metric,
                    analysis_type="overall",
                    analysis_value="overall",
                    analysis_protocol=None,
                    analysis_filter_comb=None,
                    eps=eps,
                    tool_order=common_order,
                    balanced_by=str(args.balanced_by),
                    unit_cols_candidates=unit_cols_candidates,
                )
            )

            for protocol in protocols:
                subset = df[df["protocol"].astype(str) == protocol].copy()
                if subset.empty:
                    continue
                all_rows.extend(
                    _panel_win_rate(
                        subset,
                        dataset=dataset,
                        metric=metric,
                        analysis_type="protocol",
                        analysis_value=protocol,
                        analysis_protocol=protocol,
                        analysis_filter_comb=None,
                        eps=eps,
                        tool_order=common_order,
                        balanced_by=str(args.balanced_by),
                        unit_cols_candidates=unit_cols_candidates,
                    )
                )

    out_df = pd.DataFrame(all_rows)
    if not out_df.empty:
        out_df = out_df.sort_values(
            [
                "dataset",
                "metric",
                "analysis_type",
                "analysis_value",
                "tool_a",
                "tool_b",
            ]
        ).reset_index(drop=True)

    out_dir = topic_intermediate_dir(TOPIC)
    table_path = out_dir / OUT_TABLE
    meta_path = out_dir / OUT_META

    write_parquet(out_df, table_path)
    write_json(
        {
            "topic": TOPIC,
            "rows": int(len(out_df)),
            "columns": list(out_df.columns),
            "tool_order": common_order,
            "epsilon_default": float(args.epsilon_default),
            "epsilon_zero_metrics": sorted(eps_zero_metrics),
            "balanced_by": str(args.balanced_by),
            "datasets": [name for name, _, _, _, _ in datasets],
            "metrics_by_dataset": {name: metric_list for name, _, metric_list, _, _ in datasets},
        },
        meta_path,
    )

    print(f"Wrote pairwise win-rate table: {table_path}")


if __name__ == "__main__":
    main()
