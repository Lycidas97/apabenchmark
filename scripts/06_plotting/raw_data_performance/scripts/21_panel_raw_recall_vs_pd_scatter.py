#!/usr/bin/env python3
"""Render raw recall vs predicted DE APA count scatter with x/y error bars."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import TOOL_MAP, TOOL_ORDER_CP
from _shared.io import read_table, write_tsv
from _shared.paths import topic_figure_dir, topic_intermediate_dir, topic_root
from _shared.style import apply_reference_style, figsize_from_mm, figsize_from_preset, save_figure, tool_palette

TOPIC = "raw_data_performance"
FILTER_CONFIG_PATH = topic_root(TOPIC) / "config" / "recall_gt_filter.json"
PLOT_CONFIG_PATH = topic_root(TOPIC) / "config" / "plot_panels.json"
DEFAULT_RAW_INPUT = "raw_data_recall_gt_filtered.parquet"
DEFAULT_PD_INPUT = "raw_data_performance_prepared.parquet"
OUTPUT_NAME = "raw_recall_vs_pd_scatter_errorbar.pdf"
POINTS_NAME = "raw_recall_vs_pd_scatter_points.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scatter raw recall vs predicted DE APA count with x/y error bars."
    )
    parser.add_argument("--filter-config", default=str(FILTER_CONFIG_PATH))
    parser.add_argument("--plot-config", default=str(PLOT_CONFIG_PATH))
    return parser.parse_args()


def load_filter_pairs(filter_config_path: Path) -> set[tuple[str, str]]:
    cfg = json.loads(filter_config_path.read_text(encoding="utf-8"))
    pairs = set()
    for item in cfg.get("filter_pairs", []):
        f1 = str(item.get("filter_type_1", "")).strip()
        f2 = str(item.get("filter_type_2", "")).strip()
        if f1 and f2:
            pairs.add((f1, f2))
    if not pairs:
        raise ValueError(f"No filter_pairs configured in: {filter_config_path}")
    return pairs


def _pair_key(f1: str, f2: str) -> str:
    return f"{f1}|{f2}"


def _format_filter_pair_label(pair_key: str) -> str:
    parts = str(pair_key).split("|", 1)
    if len(parts) != 2:
        return str(pair_key)

    left, right = parts
    left_map = {
        "dexseq_0.05": "DexSeq",
        "fisher_0.05": "Fisher",
    }
    left_label = left_map.get(left, left.replace("_", " "))
    if right.startswith("dexseq_log2fc_"):
        right_label = "log2FC > " + right.replace("dexseq_log2fc_", "")
    elif right.startswith("MPRO_"):
        right_label = "MPRO > " + right.replace("MPRO_", "")
    else:
        right_label = right.replace("_", " ")
    return f"{left_label} | {right_label}"


def _resolve_figsize(scatter_cfg: dict, figsize_preset: str) -> tuple[float, float]:
    mm = scatter_cfg.get("figsize_mm")
    if isinstance(mm, list) and len(mm) == 2:
        return figsize_from_mm(float(mm[0]), float(mm[1]))
    return figsize_from_preset(figsize_preset)


def _aggregate_metric_cfg(
    df: pd.DataFrame,
    value_col: str,
    prefix: str,
    groupby_cols: list[str],
    errorbar_stat: str,
) -> pd.DataFrame:
    if errorbar_stat not in {"std", "sem"}:
        raise ValueError(f"Unsupported errorbar_stat={errorbar_stat!r}. Use 'std' or 'sem'.")

    grouped = (
        df.groupby(groupby_cols, as_index=False)[value_col]
        .agg(["mean", "std", "count"])
        .rename(
            columns={
                "mean": f"{prefix}_mean",
                "std": f"{prefix}_sd",
                "count": f"{prefix}_n",
            }
        )
    )
    grouped[f"{prefix}_sd"] = grouped[f"{prefix}_sd"].fillna(0.0)
    if errorbar_stat == "sem":
        n = grouped[f"{prefix}_n"].clip(lower=1)
        grouped[f"{prefix}_sd"] = grouped[f"{prefix}_sd"] / n.pow(0.5)
    return grouped


def _load_recall_with_pd(raw_input: Path, pd_input: Path, pd_metric_domain: str) -> pd.DataFrame:
    raw = read_table(raw_input)
    pd_src = read_table(pd_input)

    required_raw = {"metric_domain", "tool", "sample", "sample_id", "pair_label", "match_type", "filter_type_1", "filter_type_2", "recall", "gt_count", "raw_run_id"}
    required_pd = required_raw | {"pd_count"}
    missing_raw = sorted(required_raw - set(raw.columns))
    missing_pd = sorted(required_pd - set(pd_src.columns))
    if missing_raw:
        raise KeyError(f"Missing required raw recall columns: {missing_raw}")
    if missing_pd:
        raise KeyError(f"Missing required prepared PD columns: {missing_pd}")

    pd_src = pd_src[pd_src["metric_domain"] == pd_metric_domain].copy()
    keys = [
        "metric_domain",
        "tool",
        "protocol",
        "sample",
        "sample_id",
        "pair_label",
        "match_type",
        "filter_type_1",
        "filter_type_2",
        "raw_run_id",
    ]
    keys = [col for col in keys if col in raw.columns and col in pd_src.columns]
    if raw.duplicated(keys).any() or pd_src.duplicated(keys).any():
        raise ValueError(f"Cannot one-to-one join recall and pd_count using keys: {keys}")

    out = raw.merge(pd_src[keys + ["pd_count"]], on=keys, how="left", validate="one_to_one")
    out["recall"] = pd.to_numeric(out["recall"], errors="coerce")
    out["pd_count"] = pd.to_numeric(out["pd_count"], errors="coerce")
    out = out[out["recall"].notna() & out["pd_count"].notna()].copy()
    if out.empty:
        raise ValueError("No rows remained after joining recall with pd_count.")
    return out


def _auto_limits(
    points: pd.DataFrame,
    *,
    x_col: str,
    y_col: str,
    xerr_col: str,
    yerr_col: str,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float | None,
    margin_frac: float,
    x_min_span: float,
    y_min_span: float,
) -> tuple[float, float, float, float]:
    x = pd.to_numeric(points[x_col], errors="coerce")
    y = pd.to_numeric(points[y_col], errors="coerce")
    xerr = pd.to_numeric(points[xerr_col], errors="coerce").fillna(0.0)
    yerr = pd.to_numeric(points[yerr_col], errors="coerce").fillna(0.0)
    valid = x.notna() & y.notna()
    if int(valid.sum()) == 0:
        return x_min, x_max, y_min, y_max if y_max is not None else y_min + y_min_span

    def _bound(lo: float, hi: float, clip_lo: float, clip_hi: float | None, min_span: float) -> tuple[float, float]:
        span = max(hi - lo, min_span)
        pad = span * float(margin_frac)
        a = max(clip_lo, lo - pad)
        b = hi + pad if clip_hi is None else min(clip_hi, hi + pad)
        if b - a < min_span:
            c = 0.5 * (a + b)
            a = max(clip_lo, c - min_span / 2.0)
            b = c + min_span / 2.0 if clip_hi is None else min(clip_hi, c + min_span / 2.0)
        return a, b

    x0, x1 = _bound(float((x[valid] - xerr[valid]).min()), float((x[valid] + xerr[valid]).max()), x_min, x_max, x_min_span)
    y0, y1 = _bound(float((y[valid] - yerr[valid]).min()), float((y[valid] + yerr[valid]).max()), y_min, y_max, y_min_span)
    return x0, x1, y0, y1


def main() -> None:
    args = parse_args()
    apply_reference_style()

    filter_pairs = load_filter_pairs(Path(args.filter_config))
    plot_cfg = json.loads(Path(args.plot_config).read_text(encoding="utf-8"))
    scatter_cfg = plot_cfg.get("recall_vs_pd_scatter", {})
    base_scatter_cfg = plot_cfg.get("scatter", {})

    raw_input = topic_intermediate_dir(TOPIC) / str(scatter_cfg.get("raw_input_name", DEFAULT_RAW_INPUT))
    pd_input = topic_intermediate_dir(TOPIC) / str(scatter_cfg.get("pd_input_name", DEFAULT_PD_INPUT))
    pd_metric_domain = str(scatter_cfg.get("pd_metric_domain", "raw_de_apa_performance"))
    output_figure = str(scatter_cfg.get("output_figure", OUTPUT_NAME))
    output_points = str(scatter_cfg.get("output_points", POINTS_NAME))

    figsize_preset = str(scatter_cfg.get("figsize_preset", base_scatter_cfg.get("figsize_preset", "single_column")))
    x_label = str(scatter_cfg.get("x_label", "Raw Recall (mean ± SD)"))
    y_label = str(scatter_cfg.get("y_label", "Predicted DE APA count (mean ± SD)"))
    marker_size = float(scatter_cfg.get("marker_size", base_scatter_cfg.get("marker_size", 4)))
    capsize = float(scatter_cfg.get("capsize", base_scatter_cfg.get("capsize", 2)))
    show_legend = bool(scatter_cfg.get("show_legend", base_scatter_cfg.get("show_legend", True)))
    legend_loc = str(scatter_cfg.get("legend_loc", base_scatter_cfg.get("legend_loc", "center left")))
    legend_bbox_anchor = scatter_cfg.get("legend_bbox_anchor", base_scatter_cfg.get("legend_bbox_anchor", [0.875, 0.5]))
    legend_fontsize = float(scatter_cfg.get("legend_fontsize", base_scatter_cfg.get("legend_fontsize", 7)))
    legend_title_fontsize = float(scatter_cfg.get("legend_title_fontsize", legend_fontsize))
    legend_title = str(scatter_cfg.get("legend_title", "")).strip()
    legend_ncol = int(scatter_cfg.get("legend_ncol", 1))
    show_correlation = bool(scatter_cfg.get("show_correlation", base_scatter_cfg.get("show_correlation", True)))
    correlation_method = str(scatter_cfg.get("correlation_method", base_scatter_cfg.get("correlation_method", "spearman"))).strip().lower()
    correlation_xy = scatter_cfg.get("correlation_xy", base_scatter_cfg.get("correlation_xy", [0.03, 0.97]))
    correlation_fontsize = float(scatter_cfg.get("correlation_fontsize", base_scatter_cfg.get("correlation_fontsize", 7)))
    axis_label_fontsize = scatter_cfg.get("axis_label_fontsize", base_scatter_cfg.get("axis_label_fontsize", 7))
    tick_label_fontsize = scatter_cfg.get("tick_label_fontsize", base_scatter_cfg.get("tick_label_fontsize", 7))
    panel_title_fontsize = float(scatter_cfg.get("panel_title_fontsize", base_scatter_cfg.get("panel_title_fontsize", 7)))
    panel_layout = scatter_cfg.get("panel_layout", base_scatter_cfg.get("panel_layout", [1, 3]))
    panel_nrows, panel_ncols = int(panel_layout[0]), int(panel_layout[1])
    panel_wspace = float(scatter_cfg.get("panel_wspace", base_scatter_cfg.get("panel_wspace", 0.1)))
    panel_hspace = float(scatter_cfg.get("panel_hspace", base_scatter_cfg.get("panel_hspace", 0.28)))
    figure_right = float(scatter_cfg.get("figure_right", base_scatter_cfg.get("figure_right", 0.86)))
    figure_bottom = float(scatter_cfg.get("figure_bottom", base_scatter_cfg.get("figure_bottom", 0.2)))
    figure_top = float(scatter_cfg.get("figure_top", base_scatter_cfg.get("figure_top", 0.88)))
    errorbar_stat = str(scatter_cfg.get("errorbar_stat", base_scatter_cfg.get("errorbar_stat", "std")))
    groupby_cols = [str(v) for v in scatter_cfg.get("groupby", base_scatter_cfg.get("groupby", ["tool"]))]
    x_min = float(scatter_cfg.get("xmin", 0.0))
    x_max = float(scatter_cfg.get("xmax", 1.02))
    y_min = float(scatter_cfg.get("ymin", 0.0))
    y_max_value = scatter_cfg.get("ymax")
    y_max = None if y_max_value is None else float(y_max_value)
    axis_margin_frac = float(scatter_cfg.get("axis_margin_frac", 0.08))
    x_min_span = float(scatter_cfg.get("x_min_span", scatter_cfg.get("axis_min_span", 0.14)))
    y_min_span = float(scatter_cfg.get("y_min_span", 100.0))
    share_y_axis = bool(scatter_cfg.get("share_y_axis", True))

    data = _load_recall_with_pd(raw_input, pd_input, pd_metric_domain)
    pair_keys = {_pair_key(f1, f2) for f1, f2 in filter_pairs}
    data_key = data["filter_type_1"].astype(str) + "|" + data["filter_type_2"].astype(str)
    data = data[data_key.isin(pair_keys)].copy()
    if data.empty:
        raise ValueError("Recall/PD table has no rows for configured filter pairs.")

    panel_points: list[pd.DataFrame] = []
    for f1, f2 in sorted(filter_pairs):
        sub = data[
            (data["filter_type_1"].astype(str) == str(f1))
            & (data["filter_type_2"].astype(str) == str(f2))
        ].copy()
        if sub.empty:
            continue
        recall_stats = _aggregate_metric_cfg(sub, "recall", "raw_recall", groupby_cols, errorbar_stat)
        pd_stats = _aggregate_metric_cfg(sub, "pd_count", "pd_count", groupby_cols, errorbar_stat)
        recall_stats["tool"] = recall_stats["tool"].map(TOOL_MAP).fillna(recall_stats["tool"])
        pd_stats["tool"] = pd_stats["tool"].map(TOOL_MAP).fillna(pd_stats["tool"])
        merge_cols = [col for col in groupby_cols if col in recall_stats.columns and col in pd_stats.columns]
        points = pd.merge(recall_stats, pd_stats, on=merge_cols, how="inner")
        if points.empty:
            continue
        points["filter_type_1"] = str(f1)
        points["filter_type_2"] = str(f2)
        points["pair_key"] = _pair_key(str(f1), str(f2))
        panel_points.append(points)

    if not panel_points:
        raise ValueError("No overlapping recall/PD groups after aggregation.")

    points_all = pd.concat(panel_points, ignore_index=True, sort=False)
    present_tools = set(points_all["tool"].astype(str))
    legend_order = [t for t in TOOL_ORDER_CP if t in present_tools]
    legend_order.extend(sorted(t for t in present_tools if t not in legend_order))
    points_all["tool"] = pd.Categorical(points_all["tool"], categories=legend_order, ordered=True)
    points_all = points_all.sort_values(["pair_key", "tool"]).reset_index(drop=True)
    write_tsv(points_all, topic_intermediate_dir(TOPIC) / output_points)

    panel_keys = points_all["pair_key"].dropna().astype(str).unique().tolist()
    nrows, ncols = panel_nrows, panel_ncols
    while nrows * ncols < len(panel_keys):
        ncols += 1

    global_limits = _auto_limits(
        points_all,
        x_col="raw_recall_mean",
        y_col="pd_count_mean",
        xerr_col="raw_recall_sd",
        yerr_col="pd_count_sd",
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
        margin_frac=axis_margin_frac,
        x_min_span=x_min_span,
        y_min_span=y_min_span,
    )

    fig, axes = plt.subplots(nrows, ncols, figsize=_resolve_figsize(scatter_cfg, figsize_preset), sharey=share_y_axis)
    fig.subplots_adjust(wspace=panel_wspace, hspace=panel_hspace, right=figure_right, bottom=figure_bottom, top=figure_top)
    axes_flat = axes.reshape(-1) if isinstance(axes, np.ndarray) else np.array([axes], dtype=object)
    colors = tool_palette(legend_order)

    for idx, ax in enumerate(axes_flat):
        if idx >= len(panel_keys):
            ax.axis("off")
            continue
        pair_key = panel_keys[idx]
        panel_df = points_all[points_all["pair_key"].astype(str) == pair_key].copy().sort_values("tool")
        for _, row in panel_df.iterrows():
            tool = str(row["tool"])
            c = colors.get(tool, "#929292")
            ax.errorbar(
                float(row["raw_recall_mean"]),
                float(row["pd_count_mean"]),
                xerr=float(row["raw_recall_sd"]),
                yerr=float(row["pd_count_sd"]),
                fmt="o",
                color=c,
                ecolor=c,
                elinewidth=0.8,
                capsize=capsize,
                markersize=marker_size,
                alpha=0.95,
            )

        px0, px1, _, _ = _auto_limits(
            panel_df,
            x_col="raw_recall_mean",
            y_col="pd_count_mean",
            xerr_col="raw_recall_sd",
            yerr_col="pd_count_sd",
            x_min=x_min,
            x_max=x_max,
            y_min=y_min,
            y_max=y_max,
            margin_frac=axis_margin_frac,
            x_min_span=x_min_span,
            y_min_span=y_min_span,
        )
        ax.set_xlim(px0, px1)
        ax.set_ylim(global_limits[2], global_limits[3])
        ax.set_xlabel(x_label if idx == (len(panel_keys) // 2) else "", fontsize=axis_label_fontsize)
        if idx % ncols == 0:
            ax.set_ylabel(y_label, fontsize=axis_label_fontsize)
        else:
            ax.set_ylabel("")
            ax.tick_params(axis="y", labelleft=False)
        if tick_label_fontsize is not None:
            ax.tick_params(axis="both", labelsize=float(tick_label_fontsize))
        ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.7)
        ax.set_title(_format_filter_pair_label(pair_key), fontsize=panel_title_fontsize)

        if show_correlation and len(panel_df) >= 2:
            corr_method = correlation_method if correlation_method in {"pearson", "spearman"} else "spearman"
            x_ser = pd.to_numeric(panel_df["raw_recall_mean"], errors="coerce")
            y_ser = pd.to_numeric(panel_df["pd_count_mean"], errors="coerce")
            valid = x_ser.notna() & y_ser.notna()
            if int(valid.sum()) >= 2:
                r = x_ser[valid].corr(y_ser[valid], method=corr_method)
                if not np.isnan(r):
                    corr_xy = correlation_xy if isinstance(correlation_xy, list) and len(correlation_xy) == 2 else [0.03, 0.97]
                    ax.text(
                        float(corr_xy[0]),
                        float(corr_xy[1]),
                        f"{corr_method.title()} r = {r:.3f}",
                        transform=ax.transAxes,
                        fontsize=correlation_fontsize,
                        va="top",
                        ha="left",
                        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.7, "pad": 1.5},
                    )

    if show_legend:
        legend_bbox = legend_bbox_anchor if isinstance(legend_bbox_anchor, list) and len(legend_bbox_anchor) == 2 else [1.02, 0.5]
        handles = [
            Line2D([0], [0], marker="o", linestyle="", color=colors.get(tool, "#929292"), label=tool, markersize=4)
            for tool in legend_order
        ]
        fig.legend(
            handles=handles,
            labels=legend_order,
            loc=legend_loc,
            bbox_to_anchor=(float(legend_bbox[0]), float(legend_bbox[1])),
            frameon=False,
            fontsize=legend_fontsize,
            title=(legend_title if legend_title else None),
            title_fontsize=legend_title_fontsize,
            ncol=max(1, legend_ncol),
        )

    out_fig = topic_figure_dir(TOPIC) / output_figure
    save_figure(fig, out_fig)
    plt.close(fig)
    print(f"Wrote scatter data: {topic_intermediate_dir(TOPIC) / output_points}")
    print(f"Wrote figure: {out_fig}")


if __name__ == "__main__":
    main()
