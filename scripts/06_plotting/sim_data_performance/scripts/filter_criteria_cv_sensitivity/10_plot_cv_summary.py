#!/usr/bin/env python3
"""Plot DE-APA criterion CV sensitivity summary panels."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[2]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.style import apply_reference_style  # noqa: E402
from _shared.constants import TOOL_ORDER  # noqa: E402

OUTPUT_DIR = SCRIPT_DIR / "output"
FIG_DIR = SCRIPT_DIR / "figures"
TOPIC_DIR = SCRIPT_DIR.parents[1]
FILTER_GRID_TABLE = TOPIC_DIR / "data" / "intermediate" / "sim_data_apa_residuals_by_filter.parquet"
TIER_ORDER = ["strict", "moderate", "lenient"]
TIER_LABELS = {"strict": "S", "moderate": "M", "lenient": "L"}
SPLIT_TYPE_LABELS = {
    "random_half_split": "Random",
    "leave_protocol_out": "LPO",
}
FILTER_TYPE_1_MAP = {
    "dexseq_0.05": "DEXSeq",
    "fisher_0.05": "Fisher",
    "wilcox_PDUI_0.05": "Wilcoxon-PDUI",
    "wilcox_PPUI_0.05": "Wilcoxon-PPUI",
    "wilcox_RWUI_0.05": "Wilcoxon-RWUI",
    "wilcox_DWUI_0.05": "Wilcoxon-DWUI",
}
FILTER_TYPE_1_ORDER = [
    "dexseq_0.05",
    "fisher_0.05",
    "wilcox_PDUI_0.05",
    "wilcox_PPUI_0.05",
    "wilcox_RWUI_0.05",
    "wilcox_DWUI_0.05",
]
FILTER_TYPE_2_ORDER = [
    "PDUI_0.1",
    "PDUI_0.2",
    "PDUI_0.3",
    "PDUI_0.4",
    "PDUI_0.5",
    "PPUI_0.1",
    "PPUI_0.2",
    "PPUI_0.3",
    "PPUI_0.4",
    "PPUI_0.5",
    "RWUI_0.1",
    "RWUI_0.2",
    "RWUI_0.3",
    "RWUI_0.4",
    "RWUI_0.5",
    "DWUI_0.1",
    "DWUI_0.2",
    "DWUI_0.3",
    "DWUI_0.4",
    "DWUI_0.5",
    "MPRO_0.1",
    "MPRO_0.2",
    "MPRO_0.3",
    "MPRO_0.4",
    "MPRO_0.5",
    "dexseq_log2fc_0.25",
    "dexseq_log2fc_0.5",
    "dexseq_log2fc_0.75",
    "dexseq_log2fc_1",
    "dexseq_log2fc_1.25",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=FIG_DIR)
    parser.add_argument("--output-name", default="filter_criteria_cv_sensitivity.pdf")
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing required table: {path}")
    return pd.read_csv(path, sep="\t")


def read_reference_criteria(input_dir: Path) -> pd.DataFrame:
    path = input_dir / "cv_metadata.json"
    if not path.exists():
        return pd.DataFrame(columns=["filter_type_1", "filter_type_2", "tier"])
    metadata = json.loads(path.read_text(encoding="utf-8"))
    return pd.DataFrame(metadata.get("reference_criteria", []))


def filter_family(name: str) -> str:
    if name.startswith("dexseq_log2fc"):
        return "log2FC"
    if name.startswith("MPRO"):
        return "MPRO"
    return f"\u0394{name.split('_', 1)[0]}"


def filter_threshold(name: str) -> str:
    raw = name.rsplit("_", 1)[-1]
    if raw == "1":
        return "1.0"
    return raw


def read_filter_grid(input_dir: Path, selected: pd.DataFrame) -> pd.DataFrame:
    columns = ["filter_type_1", "filter_type_2", "filter_type_comb"]
    if FILTER_GRID_TABLE.exists():
        grid = pd.read_parquet(FILTER_GRID_TABLE, columns=columns).drop_duplicates()
    elif (input_dir / "cv_filter_grid.tsv").exists():
        grid = read_tsv(input_dir / "cv_filter_grid.tsv")
    else:
        grid = selected[columns].drop_duplicates()
    return grid.dropna(subset=columns).copy()


def panel_design(ax: plt.Axes) -> None:
    ax.axis("off")
    labels = ["Simulated data", "Grouped split", "Train: select", "Held-out: rank", "Stability"]
    xs = [0.11, 0.31, 0.51, 0.71, 0.89]
    y = 0.55
    box_w = 0.14
    box_h = 0.34
    for label, x in zip(labels, xs):
        ax.add_patch(
            FancyBboxPatch(
                (x - box_w / 2, y - box_h / 2),
                box_w,
                box_h,
                boxstyle="round,pad=0.01",
                facecolor="white",
                edgecolor="#333333",
                linewidth=0.8,
                transform=ax.transAxes,
                clip_on=False,
            )
        )
        ax.text(
            x,
            y,
            label,
            ha="center",
            va="center",
            fontsize=8,
            transform=ax.transAxes,
        )
    for start_x, end_x in [(0.19, 0.23), (0.39, 0.43), (0.59, 0.63), (0.79, 0.82)]:
        start = (start_x, y)
        end = (end_x, y)
        ax.annotate("", xy=end, xytext=start, xycoords="axes fraction", arrowprops={"arrowstyle": "->", "lw": 0.8})


def panel_selection_frequency(
    ax: plt.Axes,
    selected: pd.DataFrame,
    filter_grid: pd.DataFrame,
    reference_criteria: pd.DataFrame,
) -> None:
    df = selected[selected["split_type"] == "random_half_split"].copy()
    if df.empty:
        ax.text(0.5, 0.5, "No selection table", ha="center", va="center")
        ax.axis("off")
        return
    n_splits = max(1, df["split_id"].nunique())
    freq = (df.groupby("filter_type_comb").size() / n_splits).reset_index(name="frequency")
    shared_grid = filter_grid[
        ~filter_grid["filter_type_1"].astype(str).str.contains("scmapa|dapars", case=False, na=False)
        & ~filter_grid["filter_type_2"].astype(str).str.contains("scmapa|dapars", case=False, na=False)
        & ~filter_grid["filter_type_comb"].astype(str).str.contains("scmapa|dapars", case=False, na=False)
    ].copy()
    plot_df = shared_grid.merge(freq, on="filter_type_comb", how="left")
    plot_df["frequency"] = plot_df["frequency"].fillna(0.0)

    observed_f1 = set(plot_df["filter_type_1"].astype(str))
    observed_f2 = set(plot_df["filter_type_2"].astype(str))
    f1_order = [name for name in FILTER_TYPE_1_ORDER if name in observed_f1]
    f1_order.extend(sorted(observed_f1.difference(f1_order)))
    f2_order = [name for name in FILTER_TYPE_2_ORDER if name in observed_f2]
    f2_order.extend(sorted(observed_f2.difference(f2_order)))

    pivot = (
        plot_df.pivot_table(
            index="filter_type_1",
            columns="filter_type_2",
            values="frequency",
            aggfunc="max",
            fill_value=0.0,
            observed=False,
        )
        .reindex(index=f1_order, columns=f2_order)
        .fillna(0.0)
    )
    sns.heatmap(
        pivot,
        cmap="Blues",
        vmin=0,
        vmax=1,
        linewidths=0.35,
        linecolor="#d9d9d9",
        cbar_kws={"label": "Selection frequency", "fraction": 0.035, "pad": 0.015},
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Statistical test")
    ax.set_xticks(np.arange(len(pivot.columns)) + 0.5)
    ax.set_yticks(np.arange(len(pivot.index)) + 0.5)
    ax.set_xticklabels([filter_threshold(str(x)) for x in pivot.columns], rotation=90, ha="center", fontsize=8)
    ax.set_yticklabels([FILTER_TYPE_1_MAP.get(str(x), str(x)) for x in pivot.index], rotation=0, fontsize=8)
    ax.tick_params(axis="x", pad=1)
    ax.tick_params(axis="y", pad=2)
    families = [filter_family(str(col)) for col in pivot.columns]
    group_start = 0
    for idx in range(1, len(families) + 1):
        if idx == len(families) or families[idx] != families[group_start]:
            center = (group_start + idx) / 2
            family_width = idx - group_start
            ax.add_patch(
                Rectangle(
                    (group_start, -0.27),
                    family_width,
                    0.12,
                    fill=False,
                    edgecolor="#bdbdbd",
                    linewidth=0.45,
                    transform=ax.get_xaxis_transform(),
                    clip_on=False,
                )
            )
            ax.text(
                center,
                -0.21,
                families[group_start],
                ha="center",
                va="center",
                fontsize=8,
                transform=ax.get_xaxis_transform(),
                clip_on=False,
            )
            if idx < len(families):
                ax.axvspan(idx - 0.045, idx + 0.045, color="white", zorder=4, lw=0)
                ax.axvline(idx - 0.045, color="#c7c7c7", lw=0.5, zorder=5)
                ax.axvline(idx + 0.045, color="#c7c7c7", lw=0.5, zorder=5)
            group_start = idx

    row_index = {str(name): i for i, name in enumerate(pivot.index)}
    col_index = {str(name): i for i, name in enumerate(pivot.columns)}
    for _, row in reference_criteria.iterrows():
        r = row_index.get(str(row.get("filter_type_1")))
        c = col_index.get(str(row.get("filter_type_2")))
        label = TIER_LABELS.get(str(row.get("tier")), "")
        if r is None or c is None or not label:
            continue
        value = float(pivot.iloc[r, c])
        text_color = "white" if value >= 0.5 else "black"
        ax.add_patch(Rectangle((c, r), 1, 1, fill=False, edgecolor="black", linewidth=1.0))
        ax.text(c + 0.5, r + 0.5, label, ha="center", va="center", fontsize=9, fontweight="bold", color=text_color)


def panel_criterion_concordance(ax: plt.Axes, criterion_concordance: pd.DataFrame):
    df = criterion_concordance.copy()
    if df.empty:
        ax.text(0.5, 0.5, "No criterion concordance table", ha="center", va="center")
        ax.axis("off")
        return [], []
    value_col = "spearman_rho"
    if value_col not in df.columns:
        raise KeyError("Missing spearman_rho in cv_criterion_concordance.tsv. Rerun 00_prepare_cv_tables.py.")
    df["split_type_label"] = df["split_type"].map(SPLIT_TYPE_LABELS).fillna(df["split_type"])
    sns.boxplot(
        data=df,
        x="tier",
        y=value_col,
        hue="split_type_label",
        order=TIER_ORDER,
        fliersize=0,
        linewidth=1.1,
        ax=ax,
    )
    sns.stripplot(
        data=df,
        x="tier",
        y=value_col,
        hue="split_type_label",
        order=TIER_ORDER,
        dodge=True,
        jitter=0.025,
        size=1.4,
        alpha=0.35,
        palette={"Random": "#555555", "LPO": "#333333"},
        legend=False,
        ax=ax,
    )
    ax.axhline(0, color="#333333", lw=0.6)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("")
    ax.set_ylabel("Criterion concordance")
    ax.tick_params(axis="x", labelrotation=30, labelsize=7)
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.get_legend()
    if legend is not None:
        legend.remove()
    return handles, labels


def panel_rank_frequency(ax: plt.Axes, metrics: pd.DataFrame) -> None:
    df = metrics[(metrics["split_type"] == "random_half_split") & (metrics["tier"] == "moderate")].copy()
    if df.empty:
        ax.text(0.5, 0.5, "No rank table", ha="center", va="center")
        ax.axis("off")
        return
    freq = (
        df.groupby(["tool", "rank"]).size()
        / df.groupby("tool").size()
    ).reset_index(name="frequency")
    pivot = freq.pivot(index="tool", columns="rank", values="frequency").fillna(0)
    observed_tools = set(pivot.index.astype(str))
    tool_order = [tool for tool in TOOL_ORDER if tool in observed_tools]
    tool_order.extend(sorted(observed_tools.difference(tool_order)))
    pivot = pivot.reindex(index=tool_order)
    sns.heatmap(
        pivot,
        cmap="Reds",
        vmin=0,
        vmax=1,
        linewidths=0.35,
        linecolor="#d9d9d9",
        cbar_kws={"label": "Frequency", "fraction": 0.06, "pad": 0.02},
        ax=ax,
    )
    ax.set_xlabel("Held-out F1 rank")
    ax.set_ylabel("")
    ax.tick_params(axis="both", labelsize=7)


def panel_concordance(ax: plt.Axes, concordance: pd.DataFrame):
    df = concordance.melt(
        id_vars=["split_id", "split_type", "tier"],
        value_vars=["spearman_rho", "kendall_tau"],
        var_name="metric",
        value_name="value",
    ).dropna()
    if df.empty:
        ax.text(0.5, 0.5, "No concordance table", ha="center", va="center")
        ax.axis("off")
        return [], []
    df["split_type_label"] = df["split_type"].map(SPLIT_TYPE_LABELS).fillna(df["split_type"])
    df["metric_label"] = df["metric"].replace({"spearman_rho": "Spearman", "kendall_tau": "Kendall"})
    sns.boxplot(data=df, x="split_type_label", y="value", hue="metric_label", fliersize=0, linewidth=1.1, ax=ax)
    sns.stripplot(
        data=df[df["split_type"] == "leave_protocol_out"],
        x="split_type_label",
        y="value",
        hue="metric_label",
        dodge=True,
        jitter=0.035,
        size=2.5,
        alpha=0.8,
        palette={"Spearman": "#333333", "Kendall": "#333333"},
        legend=False,
        ax=ax,
    )
    ax.set_ylim(0, 1.02)
    ax.axhline(0, color="#333333", lw=0.6)
    ax.set_xlabel("")
    ax.set_ylabel("Rank concordance")
    ax.tick_params(axis="x", labelsize=7)
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.get_legend()
    if legend is not None:
        legend.remove()
    return handles, labels


def main() -> None:
    args = parse_args()
    apply_reference_style()
    selected = read_tsv(args.input_dir / "cv_selected_criteria.tsv")
    metrics = read_tsv(args.input_dir / "cv_tool_metrics.tsv")
    criterion_concordance = read_tsv(args.input_dir / "cv_criterion_concordance.tsv")
    concordance = read_tsv(args.input_dir / "cv_rank_concordance.tsv")
    filter_grid = read_filter_grid(args.input_dir, selected)
    reference_criteria = read_reference_criteria(args.input_dir)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(178 / 25.4, 202 / 25.4))
    outer = GridSpec(
        4,
        1,
        figure=fig,
        height_ratios=[0.25, 1.52, 0.30, 1.12],
        hspace=0.16,
    )
    panel_design(fig.add_subplot(outer[0, 0]))
    panel_selection_frequency(fig.add_subplot(outer[1, 0]), selected, filter_grid, reference_criteria)
    spacer = fig.add_subplot(outer[2, 0])
    spacer.axis("off")

    bottom = outer[3, 0].subgridspec(
        2,
        4,
        height_ratios=[0.18, 1.0],
        width_ratios=[0.88, 0.88, 1.35, 1.35],
        hspace=0.03,
        wspace=0.58,
    )
    criterion_legend_ax = fig.add_subplot(bottom[0, 0])
    concord_legend_ax = fig.add_subplot(bottom[0, 1])
    for legend_ax in (criterion_legend_ax, concord_legend_ax):
        legend_ax.axis("off")

    criterion_handles, criterion_labels = panel_criterion_concordance(
        fig.add_subplot(bottom[1, 0]), criterion_concordance
    )
    concord_handles, concord_labels = panel_concordance(fig.add_subplot(bottom[1, 1]), concordance)
    panel_rank_frequency(fig.add_subplot(bottom[:, 2:4]), metrics)
    if criterion_handles:
        criterion_legend_ax.legend(
            criterion_handles, criterion_labels, title="", fontsize=7, frameon=False, loc="center left"
        )
    if concord_handles:
        concord_legend_ax.legend(concord_handles, concord_labels, title="", fontsize=7, frameon=False, loc="center left")

    out_path = args.output_dir / args.output_name
    fig.savefig(out_path, bbox_inches="tight", dpi=plt.rcParams.get("savefig.dpi", 300))
    plt.close(fig)
    print(f"Wrote figure: {out_path}")


if __name__ == "__main__":
    main()
