#!/usr/bin/env python3
"""Render max residual vs ground-truth recall scatter panel."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import PALETTE_HEX
from _shared.paths import topic_figure_dir
from _shared.style import apply_reference_style, save_figure

from _plot_helpers import TOPIC, format_filter_type_comb, load_apa_residual_by_filter_table
from _plot_config import apply_runtime_rcparams, figsize_for_figure

OUTPUT_NAME = "max_residuals_vs_gt_recall_with_arrows.pdf"


def main() -> None:
    apply_reference_style()
    apply_runtime_rcparams()
    df = load_apa_residual_by_filter_table()

    required = {"filter_type_comb", "f1_residuals", "gt_recall"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise KeyError(f"Missing required columns for scatter panel: {missing}")
    if df.empty:
        raise ValueError("Residual-by-filter table is empty.")

    fig, ax = plt.subplots(figsize=figsize_for_figure(OUTPUT_NAME, "residual_scatter_main"))

    sns.scatterplot(
        data=df,
        x="gt_recall",
        y="f1_residuals",
        color="grey",
        alpha=0.5,
        label="_nolegend_",
        s=14,
        ax=ax,
    )

    y_min = float(df["gt_recall"].min())
    y_max = float(df["gt_recall"].max())
    third_1 = y_min + (y_max - y_min) / 3.0
    third_2 = y_min + 2.0 * (y_max - y_min) / 3.0
    ax.axvline(x=third_1, color="black", linestyle="--", linewidth=0.8)
    ax.axvline(x=third_2, color="black", linestyle="--", linewidth=0.8)

    sections = [
        df[df["gt_recall"] <= third_1],
        df[(df["gt_recall"] > third_1) & (df["gt_recall"] <= third_2)],
        df[df["gt_recall"] > third_2],
    ]

    point_colors = [PALETTE_HEX[0], PALETTE_HEX[1], PALETTE_HEX[2], PALETTE_HEX[3], PALETTE_HEX[4]]
    x_span = max(1e-9, float(df["gt_recall"].max() - df["gt_recall"].min()))
    y_span = max(1e-9, float(df["f1_residuals"].max() - df["f1_residuals"].min()))
    color_index = 0
    for section in sections:
        if section.empty:
            continue
        # Select the true top residual in each section.
        point = section.nlargest(1, "f1_residuals").iloc[0]
        x_coord = float(point["gt_recall"])
        y_coord = float(point["f1_residuals"])
        label = format_filter_type_comb(str(point["filter_type_comb"]))
        color = point_colors[color_index % len(point_colors)]
        color_index += 1

        ax.scatter(x_coord, y_coord, color=color, label=label, s=16, zorder=5)
        ax.annotate(
            text="",
            xy=(x_coord, y_coord),
            xytext=(x_coord - 0.03 * x_span, y_coord + 0.04 * y_span),
            arrowprops={"arrowstyle": "->", "color": color, "linewidth": 1.2},
        )

    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.set_xlabel("Ground Truth Recall")
    ax.set_ylabel("Mean Residual")
    ax.legend(
        bbox_to_anchor=(0.5, -0.25),
        loc="lower center",
        borderaxespad=0.0,
        fontsize=8,
        ncols=1,
        frameon=False,
    )

    out_path = topic_figure_dir(TOPIC) / OUTPUT_NAME
    save_figure(fig, out_path)
    plt.close(fig)
    print(f"Wrote figure: {out_path}")


if __name__ == "__main__":
    main()
