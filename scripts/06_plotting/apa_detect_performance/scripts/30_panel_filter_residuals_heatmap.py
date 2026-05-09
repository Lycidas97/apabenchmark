#!/usr/bin/env python3
"""Render APA filter residual heatmap panel."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[2]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.io import read_table
from _shared.paths import topic_figure_dir, topic_intermediate_dir
from _shared.style import apply_reference_style, draw_placeholder, figsize_from_preset, save_figure

TOPIC = "apa_detect_performance"
INPUT_NAME = "apa_detect_performance_prepared.tsv"
OUTPUT_NAME = "filter_residuals_heatmap.pdf"


def pick_residual_column(columns: list[str]) -> str | None:
    candidates = ["max_residual", "residual", "residuals"]
    for col in candidates:
        if col in columns:
            return col
    return None


def main() -> None:
    apply_reference_style()
    df = read_table(topic_intermediate_dir(TOPIC) / INPUT_NAME)
    out_path = topic_figure_dir(TOPIC) / OUTPUT_NAME

    fig, ax = plt.subplots(figsize=figsize_from_preset("single_column"))

    residual_col = pick_residual_column(list(df.columns))
    if residual_col and {"filter_type_1", "filter_type_2"}.issubset(df.columns):
        pivot = df.pivot_table(
            index="filter_type_1",
            columns="filter_type_2",
            values=residual_col,
            aggfunc="mean",
        )
        sns.heatmap(pivot, cmap="Blues", ax=ax, cbar_kws={"label": residual_col})
        ax.set_xlabel("filter_type_2")
        ax.set_ylabel("filter_type_1")
        ax.set_title("Filter Residual Heatmap", loc="left", fontsize=7)
    else:
        draw_placeholder(
            ax,
            "Filter Residual Heatmap",
            "Required columns not found in intermediate table.\n"
            "Need: filter_type_1, filter_type_2, and one residual column.",
        )

    save_figure(fig, out_path)
    plt.close(fig)
    print(f"Wrote figure: {out_path}")


if __name__ == "__main__":
    main()
