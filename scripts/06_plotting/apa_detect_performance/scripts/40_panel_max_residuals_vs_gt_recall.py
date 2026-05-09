#!/usr/bin/env python3
"""Render max residual vs recall panel."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd
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
OUTPUT_NAME = "max_residuals_vs_gt_recall_with_arrows.pdf"


def pick_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for col in candidates:
        if col in df.columns:
            return col
    return None


def main() -> None:
    apply_reference_style()
    df = read_table(topic_intermediate_dir(TOPIC) / INPUT_NAME)
    out_path = topic_figure_dir(TOPIC) / OUTPUT_NAME

    fig, ax = plt.subplots(figsize=figsize_from_preset("residual_scatter"))

    x_col = pick_column(df, ["gt_recall", "recall"])
    y_col = pick_column(df, ["max_residual", "residual", "residuals"])
    hue_col = "match_type" if "match_type" in df.columns else None

    if x_col and y_col:
        sns.scatterplot(data=df, x=x_col, y=y_col, hue=hue_col, alpha=0.7, s=10, ax=ax)
        q1 = df[y_col].quantile(1 / 3)
        q2 = df[y_col].quantile(2 / 3)
        ax.axhline(y=q1, color="black", linestyle="--", linewidth=0.8)
        ax.axhline(y=q2, color="black", linestyle="--", linewidth=0.8)
        ax.grid(True, linestyle=":", linewidth=0.3)
        ax.set_title("Max Residual vs Recall", loc="left", fontsize=7)
    else:
        draw_placeholder(
            ax,
            "Max Residual vs Recall",
            "Required columns not found in intermediate table.\n"
            "Need one recall column and one residual column.",
        )

    save_figure(fig, out_path)
    plt.close(fig)
    print(f"Wrote figure: {out_path}")


if __name__ == "__main__":
    main()
