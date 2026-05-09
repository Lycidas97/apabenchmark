#!/usr/bin/env python3
"""Render AIC/BIC comparison panel for peak model params."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.io import read_table
from _shared.paths import topic_figure_dir, topic_intermediate_dir
from _shared.style import apply_reference_style, draw_placeholder, figsize_from_preset, palette, save_figure

TOPIC = "peak_model_params"
INPUT_NAME = "peak_model_params_prepared.tsv"
OUTPUT_NAME = "AIC_BIC_boxplot.pdf"
FLIER_PROPS = {
    "markersize": 0.6,
    "markeredgewidth": 0.15,
    "markerfacecolor": "#666666",
    "markeredgecolor": "#666666",
    "alpha": 0.5,
}


def main() -> None:
    apply_reference_style()
    df = read_table(topic_intermediate_dir(TOPIC) / INPUT_NAME)

    fig, ax = plt.subplots(figsize=figsize_from_preset("small_square"))

    if {"delta_AIC", "delta_BIC"}.issubset(df.columns):
        sns.boxplot(
            data=df[["delta_AIC", "delta_BIC"]],
            palette=palette(2),
            flierprops=FLIER_PROPS,
            ax=ax,
        )
        ax.axhline(0, color="gray", linestyle="--")
        ax.set_ylabel("Delta score")
        ax.set_xlabel("")
        ax.set_title("AIC/BIC Comparison", loc="left", fontsize=7)
    else:
        draw_placeholder(
            ax,
            "AIC/BIC Comparison",
            "Required columns not found in intermediate table.\n"
            "Need: delta_AIC and delta_BIC.",
        )

    out_path = topic_figure_dir(TOPIC) / OUTPUT_NAME
    save_figure(fig, out_path)
    plt.close(fig)
    print(f"Wrote figure: {out_path}")


if __name__ == "__main__":
    main()
