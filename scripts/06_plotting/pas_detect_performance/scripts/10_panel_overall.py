#!/usr/bin/env python3
"""Render overall PAS detect performance panel."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[2]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import TOOL_ORDER
from _shared.io import read_table
from _shared.paths import topic_figure_dir, topic_intermediate_dir
from _shared.style import apply_reference_style, figsize_from_preset, palette, save_figure

TOPIC = "pas_detect_performance"
INPUT_NAME = "pas_detect_performance_prepared.tsv"
OUTPUT_NAME = "pas_detect_performance.pdf"
METRICS = ["precision", "recall", "f1"]


def main() -> None:
    apply_reference_style()
    df_path = topic_intermediate_dir(TOPIC) / INPUT_NAME
    df = read_table(df_path)

    fig, axes = plt.subplots(3, 1, figsize=figsize_from_preset("single_column"))
    fig.subplots_adjust(hspace=0.25)

    outlier_props = {
        "flierprops": {
            "marker": "o",
            "markersize": 0.1,
            "markerfacecolor": "gray",
            "markeredgecolor": "gray",
        }
    }

    for idx, metric in enumerate(METRICS):
        ax = axes[idx]
        sns.boxplot(
            data=df,
            x="tool",
            y=metric,
            hue="match_type",
            order=TOOL_ORDER,
            palette=palette(),
            dodge=True,
            width=0.8,
            ax=ax,
            **outlier_props,
        )
        ax.grid(linestyle="--", alpha=0.6, linewidth=0.5)
        ax.set_ylabel(metric.capitalize())
        if idx == 0:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(
                handles,
                labels,
                title="Match Type",
                title_fontsize=6,
                fontsize=6,
                frameon=False,
                bbox_to_anchor=(1.02, 0.5),
                loc="center left",
                borderaxespad=0.0,
            )
        else:
            leg = ax.get_legend()
            if leg:
                leg.remove()
        if idx != len(METRICS) - 1:
            ax.set_xlabel("")
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Tool")

    out_path = topic_figure_dir(TOPIC) / OUTPUT_NAME
    save_figure(fig, out_path)
    plt.close(fig)
    print(f"Wrote figure: {out_path}")


if __name__ == "__main__":
    main()
