#!/usr/bin/env python3
"""Render PAS detect metric-by-protocol panels."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[2]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import PROTOCOL_ORDER, TOOL_ORDER
from _shared.io import read_table
from _shared.paths import topic_figure_dir, topic_intermediate_dir
from _shared.style import apply_reference_style, figsize_from_preset, palette, save_figure

TOPIC = "pas_detect_performance"
INPUT_NAME = "pas_detect_performance_prepared.tsv"
METRICS = ["precision", "recall", "f1"]


def render_metric(df, metric: str) -> None:
    fig, axes = plt.subplots(
        nrows=len(PROTOCOL_ORDER),
        ncols=1,
        figsize=figsize_from_preset("narrow_tall"),
        dpi=300,
    )

    if len(PROTOCOL_ORDER) == 1:
        axes = [axes]

    for idx, (protocol, ax) in enumerate(zip(PROTOCOL_ORDER, axes)):
        plot_data = df[(df["protocol"] == protocol) & (df[metric].notna())]
        sns.boxplot(
            data=plot_data,
            x="tool",
            y=metric,
            hue="match_type",
            order=TOOL_ORDER,
            palette=palette(),
            ax=ax,
            dodge=True,
            width=0.7,
            flierprops={
                "marker": "o",
                "markersize": 0.5,
                "markerfacecolor": "gray",
                "markeredgecolor": "gray",
            },
            linewidth=0.4,
        )
        ax.set_title(protocol, loc="left", pad=2, fontsize=6)
        ax.set_ylabel(metric.capitalize() if idx == 0 else "", labelpad=2)
        ax.set_xlabel("")
        ax.tick_params(axis="both", labelsize=6, pad=1)
        ax.grid(True, linestyle=":", linewidth=0.3)

        if idx != len(PROTOCOL_ORDER) - 1:
            ax.set_xticklabels([])
        else:
            plt.setp(ax.get_xticklabels(), rotation=90)

        leg = ax.get_legend()
        if leg:
            leg.remove()

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.93),
        ncol=max(1, len(labels)),
        title="Match Type",
        frameon=False,
        title_fontsize=7,
        fontsize=6,
        columnspacing=0.8,
        handletextpad=0.3,
    )

    fig.subplots_adjust(hspace=0.25, left=0.12, right=0.88, bottom=0.08)
    out_path = topic_figure_dir(TOPIC) / f"pas_detect_{metric}_by_protocol.pdf"
    save_figure(fig, out_path)
    plt.close(fig)
    print(f"Wrote figure: {out_path}")


def main() -> None:
    apply_reference_style()
    df = read_table(topic_intermediate_dir(TOPIC) / INPUT_NAME)
    for metric in METRICS:
        render_metric(df, metric)


if __name__ == "__main__":
    main()
