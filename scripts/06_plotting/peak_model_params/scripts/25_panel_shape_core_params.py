#!/usr/bin/env python3
"""Render protocol-wise boxplots for selected shape core parameters."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import pandas as pd
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import PROTOCOL_ORDER
from _shared.io import read_table
from _shared.paths import topic_figure_dir, topic_intermediate_dir
from _shared.style import (
    apply_reference_style,
    draw_placeholder,
    figsize_from_mm,
    figsize_from_preset,
    palette,
    save_figure,
)

TOPIC = "peak_model_params"
INPUT_NAME = "peak_model_params_prepared.tsv"
OUTPUT_NAME = "shape_core_params.pdf"
PARAMETERS = [
    "shape_per_peak_upstream_asymmetry_log2_mean",
    "shape_per_peak_effective_width_q10_q90_mean",
    "shape_per_peak_average_distance_mean",
]
DEFAULT_FIG_PRESET = "mid_wide"
FLIER_PROPS = {
    "markersize": 0.6,
    "markeredgewidth": 0.15,
    "markerfacecolor": "#666666",
    "markeredgecolor": "#666666",
    "alpha": 0.5,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render protocol-wise boxplots for selected shape core parameters"
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=topic_intermediate_dir(TOPIC) / INPUT_NAME,
        help="Prepared input table",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=topic_figure_dir(TOPIC) / OUTPUT_NAME,
        help="Output figure path",
    )
    parser.add_argument(
        "--fig-preset",
        type=str,
        default=DEFAULT_FIG_PRESET,
        help="Figure size preset from _shared.style.SIZE_PRESETS_MM",
    )
    parser.add_argument("--width-mm", type=float, default=None, help="Manual figure width in mm")
    parser.add_argument("--height-mm", type=float, default=None, help="Manual figure height in mm")
    parser.add_argument("--wspace", type=float, default=0.08, help="Subplot width spacing")
    return parser.parse_args()


def resolve_figsize(args: argparse.Namespace) -> tuple[float, float]:
    if (args.width_mm is None) ^ (args.height_mm is None):
        raise ValueError("width-mm and height-mm must be provided together")
    if args.width_mm is not None and args.height_mm is not None:
        return figsize_from_mm(args.width_mm, args.height_mm)
    return figsize_from_preset(args.fig_preset)


def main() -> None:
    args = parse_args()
    apply_reference_style()

    df = read_table(args.input)
    fig, axes = plt.subplots(1, len(PARAMETERS), figsize=resolve_figsize(args), sharey=True)

    if len(PARAMETERS) == 1:
        axes = [axes]

    required = set(PARAMETERS + ["protocol"])
    if not required.issubset(df.columns):
        draw_placeholder(
            axes[0],
            "Shape Core Parameters",
            "Required columns are missing in intermediate table.",
        )
        for ax in axes[1:]:
            ax.axis("off")
    else:
        work = df.copy()
        work["protocol"] = pd.Categorical(work["protocol"], categories=PROTOCOL_ORDER, ordered=True)
        work = work.sort_values("protocol")

        for i, (ax, param) in enumerate(zip(axes, PARAMETERS)):
            sns.boxplot(
                data=work,
                x=param,
                y="protocol",
                hue="protocol",
                palette=palette(),
                flierprops=FLIER_PROPS,
                legend=False,
                ax=ax,
                orient="h",
            )
            # Match reference-style micro tweaks: inward ticks + numeric minor ticks.
            ax.tick_params(axis="both", which="both", direction="in")
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.tick_params(axis="x", which="minor")
            # Keep inward ticks from touching boxes near plot edges.
            ax.margins(y=0.06, x=0.07)
            ax.set_xlabel(param)
            ax.set_ylabel("" if i > 0 else "")
            if i > 0:
                ax.tick_params(axis="y", which="both", length=0)
            else:
                ax.tick_params(axis="y", which="minor", length=0)

        fig.subplots_adjust(wspace=max(0.0, float(args.wspace)))

    save_figure(fig, args.output)
    plt.close(fig)
    print(f"Wrote figure: {args.output}")


if __name__ == "__main__":
    main()
