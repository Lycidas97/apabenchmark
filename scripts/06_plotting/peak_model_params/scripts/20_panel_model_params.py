#!/usr/bin/env python3
"""Render protocol-wise custom model parameter panel."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import PROTOCOL_ORDER
from _shared.io import read_table
from _shared.paths import topic_figure_dir, topic_intermediate_dir
from _shared.style import apply_reference_style, draw_placeholder, figsize_from_preset, palette, save_figure

TOPIC = "peak_model_params"
INPUT_NAME = "peak_model_params_prepared.tsv"
OUTPUT_NAME = "model_params.pdf"
PARAMETERS = ["custom_mu1", "custom_sigma1", "custom_sigma2", "custom_a", "custom_b"]
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
    out_path = topic_figure_dir(TOPIC) / OUTPUT_NAME

    fig, axes = plt.subplots(1, len(PARAMETERS), figsize=figsize_from_preset("mid_wide"), sharey=True)

    if len(PARAMETERS) == 1:
        axes = [axes]

    required = set(PARAMETERS + ["protocol"])
    if not required.issubset(df.columns):
        draw_placeholder(
            axes[0],
            "Model Parameters",
            "Required columns are missing in intermediate table.",
        )
        for ax in axes[1:]:
            ax.axis("off")
    else:
        df = df[df["custom_a"] > -1000].copy() if "custom_a" in df.columns else df.copy()
        df["protocol"] = pd.Categorical(df["protocol"], categories=PROTOCOL_ORDER, ordered=True)
        df = df.sort_values("protocol")

        for i, (ax, param) in enumerate(zip(axes, PARAMETERS)):
            sns.boxplot(
                data=df,
                x=param,
                y="protocol",
                hue="protocol",
                palette=palette(),
                flierprops=FLIER_PROPS,
                legend=False,
                ax=ax,
                orient="h",
            )
            ax.set_xlabel(param)
            ax.set_ylabel("" if i > 0 else "")
            if i > 0:
                ax.tick_params(axis="y", which="both", length=0)

        fig.subplots_adjust(wspace=0.08)

    save_figure(fig, out_path)
    plt.close(fig)
    print(f"Wrote figure: {out_path}")


if __name__ == "__main__":
    main()
