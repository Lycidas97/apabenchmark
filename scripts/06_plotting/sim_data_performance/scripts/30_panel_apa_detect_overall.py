#!/usr/bin/env python3
"""Render APA detect overall panel for sim data performance."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.paths import topic_figure_dir
from _shared.style import apply_reference_style, save_figure

from _plot_helpers import (
    TOPIC,
    apply_pd_hatching,
    available_tool_order,
    black_boxplot_props,
    build_apa_plot_table,
    build_hatched_legend_handles,
    match_type_palette,
    resolve_match_type_order,
    trim_metric_to_quantile,
)
from _plot_config import apply_runtime_rcparams, cfg, cfg_bool, cfg_float, figsize_for_figure

OUTPUT_NAME = "apa_detect_performance.pdf"
METRICS = ["precision", "recall", "f1"]


def main() -> None:
    apply_reference_style()
    apply_runtime_rcparams()
    df = build_apa_plot_table()

    hue_order = resolve_match_type_order(df["match_type"].dropna().astype(str).tolist())
    tool_order = available_tool_order(df)

    fig, axes = plt.subplots(3, 1, figsize=figsize_for_figure(OUTPUT_NAME, "full_wide"))
    fig.subplots_adjust(hspace=cfg_float("layout.overall_hspace", 0.25))

    for idx, metric in enumerate(METRICS):
        ax = axes[idx]
        plot_data = trim_metric_to_quantile(df, metric)
        sns.boxplot(
            data=plot_data,
            x="tool",
            y=metric,
            hue="match_type",
            hue_order=hue_order,
            order=tool_order,
            palette=match_type_palette(hue_order),
            dodge=True,
            width=cfg_float("boxplot.overall_width", 0.8),
            ax=ax,
            showfliers=cfg_bool("boxplot.showfliers", False),
            **black_boxplot_props(cfg_float("boxplot.overall_linewidth", 0.5)),
        )
        apply_pd_hatching(ax, hue_order)
        ax.grid(
            linestyle=str(cfg("grid.overall_linestyle", "--")),
            alpha=cfg_float("grid.overall_alpha", 0.6),
            linewidth=cfg_float("grid.overall_linewidth", 0.5),
        )
        ax.set_ylabel(metric.capitalize())

        if idx == 0:
            labels = [name for name in hue_order]
            ax.legend(
                build_hatched_legend_handles(hue_order),
                labels,
                title="Match Type",
                title_fontsize=8,
                fontsize=8,
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
