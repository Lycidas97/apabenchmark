#!/usr/bin/env python3
"""Render APA detect metric-by-protocol panels for sim data performance."""

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
    available_protocol_order,
    available_tool_order,
    black_boxplot_props,
    build_apa_plot_table,
    build_hatched_legend_handles,
    match_type_palette,
    resolve_match_type_order,
    trim_metric_to_quantile,
)
from _plot_config import (
    apply_runtime_rcparams,
    cfg,
    cfg_bool,
    cfg_float,
    cfg_int,
    figsize_for_figure,
)

METRICS = ["precision", "recall", "f1"]


def render_metric(df, metric: str) -> None:
    output_name = f"apa_detect_{metric}_by_protocol.pdf"
    protocol_order = available_protocol_order(df)
    hue_order = resolve_match_type_order(df["match_type"].dropna().astype(str).tolist())
    tool_order = available_tool_order(df)

    fig, axes = plt.subplots(
        nrows=len(protocol_order),
        ncols=1,
        figsize=figsize_for_figure(output_name, "narrow_tall"),
        dpi=cfg_int("render.per_figure_dpi", 300),
    )
    if len(protocol_order) == 1:
        axes = [axes]

    for idx, (protocol, ax) in enumerate(zip(protocol_order, axes)):
        plot_data = trim_metric_to_quantile(df[df["protocol"] == protocol], metric)
        sns.boxplot(
            data=plot_data,
            x="tool",
            y=metric,
            hue="match_type",
            hue_order=hue_order,
            order=tool_order,
            palette=match_type_palette(hue_order),
            ax=ax,
            dodge=True,
            width=cfg_float("boxplot.by_protocol_width", 0.7),
            showfliers=cfg_bool("boxplot.showfliers", False),
            **black_boxplot_props(cfg_float("boxplot.by_protocol_linewidth", 0.4)),
        )
        apply_pd_hatching(ax, hue_order)

        ax.set_title(protocol, loc="left", pad=2, fontsize=8)
        ax.set_ylabel(metric.capitalize() if idx == 0 else "", labelpad=2)
        ax.set_xlabel("")
        ax.tick_params(axis="both", labelsize=8, pad=1)
        ax.grid(
            True,
            linestyle=str(cfg("grid.by_protocol_linestyle", ":")),
            linewidth=cfg_float("grid.by_protocol_linewidth", 0.3),
        )

        if idx != len(protocol_order) - 1:
            ax.set_xticklabels([])
        else:
            plt.setp(ax.get_xticklabels(), rotation=90)

        leg = ax.get_legend()
        if leg:
            leg.remove()

    fig.legend(
        build_hatched_legend_handles(hue_order),
        hue_order,
        loc="upper center",
        bbox_to_anchor=(0.5, cfg_float("layout.apa_by_protocol_legend_y", 0.93)),
        ncol=max(1, len(hue_order)),
        title="Match Type",
        frameon=False,
        title_fontsize=9,
        fontsize=8,
        columnspacing=0.8,
        handletextpad=0.3,
    )

    fig.subplots_adjust(
        hspace=cfg_float("layout.by_protocol_hspace", 0.25),
        left=cfg_float("layout.by_protocol_left", 0.12),
        right=cfg_float("layout.by_protocol_right", 0.88),
        top=cfg_float("layout.by_protocol_top", 0.82),
        bottom=cfg_float("layout.by_protocol_bottom", 0.08),
    )
    out_path = topic_figure_dir(TOPIC) / output_name
    save_figure(fig, out_path)
    plt.close(fig)
    print(f"Wrote figure: {out_path}")


def main() -> None:
    apply_reference_style()
    apply_runtime_rcparams()
    df = build_apa_plot_table()
    for metric in METRICS:
        render_metric(df, metric)


if __name__ == "__main__":
    main()
