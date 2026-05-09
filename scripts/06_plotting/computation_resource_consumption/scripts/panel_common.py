#!/usr/bin/env python3
"""Shared helpers for computation resource consumption panels."""

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

from _shared.constants import TOOL_ORDER_CP
from _shared.io import read_table
from _shared.paths import topic_figure_dir, topic_intermediate_dir
from _shared.style import apply_reference_style, figsize_from_mm, save_figure, tool_palette

TOPIC = "computation_resource_consumption"
INPUT_NAME = "cp_resource_consumption_prepared.tsv"
X_PANEL_SPECS = [
    ("peak_size", "Peak Size", [250, 500, 1000, 2000, 4000]),
    ("gene_number", "Gene Number", [1000, 2000, 4000, 8000, 16000]),
    ("barcode_number", "Barcode Number", [1000, 2000, 4000, 8000, 16000]),
    ("read_length", "Read Length", [40, 70, 100, 130, 160]),
]


def load_prepared_table() -> pd.DataFrame:
    input_path = topic_intermediate_dir(TOPIC) / INPUT_NAME
    if not input_path.exists():
        raise FileNotFoundError(
            f"Missing prepared table: {input_path}. "
            "Run 00_prepare_data.py first."
        )
    return read_table(input_path)


def resolve_tool_order(df: pd.DataFrame) -> list[str]:
    available = set(df["tool"].dropna().unique().tolist())
    ordered = [tool for tool in TOOL_ORDER_CP if tool in available]
    extras = sorted(available.difference(ordered))
    return ordered + extras


def infer_baseline_values(df: pd.DataFrame) -> dict[str, float | int]:
    baseline: dict[str, float | int] = {}
    for column, _, _ in X_PANEL_SPECS:
        counts = df[column].value_counts(dropna=True)
        if counts.empty:
            raise ValueError(f"Cannot infer baseline for {column}: no non-null values.")
        top_count = counts.iloc[0]
        top_values = counts[counts == top_count].index.tolist()
        if len(top_values) != 1:
            raise ValueError(
                f"Cannot infer baseline for {column}: tied most frequent values {top_values}."
            )
        baseline[column] = top_values[0]
    return baseline


def select_panel_df(
    df: pd.DataFrame,
    x_var: str,
    baseline: dict[str, float | int],
) -> pd.DataFrame:
    panel_df = df
    for column, _, _ in X_PANEL_SPECS:
        if column == x_var:
            continue
        panel_df = panel_df[panel_df[column] == baseline[column]]
    return panel_df


def render_metric_panel(df: pd.DataFrame, metric: str, y_label: str, output_name: str) -> None:
    apply_reference_style()
    tool_order = resolve_tool_order(df)
    color_map = tool_palette(tool_order)
    baseline = infer_baseline_values(df)

    fig, axes = plt.subplots(1, len(X_PANEL_SPECS), figsize=figsize_from_mm(180, 40), sharey=True)
    for idx, (x_var, x_label, x_ticks) in enumerate(X_PANEL_SPECS):
        ax = axes[idx]
        show_legend = idx == len(X_PANEL_SPECS) - 1
        panel_df = select_panel_df(df, x_var, baseline)
        if panel_df.empty:
            raise ValueError(
                f"No records available for panel {x_var} after applying baseline filter {baseline}."
            )
        sns.lineplot(
            data=panel_df,
            x=x_var,
            y=metric,
            hue="tool",
            hue_order=tool_order,
            palette=color_map,
            err_style="bars",
            linewidth=0.5,
            legend=show_legend,
            ax=ax,
        )
        ax.tick_params(which="minor", bottom=False)
        ax.set_xlabel(x_label)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticks, rotation=90)

        if show_legend:
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                ax.legend(
                    handles,
                    labels,
                    title="",
                    loc="center left",
                    bbox_to_anchor=(1, 0.5),
                    frameon=False,
                )
        else:
            leg = ax.get_legend()
            if leg:
                leg.remove()

    axes[0].set_ylabel(y_label)
    fig.subplots_adjust(wspace=0.05)

    out_path = topic_figure_dir(TOPIC) / output_name
    save_figure(fig, out_path)
    plt.close(fig)
    print(f"Wrote figure: {out_path}")
