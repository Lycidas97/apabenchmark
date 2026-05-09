#!/usr/bin/env python3
"""Render raw-data recall boxplot panel."""

from __future__ import annotations

import json
import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import PALETTE_HEX, TOOL_ORDER
from _shared.io import read_table
from _shared.paths import topic_figure_dir, topic_intermediate_dir, topic_root
from _shared.style import apply_reference_style, figsize_from_preset, save_figure

TOPIC = "raw_data_performance"
INPUT_NAME = "raw_data_recall_gt_filtered.parquet"
FILTER_CONFIG_PATH = topic_root(TOPIC) / "config" / "recall_gt_filter.json"
PLOT_CONFIG_PATH = topic_root(TOPIC) / "config" / "plot_panels.json"
OUTPUT_NAME = "raw_recall_boxplot.pdf"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render raw recall boxplot panel.")
    parser.add_argument(
        "--filter-config",
        default=str(FILTER_CONFIG_PATH),
        help="Path to recall_gt_filter.json.",
    )
    parser.add_argument(
        "--plot-config",
        default=str(PLOT_CONFIG_PATH),
        help="Path to plot_panels.json.",
    )
    return parser.parse_args()


def load_filter_pair_order(filter_config_path: Path) -> list[str]:
    cfg = json.loads(filter_config_path.read_text(encoding="utf-8"))
    out = []
    for item in cfg.get("filter_pairs", []):
        f1 = str(item.get("filter_type_1", "")).strip()
        f2 = str(item.get("filter_type_2", "")).strip()
        if f1 and f2:
            out.append(f"{f1} | {f2}")
    return out


def _with_suffix(file_name: str, suffix: str) -> str:
    p = Path(file_name)
    return f"{p.stem}_{suffix}{p.suffix}"


def _species_mask(df, species: str):
    if "sample" not in df.columns:
        return None
    sample = df["sample"].astype(str).str.lower()
    return sample.str.contains(str(species).lower(), na=False)


def _render_boxplot(
    *,
    df,
    filter_pair_order: list[str],
    figsize_preset: str,
    ymin: float,
    ymax: float,
    x_label: str,
    y_label: str,
    x_tick_rotation: float,
    legend_title: str,
    output_name: str,
) -> None:
    present_tools = set(df["tool"].dropna().astype(str).unique())
    tool_order = [name for name in TOOL_ORDER if name in present_tools]
    tool_order.extend(sorted(name for name in present_tools if name not in tool_order))

    present_pairs = set(df["filter_pair"].dropna().astype(str).unique())
    local_pair_order = [name for name in filter_pair_order if name in present_pairs]
    local_pair_order.extend(sorted(name for name in present_pairs if name not in local_pair_order))

    fig, ax = plt.subplots(figsize=figsize_from_preset(figsize_preset))
    palette = {name: PALETTE_HEX[idx % len(PALETTE_HEX)] for idx, name in enumerate(local_pair_order)}

    sns.boxplot(
        data=df,
        x="tool",
        y="recall",
        hue="filter_pair",
        order=tool_order,
        hue_order=local_pair_order,
        palette=palette,
        width=0.75,
        dodge=True,
        linewidth=0.5,
        flierprops={
            "marker": "o",
            "markersize": 1.0,
            "markerfacecolor": "gray",
            "markeredgecolor": "gray",
        },
        ax=ax,
    )

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_ylim(ymin, ymax)
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.6)
    ax.tick_params(axis="x", rotation=x_tick_rotation)
    ax.legend(
        title=legend_title,
        frameon=False,
        fontsize=6,
        title_fontsize=6,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0,
    )

    out_path = topic_figure_dir(TOPIC) / output_name
    save_figure(fig, out_path)
    plt.close(fig)
    print(f"Wrote figure: {out_path}")


def main() -> None:
    args = parse_args()
    apply_reference_style()

    plot_cfg = json.loads(Path(args.plot_config).read_text(encoding="utf-8"))
    box_cfg = plot_cfg.get("boxplot", {})

    input_name = str(box_cfg.get("input_name", INPUT_NAME))
    output_name = str(box_cfg.get("output_name", OUTPUT_NAME))
    figsize_preset = str(box_cfg.get("figsize_preset", "full_wide"))
    ymin = float(box_cfg.get("ymin", 0.0))
    ymax = float(box_cfg.get("ymax", 1.02))
    x_label = str(box_cfg.get("x_label", "Tool"))
    y_label = str(box_cfg.get("y_label", "Raw Recall"))
    x_tick_rotation = float(box_cfg.get("x_tick_rotation", 0))
    legend_title = str(box_cfg.get("legend_title", "Filter Pair"))

    path = topic_intermediate_dir(TOPIC) / input_name
    if not path.exists():
        raise FileNotFoundError(f"Missing filtered table: {path}. Run 01_prepare_recall_gt_data.py first.")

    df = read_table(path)
    required = {"tool", "recall", "filter_pair"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    df = df.copy()
    df["recall"] = df["recall"].astype(float)
    df = df[df["recall"].notna()]
    if df.empty:
        raise ValueError("No recall rows available after filtering.")

    filter_pair_order = load_filter_pair_order(Path(args.filter_config))
    _render_boxplot(
        df=df,
        filter_pair_order=filter_pair_order,
        figsize_preset=figsize_preset,
        ymin=ymin,
        ymax=ymax,
        x_label=x_label,
        y_label=y_label,
        x_tick_rotation=x_tick_rotation,
        legend_title=legend_title,
        output_name=output_name,
    )

    for species in ["human", "mouse"]:
        mask = _species_mask(df, species)
        if mask is None:
            print(f"Skip species split ({species}): missing 'sample' column.")
            break
        sub = df[mask].copy()
        if sub.empty:
            print(f"Skip species split ({species}): no matched rows.")
            continue
        _render_boxplot(
            df=sub,
            filter_pair_order=filter_pair_order,
            figsize_preset=figsize_preset,
            ymin=ymin,
            ymax=ymax,
            x_label=x_label,
            y_label=y_label,
            x_tick_rotation=x_tick_rotation,
            legend_title=legend_title,
            output_name=_with_suffix(output_name, species),
        )


if __name__ == "__main__":
    main()
