#!/usr/bin/env python3
"""Render recall boxplot for dedicated raw DaPars2/scMAPA-aligned pipeline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import PALETTE_HEX, TOOL_ORDER_CP, TOOL_MAP
from _shared.io import read_table, write_tsv
from _shared.paths import topic_figure_dir, topic_intermediate_dir, topic_root
from _shared.style import apply_reference_style, figsize_from_preset, save_figure, tool_palette

TOPIC = "raw_data_performance"
DEFAULT_INPUT = topic_intermediate_dir(TOPIC) / "raw_data_performance_prepared.parquet"
DEFAULT_FILTER_CONFIG = topic_root(TOPIC) / "config" / "recall_gt_filter.json"
DEFAULT_DAPARS2_CONFIG = topic_root(TOPIC) / "config" / "dapars2_pipeline.json"
DEFAULT_PLOT_CONFIG = topic_root(TOPIC) / "config" / "plot_panels.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render recall boxplot for dedicated dapars2/scmapa-aligned pipeline outputs."
    )
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Prepared metrics parquet path.")
    parser.add_argument("--filter-config", default=str(DEFAULT_FILTER_CONFIG), help="recall_gt_filter config path.")
    parser.add_argument("--dapars2-config", default=str(DEFAULT_DAPARS2_CONFIG), help="dapars2_pipeline config path.")
    parser.add_argument("--plot-config", default=str(DEFAULT_PLOT_CONFIG), help="plot_panels config path.")
    return parser.parse_args()


def _load_selected_tools(dapars2_cfg: dict) -> list[str] | None:
    mode = str(dapars2_cfg.get("tools_mode", "selected")).strip().lower()
    if mode == "all":
        return None
    tools = [str(v).strip() for v in dapars2_cfg.get("selected_tools", []) if str(v).strip()]
    return tools if tools else None


def _with_suffix(file_name: str, suffix: str) -> str:
    p = Path(file_name)
    return f"{p.stem}_{suffix}{p.suffix}"


def _species_mask(df: pd.DataFrame, species: str) -> pd.Series | None:
    if "sample" not in df.columns:
        return None
    sample = df["sample"].astype(str).str.lower()
    return sample.str.contains(str(species).lower(), na=False)


def _render_boxplot(
    *,
    df: pd.DataFrame,
    figsize_preset: str,
    ymin: float,
    ymax: float,
    x_label: str,
    y_label: str,
    x_tick_rotation: float,
    show_filter_hue: bool,
    single_filter_color_by_tool: bool,
    output_name: str,
) -> None:
    present_tools = set(df["tool"].dropna().astype(str).unique())
    tool_order = [t for t in TOOL_ORDER_CP if t in present_tools]
    tool_order.extend(sorted(t for t in present_tools if t not in tool_order))

    fig, ax = plt.subplots(figsize=figsize_from_preset(figsize_preset))

    if show_filter_hue and df["filter_pair"].nunique() > 1:
        hue_order = sorted(df["filter_pair"].dropna().astype(str).unique())
        palette = {name: PALETTE_HEX[i % len(PALETTE_HEX)] for i, name in enumerate(hue_order)}
        sns.boxplot(
            data=df,
            x="tool",
            y="recall",
            hue="filter_pair",
            order=tool_order,
            hue_order=hue_order,
            palette=palette,
            width=0.75,
            dodge=True,
            linewidth=0.5,
            flierprops={"marker": "o", "markersize": 1.0, "markerfacecolor": "gray", "markeredgecolor": "gray"},
            ax=ax,
        )
        ax.legend(title="Filter Pair", frameon=False, fontsize=6, title_fontsize=6, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    else:
        if single_filter_color_by_tool:
            sns.boxplot(
                data=df,
                x="tool",
                y="recall",
                hue="tool",
                order=tool_order,
                hue_order=tool_order,
                palette=tool_palette(tool_order),
                dodge=False,
                width=0.6,
                linewidth=0.5,
                flierprops={"marker": "o", "markersize": 1.0, "markerfacecolor": "gray", "markeredgecolor": "gray"},
                ax=ax,
            )
        else:
            sns.boxplot(
                data=df,
                x="tool",
                y="recall",
                order=tool_order,
                color=PALETTE_HEX[0],
                width=0.6,
                linewidth=0.5,
                flierprops={"marker": "o", "markersize": 1.0, "markerfacecolor": "gray", "markeredgecolor": "gray"},
                ax=ax,
            )
        leg = ax.get_legend()
        if leg:
            leg.remove()

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_ylim(ymin, ymax)
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.6)
    ax.tick_params(axis="x", rotation=x_tick_rotation)

    out_fig = topic_figure_dir(TOPIC) / output_name
    save_figure(fig, out_fig)
    plt.close(fig)
    print(f"Wrote figure: {out_fig}")


def main() -> None:
    args = parse_args()
    apply_reference_style()

    input_path = Path(args.input).expanduser().resolve()
    filter_cfg = json.loads(Path(args.filter_config).expanduser().resolve().read_text(encoding="utf-8"))
    dapars2_cfg = json.loads(Path(args.dapars2_config).expanduser().resolve().read_text(encoding="utf-8"))
    plot_cfg = json.loads(Path(args.plot_config).expanduser().resolve().read_text(encoding="utf-8"))

    panel_cfg = plot_cfg.get("dapars2_pipeline_boxplot", {})
    output_name = str(panel_cfg.get("output_name", "raw_dapars2_pipeline_recall_boxplot.pdf"))
    output_points = str(panel_cfg.get("output_points", "raw_dapars2_pipeline_recall_boxplot_points.tsv"))
    figsize_preset = str(panel_cfg.get("figsize_preset", "full_wide"))
    ymin = float(panel_cfg.get("ymin", 0.0))
    ymax = float(panel_cfg.get("ymax", 1.02))
    x_label = str(panel_cfg.get("x_label", "Tool"))
    y_label = str(panel_cfg.get("y_label", "Raw Recall (DaPars2 pipeline)"))
    x_tick_rotation = float(panel_cfg.get("x_tick_rotation", 0))
    show_filter_hue = bool(panel_cfg.get("show_filter_hue", True))
    single_filter_color_by_tool = bool(panel_cfg.get("single_filter_color_by_tool", True))

    if not input_path.exists():
        raise FileNotFoundError(f"Missing prepared metrics table: {input_path}")

    df = read_table(input_path)
    required = {"metric_domain", "tool", "recall", "gt_count", "sample_id", "pair_label"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    df = df[df["metric_domain"] == "raw_dapars2_performance"].copy()
    if df.empty:
        raise ValueError("No rows for metric_domain='raw_dapars2_performance'.")

    selected_tools = _load_selected_tools(dapars2_cfg)
    if selected_tools is not None:
        selected_disp = {TOOL_MAP.get(t, t) for t in selected_tools}
        df = df[df["tool"].astype(str).isin(selected_disp)].copy()

    df["recall"] = pd.to_numeric(df["recall"], errors="coerce")
    df["gt_count"] = pd.to_numeric(df["gt_count"], errors="coerce")
    df = df[df["recall"].notna() & df["gt_count"].notna()].copy()

    # Align with current filtering strategy for comparability.
    gt_rules = filter_cfg.get("gt_rules", {})
    recall_rules = filter_cfg.get("recall_rules", {})
    if bool(gt_rules.get("drop_gt_zero", True)):
        df = df[df["gt_count"] != 0].copy()
    min_gt = int(gt_rules.get("min_gt_count", 1))
    df = df[df["gt_count"] >= min_gt].copy()

    if bool(recall_rules.get("drop_all_tool_zero_recall_sample_pairs", False)):
        key_cols = [c for c in ["filter_type_1", "filter_type_2", "sample_id", "pair_label"] if c in df.columns]
        if {"sample_id", "pair_label"}.issubset(set(key_cols)):
            per_pair = (
                df.groupby(key_cols, dropna=False)["recall"]
                .max()
                .reset_index(name="max_recall")
            )
            zeros = per_pair[per_pair["max_recall"] == 0]
            if not zeros.empty:
                zero_keys = set(tuple(row[c] for c in key_cols) for _, row in zeros.iterrows())
                mask = df.apply(lambda r: tuple(r[c] for c in key_cols) not in zero_keys, axis=1)
                df = df[mask].copy()

    if df.empty:
        raise ValueError("No rows left after dapars2 pipeline filtering.")

    if {"filter_type_1", "filter_type_2"}.issubset(df.columns):
        df["filter_pair"] = df["filter_type_1"].astype(str) + " | " + df["filter_type_2"].astype(str)
    else:
        df["filter_pair"] = "pipeline_filter"

    out_points_path = topic_intermediate_dir(TOPIC) / output_points
    write_tsv(df, out_points_path)
    print(f"Wrote points table: {out_points_path}")

    _render_boxplot(
        df=df,
        figsize_preset=figsize_preset,
        ymin=ymin,
        ymax=ymax,
        x_label=x_label,
        y_label=y_label,
        x_tick_rotation=x_tick_rotation,
        show_filter_hue=show_filter_hue,
        single_filter_color_by_tool=single_filter_color_by_tool,
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
        out_points_species = topic_intermediate_dir(TOPIC) / _with_suffix(output_points, species)
        write_tsv(sub, out_points_species)
        print(f"Wrote points table: {out_points_species}")
        _render_boxplot(
            df=sub,
            figsize_preset=figsize_preset,
            ymin=ymin,
            ymax=ymax,
            x_label=x_label,
            y_label=y_label,
            x_tick_rotation=x_tick_rotation,
            show_filter_hue=show_filter_hue,
            single_filter_color_by_tool=single_filter_color_by_tool,
            output_name=_with_suffix(output_name, species),
        )


if __name__ == "__main__":
    main()
