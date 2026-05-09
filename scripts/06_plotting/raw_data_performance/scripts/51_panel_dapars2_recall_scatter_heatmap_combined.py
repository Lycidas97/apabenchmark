#!/usr/bin/env python3
"""Render combined DaPars2 recall scatter and pairwise win-rate heatmap."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import TOOL_ORDER_CP
from _shared.io import read_table
from _shared.paths import topic_figure_dir, topic_intermediate_dir, topic_root
from _shared.style import apply_reference_style, figsize_from_mm, tool_palette

TOPIC = "raw_data_performance"
PLOT_CONFIG_PATH = topic_root(TOPIC) / "config" / "plot_panels.json"
DEFAULT_PANEL_KEY = "dapars2_recall_scatter_heatmap_combined"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render combined DaPars2 recall scatter and win-rate heatmap.")
    parser.add_argument("--plot-config", default=str(PLOT_CONFIG_PATH), help="Path to plot_panels.json.")
    parser.add_argument("--panel-key", default=DEFAULT_PANEL_KEY, help="Config key in plot_panels.json.")
    return parser.parse_args()


def _mm_rect(x_mm: float, y_mm: float, w_mm: float, h_mm: float, width_mm: float, height_mm: float) -> list[float]:
    return [x_mm / width_mm, y_mm / height_mm, w_mm / width_mm, h_mm / height_mm]


def _ordered_tools(present_tools: set[str], configured: list[str] | None = None) -> list[str]:
    base = configured or TOOL_ORDER_CP
    ordered = [t for t in base if t in present_tools]
    ordered.extend(sorted(t for t in present_tools if t not in ordered))
    return ordered


def _panel_label(label: str) -> str:
    parts = str(label).split("|", 1)
    if len(parts) != 2:
        return str(label)
    return f"{parts[0].strip()}\n{parts[1].strip()}"


def _load_boxplot_points(name: str) -> pd.DataFrame:
    path = topic_intermediate_dir(TOPIC) / name
    if not path.exists():
        raise FileNotFoundError(f"Missing boxplot points table: {path}")
    df = read_table(path)
    required = {"tool", "recall"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise KeyError(f"Missing required boxplot columns: {missing}")
    df = df.copy()
    df["recall"] = pd.to_numeric(df["recall"], errors="coerce")
    df = df[df["recall"].notna()].copy()
    if "filter_pair" not in df.columns:
        df["filter_pair"] = "pipeline_filter"
    return df


def _load_scatter_points(name: str) -> pd.DataFrame:
    path = topic_intermediate_dir(TOPIC) / name
    if not path.exists():
        raise FileNotFoundError(f"Missing scatter points table: {path}")
    df = read_table(path)
    required = {
        "tool",
        "raw_recall_mean",
        "raw_recall_sd",
        "sim_f1_mean",
        "sim_f1_sd",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise KeyError(f"Missing required scatter columns: {missing}")
    df = df.copy()
    for col in ["raw_recall_mean", "raw_recall_sd", "sim_f1_mean", "sim_f1_sd"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[df["raw_recall_mean"].notna() & df["sim_f1_mean"].notna()].copy()
    if df.empty:
        raise ValueError(f"No usable scatter points in: {path}")
    return df


def _load_heatmap_points(name: str) -> pd.DataFrame:
    path = topic_intermediate_dir(TOPIC) / name
    if not path.exists():
        raise FileNotFoundError(f"Missing heatmap points table: {path}")
    df = read_table(path)
    required = {"tool_y", "tool_x", "win_rate"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise KeyError(f"Missing required heatmap columns: {missing}")
    df = df.copy()
    df["win_rate"] = pd.to_numeric(df["win_rate"], errors="coerce")
    if "filter_group" not in df.columns:
        df["filter_group"] = "all"
    return df[df["win_rate"].notna()].copy()


def _draw_boxplot(
    fig: plt.Figure,
    *,
    rect: list[float],
    df: pd.DataFrame,
    tool_order: list[str],
    y_label: str,
    font_size: float,
) -> None:
    ax = fig.add_axes(rect)
    sns.boxplot(
        data=df,
        x="recall",
        y="tool",
        hue="tool",
        order=tool_order,
        hue_order=tool_order,
        palette=tool_palette(tool_order),
        dodge=False,
        width=0.58,
        linewidth=0.55,
        flierprops={
            "marker": "o",
            "markersize": 1.0,
            "markerfacecolor": "gray",
            "markeredgecolor": "gray",
            "markeredgewidth": 0.2,
        },
        ax=ax,
    )
    leg = ax.get_legend()
    if leg:
        leg.remove()
    ax.set_xlabel(y_label, fontsize=font_size)
    ax.set_ylabel("")
    ax.set_xlim(0, 1.02)
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.tick_params(axis="x", labelsize=font_size, pad=2)
    ax.tick_params(axis="y", labelsize=font_size, length=0, pad=2)
    ax.grid(axis="x", linestyle="--", linewidth=0.4, alpha=0.55)


def _draw_scatter(
    fig: plt.Figure,
    *,
    rect: list[float],
    df: pd.DataFrame,
    tool_order: list[str],
    x_label: str,
    y_label: str,
    font_size: float,
    marker_size: float,
    capsize: float,
    show_legend: bool,
    show_correlation: bool,
) -> None:
    ax = fig.add_axes(rect)
    colors = tool_palette(tool_order)
    points = df.copy()
    points["tool"] = pd.Categorical(points["tool"].astype(str), categories=tool_order, ordered=True)
    points = points.sort_values("tool").reset_index(drop=True)

    for _, row in points.iterrows():
        tool = str(row["tool"])
        color = colors.get(tool, "#929292")
        ax.errorbar(
            float(row["raw_recall_mean"]),
            float(row["sim_f1_mean"]),
            xerr=float(row["raw_recall_sd"]) if pd.notna(row["raw_recall_sd"]) else 0.0,
            yerr=float(row["sim_f1_sd"]) if pd.notna(row["sim_f1_sd"]) else 0.0,
            fmt="o",
            color=color,
            ecolor=color,
            elinewidth=0.7,
            capsize=capsize,
            markersize=marker_size,
            alpha=0.95,
        )

    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel(x_label, fontsize=font_size)
    ax.set_ylabel(y_label, fontsize=font_size)
    ax.tick_params(axis="both", labelsize=font_size, pad=2)
    ax.grid(True, linestyle=":", linewidth=0.45, alpha=0.7)

    if show_correlation and len(points) >= 2:
        x_ser = pd.to_numeric(points["raw_recall_mean"], errors="coerce")
        y_ser = pd.to_numeric(points["sim_f1_mean"], errors="coerce")
        valid = x_ser.notna() & y_ser.notna()
        if int(valid.sum()) >= 2:
            r = x_ser[valid].corr(y_ser[valid], method="spearman")
            if pd.notna(r):
                ax.text(
                    0.04,
                    0.96,
                    f"Spearman\nr = {r:.2f}",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=font_size,
                )

    if show_legend:
        handles = [
            Line2D([0], [0], marker="o", linestyle="", color=colors.get(tool, "#929292"), label=tool, markersize=3)
            for tool in tool_order
            if tool in set(points["tool"].astype(str))
        ]
        ax.legend(
            handles=handles,
            loc="lower right",
            frameon=False,
            fontsize=max(5.5, font_size - 1.0),
            handlelength=0.8,
            handletextpad=0.3,
            borderaxespad=0.2,
            labelspacing=0.25,
        )


def _matrix_from_points(df: pd.DataFrame, tool_order: list[str], group: str) -> np.ndarray:
    gdf = df[df["filter_group"].astype(str) == str(group)].copy()
    mat = pd.DataFrame(np.nan, index=tool_order, columns=tool_order, dtype=float)
    for _, row in gdf.iterrows():
        y_tool = str(row["tool_y"])
        x_tool = str(row["tool_x"])
        if y_tool in mat.index and x_tool in mat.columns:
            mat.loc[y_tool, x_tool] = float(row["win_rate"])
    return mat.to_numpy(dtype=float)


def _draw_heatmap(
    fig: plt.Figure,
    *,
    x0_mm: float,
    y0_mm: float,
    width_mm: float,
    height_mm: float,
    figure_width_mm: float,
    figure_height_mm: float,
    df: pd.DataFrame,
    tool_order: list[str],
    cmap: str,
    vmin: float,
    vmax: float,
    font_size: float,
    title_font_size: float,
    group_label_h_mm: float,
    label_h_mm: float,
    y_label_w_mm: float,
    panel_gap_mm: float,
    cbar_label_w_mm: float,
    cbar_width_mm: float,
    cbar_gap_mm: float,
    cbar_label: str,
    show_y_axis_label: bool,
) -> None:
    groups = [str(g) for g in df["filter_group"].dropna().astype(str).unique()]
    if not groups:
        raise ValueError("No heatmap groups to draw.")

    cbar_total_mm = cbar_gap_mm + cbar_width_mm + cbar_label_w_mm
    heat_area_w_mm = width_mm - y_label_w_mm - cbar_total_mm
    cell_size_mm = min(
        (heat_area_w_mm - panel_gap_mm * (len(groups) - 1)) / len(groups),
        height_mm - group_label_h_mm - label_h_mm - 0.8,
    )
    heat_block_w_mm = len(groups) * cell_size_mm + panel_gap_mm * (len(groups) - 1)
    left_pad_mm = max(0.0, (heat_area_w_mm - heat_block_w_mm) / 2)
    cmap_obj = plt.get_cmap(cmap)
    norm = Normalize(vmin=vmin, vmax=vmax)

    for idx, group in enumerate(groups):
        mat = _matrix_from_points(df, tool_order, group)
        np.fill_diagonal(mat, np.nan)
        px_mm = x0_mm + y_label_w_mm + left_pad_mm + idx * (cell_size_mm + panel_gap_mm)
        py_mm = y0_mm + label_h_mm
        ax = fig.add_axes(_mm_rect(px_mm, py_mm, cell_size_mm, cell_size_mm, figure_width_mm, figure_height_mm))
        ax.set_xlim(-0.5, len(tool_order) - 0.5)
        ax.set_ylim(len(tool_order) - 0.5, -0.5)
        ax.set_facecolor("white")
        for iy in range(len(tool_order)):
            for ix in range(len(tool_order)):
                val = mat[iy, ix]
                facecolor = "white" if np.isnan(val) else cmap_obj(norm(val))
                ax.add_patch(
                    Rectangle(
                        (ix - 0.5, iy - 0.5),
                        1.0,
                        1.0,
                        facecolor=facecolor,
                        edgecolor="white",
                        linewidth=0.55,
                    )
                )
        ax.set_xticks(range(len(tool_order)))
        ax.set_yticks(range(len(tool_order)))
        ax.set_xticklabels(tool_order, rotation=90, fontsize=font_size)
        ax.set_yticklabels(tool_order if idx == 0 else [], fontsize=font_size)
        ax.tick_params(axis="x", length=0, pad=2)
        ax.tick_params(axis="y", length=0, pad=2)
        for spine in ax.spines.values():
            spine.set_visible(False)
        for iy in range(len(tool_order)):
            for ix in range(len(tool_order)):
                val = mat[iy, ix]
                if np.isnan(val):
                    continue
                ax.text(ix, iy, f"{val:.2f}", ha="center", va="center", fontsize=font_size, color="#202020")

        title_ax = fig.add_axes(
            _mm_rect(px_mm, py_mm + cell_size_mm + 0.8, cell_size_mm, group_label_h_mm, figure_width_mm, figure_height_mm)
        )
        title_ax.axis("off")
        title_ax.text(0.5, 0.0, _panel_label(group), ha="center", va="bottom", fontsize=title_font_size)

    if show_y_axis_label:
        label_ax = fig.add_axes(_mm_rect(x0_mm, y0_mm + label_h_mm, y_label_w_mm, cell_size_mm, figure_width_mm, figure_height_mm))
        label_ax.axis("off")
        label_ax.text(0.02, 0.5, "Y-axis tool beats X-axis tool", rotation=90, ha="left", va="center", fontsize=font_size)

    cax_x_mm = x0_mm + y_label_w_mm + left_pad_mm + heat_block_w_mm + cbar_gap_mm
    cax_h_mm = cell_size_mm * 0.74
    cax_y_mm = y0_mm + label_h_mm + (cell_size_mm - cax_h_mm) / 2
    cax = fig.add_axes(_mm_rect(cax_x_mm, cax_y_mm, cbar_width_mm, cax_h_mm, figure_width_mm, figure_height_mm))
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.ax.tick_params(labelsize=font_size, pad=1)
    cbar.set_label(cbar_label, fontsize=font_size, labelpad=2)


def main() -> None:
    args = parse_args()
    apply_reference_style()

    cfg = json.loads(Path(args.plot_config).read_text(encoding="utf-8"))
    panel = cfg.get(args.panel_key, {})
    if not panel:
        raise KeyError(f"Missing panel config: {args.panel_key}")

    width_mm = float(panel.get("width_mm", 178))
    height_mm = float(panel.get("height_mm", 82))
    side_margin_mm = float(panel.get("side_margin_mm", 3))
    left_margin_mm = float(panel.get("left_margin_mm", side_margin_mm))
    right_margin_mm = float(panel.get("right_margin_mm", side_margin_mm))
    panel_gap_mm = float(panel.get("panel_gap_mm", 4))
    scatter_width_frac = float(panel.get("scatter_width_frac", panel.get("boxplot_width_frac", 0.42)))
    output_figure = str(panel.get("output_figure", "raw_dapars2_pipeline_recall_scatter_heatmap_combined.pdf"))

    scatter_df = _load_scatter_points(
        str(panel.get("scatter_points", "raw_dapars2_pipeline_recall_vs_sim_f1_scatter_points.tsv"))
    )
    heat_df = _load_heatmap_points(str(panel.get("heatmap_points", "raw_dapars2_pipeline_recall_pairwise_winrate.tsv")))

    configured_order = [str(t) for t in panel.get("tool_order", TOOL_ORDER_CP)]
    present = set(scatter_df["tool"].dropna().astype(str)) | set(heat_df["tool_y"].dropna().astype(str)) | set(heat_df["tool_x"].dropna().astype(str))
    tool_order = _ordered_tools(present, configured_order)

    fig = plt.figure(figsize=figsize_from_mm(width_mm, height_mm))
    content_w_mm = width_mm - left_margin_mm - right_margin_mm
    scatter_w_mm = content_w_mm * scatter_width_frac
    heat_w_mm = content_w_mm - scatter_w_mm - panel_gap_mm
    y0_mm = float(panel.get("bottom_margin_mm", 5))
    content_h_mm = height_mm - y0_mm - float(panel.get("top_margin_mm", 4))
    heat_group_label_h_mm = float(panel.get("heatmap_group_label_height_mm", 7))
    heatmap_y_label_w_mm = float(panel.get("heatmap_y_label_width_mm", 12))
    heatmap_colorbar_label_w_mm = float(panel.get("heatmap_colorbar_label_width_mm", 8))
    heatmap_colorbar_gap_mm = float(panel.get("colorbar_gap_mm", 3))
    heatmap_colorbar_width_mm = float(panel.get("colorbar_width_mm", 4))
    heat_groups = max(1, int(heat_df["filter_group"].dropna().astype(str).nunique()))
    heatmap_panel_gap_mm = float(panel.get("heatmap_panel_gap_mm", 2))
    heatmap_label_h_mm = float(panel.get("heatmap_x_label_height_mm", 8))
    heatmap_title_gap_mm = float(panel.get("heatmap_title_gap_mm", 0.8))
    heat_area_w_mm = heat_w_mm - heatmap_y_label_w_mm - heatmap_colorbar_gap_mm - heatmap_colorbar_width_mm - heatmap_colorbar_label_w_mm
    shared_main_h_mm = min(
        (heat_area_w_mm - heatmap_panel_gap_mm * (heat_groups - 1)) / heat_groups,
        content_h_mm - heatmap_label_h_mm - heatmap_title_gap_mm - heat_group_label_h_mm,
    )
    scatter_main_h_mm = min(
        content_h_mm,
        shared_main_h_mm + float(panel.get("scatter_height_extra_mm", 0)),
    )

    _draw_scatter(
        fig,
        rect=_mm_rect(left_margin_mm, y0_mm, scatter_w_mm, scatter_main_h_mm, width_mm, height_mm),
        df=scatter_df,
        tool_order=tool_order,
        x_label=str(panel.get("scatter_x_label", "Raw Recall (DaPars2 pipeline)")),
        y_label=str(panel.get("scatter_y_label", "Sim F1")),
        font_size=float(panel.get("font_size", 7)),
        marker_size=float(panel.get("scatter_marker_size", 3.0)),
        capsize=float(panel.get("scatter_capsize", 1.4)),
        show_legend=bool(panel.get("scatter_show_legend", True)),
        show_correlation=bool(panel.get("scatter_show_correlation", True)),
    )

    _draw_heatmap(
        fig,
        x0_mm=left_margin_mm + scatter_w_mm + panel_gap_mm,
        y0_mm=y0_mm,
        width_mm=heat_w_mm,
        height_mm=content_h_mm,
        figure_width_mm=width_mm,
        figure_height_mm=height_mm,
        df=heat_df,
        tool_order=tool_order,
        cmap=str(panel.get("cmap", "RdBu")),
        vmin=float(panel.get("vmin", 0.0)),
        vmax=float(panel.get("vmax", 1.0)),
        font_size=float(panel.get("font_size", 7)),
        title_font_size=float(panel.get("title_font_size", 7)),
        group_label_h_mm=heat_group_label_h_mm,
        label_h_mm=heatmap_label_h_mm,
        y_label_w_mm=heatmap_y_label_w_mm,
        panel_gap_mm=heatmap_panel_gap_mm,
        cbar_label_w_mm=heatmap_colorbar_label_w_mm,
        cbar_width_mm=float(panel.get("colorbar_width_mm", 4)),
        cbar_gap_mm=float(panel.get("colorbar_gap_mm", 3)),
        cbar_label=str(panel.get("colorbar_label", "Win rate")),
        show_y_axis_label=bool(panel.get("heatmap_show_y_axis_label", False)),
    )

    out_fig = topic_figure_dir(TOPIC) / output_figure
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, dpi=int(panel.get("render_dpi", 450)))
    plt.close(fig)
    print(f"Wrote combined DaPars2 scatter/heatmap figure: {out_fig}")
    print(f"Combined size target: {width_mm:.1f}mm x {height_mm:.1f}mm")


if __name__ == "__main__":
    main()
