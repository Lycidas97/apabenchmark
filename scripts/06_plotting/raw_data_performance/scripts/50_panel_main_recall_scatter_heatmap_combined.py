#!/usr/bin/env python3
"""Compose main recall scatter and heatmap panels into one width-limited figure."""

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

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import TOOL_ORDER_CP
from _shared.paths import topic_figure_dir, topic_intermediate_dir, topic_root
from _shared.style import apply_reference_style, figsize_from_mm, tool_palette

TOPIC = "raw_data_performance"
PLOT_CONFIG_PATH = topic_root(TOPIC) / "config" / "plot_panels.json"
DEFAULT_PANEL_KEY = "main_recall_scatter_heatmap_combined"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compose main recall scatter and heatmap into one PDF.")
    parser.add_argument(
        "--plot-config",
        default=str(PLOT_CONFIG_PATH),
        help="Path to plot_panels.json.",
    )
    parser.add_argument(
        "--panel-key",
        default=DEFAULT_PANEL_KEY,
        help="Config key in plot_panels.json.",
    )
    return parser.parse_args()


def _format_filter_pair_label(label: str) -> str:
    parts = str(label).split("|", 1)
    if len(parts) != 2:
        return str(label)
    left, right = [p.strip() for p in parts]
    left_map = {
        "dexseq_0.05": "DexSeq",
        "fisher_0.05": "Fisher",
    }
    left_label = left_map.get(left, left.replace("_", " "))
    if right.startswith("dexseq_log2fc_"):
        right_label = "log2FC > " + right.replace("dexseq_log2fc_", "")
    elif right.startswith("MPRO_"):
        right_label = "MPRO > " + right.replace("MPRO_", "")
    else:
        right_label = right.replace("_", " ")
    return f"{left_label} | {right_label}"


def _format_filter_pair_label_wrapped(label: str) -> str:
    return str(label).replace(" | ", "\n").replace("|", "\n")


def _mm_rect(x_mm: float, y_mm: float, w_mm: float, h_mm: float, width_mm: float, height_mm: float) -> list[float]:
    return [x_mm / width_mm, y_mm / height_mm, w_mm / width_mm, h_mm / height_mm]


def _auto_axis_limits(
    points: pd.DataFrame,
    *,
    y_mean_col: str,
    y_sd_col: str,
    y_clip_max: float | None = 1.02,
    margin_frac: float = 0.08,
    min_span: float = 0.14,
    y_min_span: float | None = None,
) -> tuple[float, float, float, float]:
    x = pd.to_numeric(points["raw_recall_mean"], errors="coerce")
    y = pd.to_numeric(points[y_mean_col], errors="coerce")
    xerr = pd.to_numeric(points["raw_recall_sd"], errors="coerce").fillna(0.0)
    yerr = pd.to_numeric(points[y_sd_col], errors="coerce").fillna(0.0)
    valid = x.notna() & y.notna()
    if int(valid.sum()) == 0:
        return 0.0, 1.02, 0.0, y_clip_max if y_clip_max is not None else 1.02

    def _bound(lo: float, hi: float, clip_max: float | None, span_floor: float) -> tuple[float, float]:
        span = max(hi - lo, span_floor)
        pad = span * margin_frac
        a = max(0.0, lo - pad)
        b = hi + pad if clip_max is None else min(clip_max, hi + pad)
        if b - a < span_floor:
            c = 0.5 * (a + b)
            a = max(0.0, c - span_floor / 2.0)
            b = c + span_floor / 2.0 if clip_max is None else min(clip_max, c + span_floor / 2.0)
        return a, b

    x0, x1 = _bound(float((x[valid] - xerr[valid]).min()), float((x[valid] + xerr[valid]).max()), 1.02, min_span)
    y0, y1 = _bound(
        float((y[valid] - yerr[valid]).min()),
        float((y[valid] + yerr[valid]).max()),
        y_clip_max,
        y_min_span if y_min_span is not None else min_span,
    )
    return x0, x1, y0, y1


def _draw_native_scatter(
    fig: plt.Figure,
    *,
    points_path: Path,
    x0_mm: float,
    y0_mm: float,
    width_mm: float,
    height_mm: float,
    figure_width_mm: float,
    figure_height_mm: float,
    font_size: float,
    marker_size: float,
    capsize: float,
    y_mean_col: str = "sim_f1_mean",
    y_sd_col: str = "sim_f1_sd",
    y_label: str = "Sim F1 (mean ± SD)",
    share_y_axis: bool = False,
    x_limits: tuple[float, float] | None = None,
    y_min_span: float | None = None,
) -> None:
    df = pd.read_csv(points_path, sep="\t")
    required = {"pair_key", "tool", "raw_recall_mean", y_mean_col, "raw_recall_sd", y_sd_col}
    missing = sorted(required - set(df.columns))
    if missing:
        raise KeyError(f"Missing scatter point columns: {missing}")

    pair_keys = df["pair_key"].dropna().astype(str).unique().tolist()
    present_tools = set(df["tool"].astype(str))
    legend_order = [tool for tool in TOOL_ORDER_CP if tool in present_tools]
    legend_order.extend(sorted(tool for tool in present_tools if tool not in legend_order))
    colors = tool_palette(legend_order)

    panel_gap_mm = 1.5
    title_h_mm = 7
    label_h_mm = 7
    legend_w_mm = 18
    y_label_w_mm = 9.5
    plot_area_w_mm = width_mm - y_label_w_mm - legend_w_mm - panel_gap_mm * (len(pair_keys) - 1)
    panel_w_mm = plot_area_w_mm / len(pair_keys)
    panel_h_mm = height_mm - title_h_mm - label_h_mm - 2
    plot_y_mm = y0_mm + label_h_mm
    plot_x0_mm = x0_mm + y_label_w_mm
    shared_limits = None
    if share_y_axis:
        shared_limits = _auto_axis_limits(
            df,
            y_mean_col=y_mean_col,
            y_sd_col=y_sd_col,
            y_clip_max=None,
            y_min_span=y_min_span,
        )

    for idx, pair_key in enumerate(pair_keys):
        panel_df = df[df["pair_key"].astype(str) == pair_key].copy()
        px_mm = plot_x0_mm + idx * (panel_w_mm + panel_gap_mm)
        ax = fig.add_axes(_mm_rect(px_mm, plot_y_mm, panel_w_mm, panel_h_mm, figure_width_mm, figure_height_mm))
        for _, row in panel_df.iterrows():
            tool = str(row["tool"])
            ax.errorbar(
                float(row["raw_recall_mean"]),
                float(row[y_mean_col]),
                xerr=float(row["raw_recall_sd"]),
                yerr=float(row[y_sd_col]),
                fmt="o",
                color=colors.get(tool, "#929292"),
                ecolor=colors.get(tool, "#929292"),
                elinewidth=0.6,
                capsize=capsize,
                markersize=marker_size,
                alpha=0.95,
            )
        xmin, xmax, ymin, ymax = _auto_axis_limits(
            panel_df,
            y_mean_col=y_mean_col,
            y_sd_col=y_sd_col,
            y_clip_max=None if share_y_axis else 1.02,
            y_min_span=y_min_span,
        )
        if x_limits is not None:
            xmin, xmax = x_limits
        ax.set_xlim(xmin, xmax)
        if shared_limits is not None:
            ymin, ymax = shared_limits[2], shared_limits[3]
        ax.set_ylim(ymin, ymax)
        ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.7)
        ax.set_title(_format_filter_pair_label_wrapped(pair_key), fontsize=font_size, pad=2, linespacing=0.9)
        ax.tick_params(axis="both", labelsize=font_size, length=2, pad=2.5)
        if idx == 0:
            ax.set_ylabel(y_label, fontsize=font_size, labelpad=4)
        else:
            ax.set_yticklabels([])
        if idx == len(pair_keys) // 2:
            ax.set_xlabel("Raw Recall (mean ± SD)", fontsize=font_size, labelpad=1)

        x_ser = pd.to_numeric(panel_df["raw_recall_mean"], errors="coerce")
        y_ser = pd.to_numeric(panel_df[y_mean_col], errors="coerce")
        valid = x_ser.notna() & y_ser.notna()
        if int(valid.sum()) >= 2:
            r = x_ser[valid].corr(y_ser[valid], method="spearman")
            if not np.isnan(r):
                ax.text(
                    0.03,
                    0.97,
                    f"Spearman\nr = {r:.2f}",
                    transform=ax.transAxes,
                    fontsize=font_size,
                    va="top",
                    ha="left",
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.7, "pad": 0.8},
                )

    handles = [
        Line2D([0], [0], marker="o", linestyle="", color=colors.get(tool, "#929292"), label=tool, markersize=2.7)
        for tool in legend_order
    ]
    legend_x_mm = plot_x0_mm + plot_area_w_mm + 1.0
    leg_ax = fig.add_axes(_mm_rect(legend_x_mm, plot_y_mm, legend_w_mm, panel_h_mm, figure_width_mm, figure_height_mm))
    leg_ax.axis("off")
    leg_ax.legend(
        handles=handles,
        labels=legend_order,
        loc="center left",
        frameon=False,
        fontsize=font_size,
        handletextpad=0.15,
        borderaxespad=0,
        labelspacing=0.05,
        columnspacing=0.8,
        ncol=1,
    )


def _draw_native_heatmaps(
    fig: plt.Figure,
    *,
    points_path: Path,
    x0_mm: float,
    y0_mm: float,
    width_mm: float,
    height_mm: float,
    figure_width_mm: float,
    figure_height_mm: float,
    font_size: float,
    title_font_size: float,
    colorbar_width_mm: float,
    show_y_labels: bool,
    annotate_cells: bool,
    cmap_name: str,
    vmin: float,
    vmax: float,
    tool_order: list[str],
) -> None:
    df = pd.read_csv(points_path, sep="\t")
    required = {"filter_group", "tool_y", "tool_x", "win_rate"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise KeyError(f"Missing heatmap pairwise columns: {missing}")

    groups = df["filter_group"].dropna().astype(str).unique().tolist()
    if not groups:
        raise ValueError(f"No heatmap groups found in: {points_path}")

    cbar_gap_mm = 5
    cbar_label_w_mm = 9
    title_h_mm = 7
    label_h_mm = 10
    y_label_w_mm = 12.5 if show_y_labels else 1.5
    label_w_mm = 0.5
    inner_gap_mm = 1.5
    cbar_total_mm = colorbar_width_mm + cbar_gap_mm + cbar_label_w_mm
    panel_area_w_mm = width_mm - cbar_total_mm - y_label_w_mm
    cell_area_h_mm = height_mm - title_h_mm - label_h_mm
    cell_size_mm = min((panel_area_w_mm - inner_gap_mm * (len(groups) - 1)) / len(groups), cell_area_h_mm)
    heatmap_block_w_mm = len(groups) * cell_size_mm + (len(groups) - 1) * inner_gap_mm
    used_w_mm = y_label_w_mm + heatmap_block_w_mm + cbar_gap_mm + colorbar_width_mm + cbar_label_w_mm
    left_pad_mm = max(0.0, width_mm - used_w_mm)
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=vmin, vmax=vmax)
    for idx, group in enumerate(groups):
        gdf = df[df["filter_group"].astype(str) == group]
        present = set(gdf["tool_y"].astype(str)) | set(gdf["tool_x"].astype(str))
        tools = [tool for tool in tool_order if tool in present]
        tools.extend(sorted(tool for tool in present if tool not in tools))
        n_tools = len(tools)
        mat = np.full((n_tools, n_tools), np.nan, dtype=float)
        for _, row in gdf.iterrows():
            ty = str(row["tool_y"])
            tx = str(row["tool_x"])
            if ty in tools and tx in tools:
                mat[tools.index(ty), tools.index(tx)] = float(row["win_rate"])
        np.fill_diagonal(mat, np.nan)

        panel_x_mm = x0_mm + left_pad_mm + y_label_w_mm + idx * (cell_size_mm + inner_gap_mm)
        heat_x_mm = panel_x_mm + label_w_mm
        heat_y_mm = y0_mm + label_h_mm
        title_y_mm = heat_y_mm + cell_size_mm + 0.6
        ax = fig.add_axes(_mm_rect(heat_x_mm, heat_y_mm, cell_size_mm, cell_size_mm, figure_width_mm, figure_height_mm))
        ax.set_xlim(-0.5, n_tools - 0.5)
        ax.set_ylim(n_tools - 0.5, -0.5)
        ax.set_facecolor("white")
        for iy in range(n_tools):
            for ix in range(n_tools):
                value = mat[iy, ix]
                facecolor = "white" if np.isnan(value) else cmap(norm(value))
                ax.add_patch(
                    Rectangle(
                        (ix - 0.5, iy - 0.5),
                        1.0,
                        1.0,
                        facecolor=facecolor,
                        edgecolor="white",
                        linewidth=0.5,
                    )
                )
        labels = tools
        ax.set_xticks(range(n_tools))
        ax.set_yticks(range(n_tools))
        ax.set_xticklabels(labels, rotation=90, fontsize=font_size)
        if show_y_labels and idx == 0:
            ax.set_yticklabels(labels, fontsize=font_size)
        else:
            ax.set_yticklabels([])
        ax.tick_params(length=0, pad=2.5)
        for spine in ax.spines.values():
            spine.set_visible(False)
        if annotate_cells:
            for iy in range(n_tools):
                for ix in range(n_tools):
                    if np.isnan(mat[iy, ix]):
                        continue
                    color = "white" if mat[iy, ix] < 0.22 or mat[iy, ix] > 0.78 else "black"
                    ax.text(ix, iy, f"{mat[iy, ix]:.2f}", ha="center", va="center", fontsize=font_size, color=color)

        title_ax = fig.add_axes(_mm_rect(panel_x_mm, title_y_mm, cell_size_mm, title_h_mm, figure_width_mm, figure_height_mm))
        title_ax.axis("off")
        title_ax.text(
            0.5,
            0.0,
            _format_filter_pair_label_wrapped(group),
            ha="center",
            va="bottom",
            fontsize=title_font_size,
            linespacing=0.9,
        )

    cax_x_mm = x0_mm + left_pad_mm + y_label_w_mm + heatmap_block_w_mm + cbar_gap_mm
    cax_h_mm = cell_size_mm * 0.72
    cax_y_mm = y0_mm + label_h_mm + (cell_size_mm - cax_h_mm) / 2
    cax = fig.add_axes(_mm_rect(cax_x_mm, cax_y_mm, colorbar_width_mm, cax_h_mm, figure_width_mm, figure_height_mm))
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.ax.tick_params(labelsize=font_size, length=2, pad=2.5)
    cbar.set_label("Win rate", fontsize=font_size, labelpad=2)


def main() -> None:
    args = parse_args()
    apply_reference_style()
    cfg = json.loads(Path(args.plot_config).read_text(encoding="utf-8"))
    panel = cfg.get(args.panel_key, {})
    if not panel:
        raise KeyError(f"Missing panel config: {args.panel_key}")

    output_pdf = topic_figure_dir(TOPIC) / str(panel.get("output_figure", "raw_recall_main_scatter_heatmap_combined.pdf"))

    width_mm = float(panel.get("width_mm", 178))
    height_mm = float(panel.get("height_mm", 58))
    side_margin_mm = float(panel.get("side_margin_mm", 1))
    vertical_gap_mm = float(panel.get("vertical_gap_mm", 6))
    scatter_height_mm = float(panel.get("scatter_height_mm", 56))
    pd_scatter_height_mm = float(panel.get("pd_scatter_height_mm", scatter_height_mm))
    heatmap_height_mm = float(panel.get("heatmap_height_mm", 56))
    scatter_points = topic_intermediate_dir(TOPIC) / str(panel.get("scatter_points", "raw_recall_vs_sim_f1_scatter_points.tsv"))
    pd_scatter_points = topic_intermediate_dir(TOPIC) / str(
        panel.get("pd_scatter_points", "raw_recall_vs_pd_scatter_points.tsv")
    )
    scatter_font_size = float(panel.get("scatter_font_size", 7))
    scatter_marker_size = float(panel.get("scatter_marker_size", 3))
    scatter_capsize = float(panel.get("scatter_capsize", 1.5))
    pd_scatter_font_size = float(panel.get("pd_scatter_font_size", scatter_font_size))
    pd_scatter_marker_size = float(panel.get("pd_scatter_marker_size", scatter_marker_size))
    pd_scatter_capsize = float(panel.get("pd_scatter_capsize", scatter_capsize))
    pd_scatter_y_min_span = panel.get("pd_scatter_y_min_span")
    pd_scatter_y_min_span = None if pd_scatter_y_min_span is None else float(pd_scatter_y_min_span)
    scatter_x_limits = (
        float(panel.get("scatter_x_min", 0.0)),
        float(panel.get("scatter_x_max", 1.02)),
    )
    heatmap_points = topic_intermediate_dir(TOPIC) / str(panel.get("heatmap_points", "raw_recall_pairwise_winrate.tsv"))
    heatmap_font_size = float(panel.get("heatmap_font_size", 7))
    heatmap_title_font_size = float(panel.get("heatmap_title_font_size", heatmap_font_size))
    heatmap_show_y_labels = bool(panel.get("heatmap_show_y_labels", True))
    heatmap_annotate_cells = bool(panel.get("heatmap_annotate_cells", False))
    heatmap_colorbar_width_mm = float(panel.get("heatmap_colorbar_width_mm", 4))
    heatmap_cfg = cfg.get("heatmap_main_recall", {})
    heatmap_cmap = str(heatmap_cfg.get("cmap", "RdBu"))
    heatmap_vmin = float(heatmap_cfg.get("vmin", 0.0))
    heatmap_vmax = float(heatmap_cfg.get("vmax", 1.0))
    tool_order = [str(t) for t in heatmap_cfg.get("tool_order", TOOL_ORDER_CP)]
    expected_height = scatter_height_mm + vertical_gap_mm + pd_scatter_height_mm + vertical_gap_mm + heatmap_height_mm
    if abs(height_mm - expected_height) > 0.01:
        height_mm = expected_height

    fig = plt.figure(figsize=figsize_from_mm(width_mm, height_mm))
    heatmap_y_mm = 2.0
    pd_scatter_y_mm = heatmap_y_mm + heatmap_height_mm + vertical_gap_mm
    scatter_y_mm = pd_scatter_y_mm + pd_scatter_height_mm + vertical_gap_mm
    _draw_native_scatter(
        fig,
        points_path=scatter_points,
        x0_mm=side_margin_mm,
        y0_mm=scatter_y_mm,
        width_mm=width_mm - 2 * side_margin_mm,
        height_mm=scatter_height_mm,
        figure_width_mm=width_mm,
        figure_height_mm=height_mm,
        font_size=scatter_font_size,
        marker_size=scatter_marker_size,
        capsize=scatter_capsize,
        y_mean_col="sim_f1_mean",
        y_sd_col="sim_f1_sd",
        y_label="Sim F1 (mean ± SD)",
        share_y_axis=False,
        x_limits=scatter_x_limits,
    )

    _draw_native_scatter(
        fig,
        points_path=pd_scatter_points,
        x0_mm=side_margin_mm,
        y0_mm=pd_scatter_y_mm,
        width_mm=width_mm - 2 * side_margin_mm,
        height_mm=pd_scatter_height_mm,
        figure_width_mm=width_mm,
        figure_height_mm=height_mm,
        font_size=pd_scatter_font_size,
        marker_size=pd_scatter_marker_size,
        capsize=pd_scatter_capsize,
        y_mean_col="pd_count_mean",
        y_sd_col="pd_count_sd",
        y_label="Predicted DE APA count (mean ± SD)",
        share_y_axis=True,
        x_limits=scatter_x_limits,
        y_min_span=pd_scatter_y_min_span,
    )

    _draw_native_heatmaps(
        fig,
        points_path=heatmap_points,
        x0_mm=side_margin_mm,
        y0_mm=heatmap_y_mm,
        width_mm=width_mm - 2 * side_margin_mm,
        height_mm=heatmap_height_mm - 2.0,
        figure_width_mm=width_mm,
        figure_height_mm=height_mm,
        font_size=heatmap_font_size,
        title_font_size=heatmap_title_font_size,
        colorbar_width_mm=heatmap_colorbar_width_mm,
        show_y_labels=heatmap_show_y_labels,
        annotate_cells=heatmap_annotate_cells,
        cmap_name=heatmap_cmap,
        vmin=heatmap_vmin,
        vmax=heatmap_vmax,
        tool_order=tool_order,
    )
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_pdf, dpi=300)
    plt.close(fig)
    print(f"Wrote combined main figure: {output_pdf}")
    print(f"Combined size target: {width_mm:.1f}mm x {height_mm:.1f}mm")


if __name__ == "__main__":
    main()
