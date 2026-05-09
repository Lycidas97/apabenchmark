#!/usr/bin/env python3
"""Render pairwise win-rate heatmaps for sim_data_performance."""

from __future__ import annotations

import json
from pathlib import Path
import sys

from matplotlib.cm import ScalarMappable
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.io import read_table
from _shared.paths import topic_figure_dir, topic_intermediate_dir
from _shared.style import apply_reference_style

from _plot_helpers import TOPIC, available_protocol_order
from _plot_config import apply_runtime_rcparams, cfg_bool, cfg_float

MM = 1 / 25.4
INPUT_TABLE = "sim_data_pairwise_winrate.parquet"
INPUT_META = "sim_data_pairwise_winrate_metadata.json"


def _text_bbox_overlaps(texts: list[plt.Text], renderer) -> int:
    bboxes = []
    for t in texts:
        if not t.get_visible():
            continue
        s = str(t.get_text()).strip()
        if not s:
            continue
        bb = t.get_window_extent(renderer=renderer)
        if bb.width <= 0 or bb.height <= 0:
            continue
        bboxes.append(bb)
    overlaps = 0
    for i in range(len(bboxes)):
        for j in range(i + 1, len(bboxes)):
            if bboxes[i].overlaps(bboxes[j]):
                overlaps += 1
    return overlaps


def _title_overlap_count(ax: plt.Axes, renderer) -> int:
    if not ax.title.get_visible() or not str(ax.title.get_text()).strip():
        return 0
    tbb = ax.title.get_window_extent(renderer=renderer)
    others = []
    for t in ax.texts:
        if t is ax.title:
            continue
        others.append(t)
    others.extend(ax.get_xticklabels())
    others.extend(ax.get_yticklabels())
    count = 0
    for t in others:
        if not t.get_visible():
            continue
        s = str(t.get_text()).strip()
        if not s:
            continue
        bb = t.get_window_extent(renderer=renderer)
        if bb.width <= 0 or bb.height <= 0:
            continue
        if tbb.overlaps(bb):
            count += 1
    return count


def _detect_overlap_issues(fig: plt.Figure) -> tuple[int, dict[int, dict[str, int]]]:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    issues: dict[int, dict[str, int]] = {}
    total = 0
    for ax_idx, ax in enumerate(fig.axes):
        annots = [t for t in ax.texts if t.get_gid() == "cell_annot"]
        xticks = [t for t in ax.get_xticklabels() if t.get_visible()]
        yticks = [t for t in ax.get_yticklabels() if t.get_visible()]
        annot_ov = _text_bbox_overlaps(annots, renderer)
        xtick_ov = _text_bbox_overlaps(xticks, renderer)
        ytick_ov = _text_bbox_overlaps(yticks, renderer)
        title_ov = _title_overlap_count(ax, renderer)
        subtotal = annot_ov + xtick_ov + ytick_ov + title_ov
        if subtotal > 0:
            issues[ax_idx] = {
                "annot": annot_ov,
                "xtick": xtick_ov,
                "ytick": ytick_ov,
                "title": title_ov,
            }
            total += subtotal
    if getattr(fig, "_suptitle", None) is not None and fig._suptitle.get_visible():
        suptitle = fig._suptitle
        sbb = suptitle.get_window_extent(renderer=renderer)
        st = str(suptitle.get_text()).strip()
        if st:
            for ax in fig.axes:
                abb = ax.get_window_extent(renderer=renderer)
                if sbb.overlaps(abb):
                    total += 1
                    issues[-1] = {"suptitle": 1}
                    break
    return total, issues


def _reduce_text_sizes(fig: plt.Figure, issues: dict[int, dict[str, int]]) -> None:
    for ax_idx, detail in issues.items():
        if ax_idx == -1:
            if getattr(fig, "_suptitle", None) is not None:
                fs = float(fig._suptitle.get_fontsize())
                fig._suptitle.set_fontsize(max(7.0, fs - 0.4))
            continue
        if ax_idx < 0 or ax_idx >= len(fig.axes):
            continue
        ax = fig.axes[ax_idx]
        if detail.get("annot", 0) > 0:
            for t in ax.texts:
                if t.get_gid() == "cell_annot":
                    fs = float(t.get_fontsize())
                    t.set_fontsize(max(3.2, fs - 0.4))
        if detail.get("xtick", 0) > 0:
            for t in ax.get_xticklabels():
                fs = float(t.get_fontsize())
                t.set_fontsize(max(7.0, fs - 0.3))
        if detail.get("ytick", 0) > 0:
            for t in ax.get_yticklabels():
                fs = float(t.get_fontsize())
                t.set_fontsize(max(7.0, fs - 0.3))
        if detail.get("title", 0) > 0 and ax.title is not None:
            fs = float(ax.title.get_fontsize())
            ax.title.set_fontsize(max(6.5, fs - 0.3))


def _resolve_text_overlaps(fig: plt.Figure, out_path: Path) -> None:
    # Save-time guard: detect and reduce overlaps per figure before writing.
    total, issues = _detect_overlap_issues(fig)
    if total == 0:
        print(f"Checked overlaps: 0 ({out_path.name})")
        return
    for attempt in range(1, 6):
        _reduce_text_sizes(fig, issues)
        total, issues = _detect_overlap_issues(fig)
        if total == 0:
            print(f"Checked overlaps: fixed in {attempt} step(s) ({out_path.name})")
            return
    print(f"Checked overlaps: remaining={total} ({out_path.name})")


def _slug(text: str) -> str:
    out = []
    for ch in text.lower():
        if ch.isalnum():
            out.append(ch)
        else:
            out.append("_")
    s = "".join(out)
    while "__" in s:
        s = s.replace("__", "_")
    return s.strip("_")


def _dataset_label(name: str) -> str:
    mapping = {
        "pas_detect": "PAS Detect",
        "pas_quantification": "PAS Quantification",
        "apa_detect": "APA Detect",
        "apa_detect_single_filter": "APA Detect (Single Filter)",
        "apa_detect_filter_comb": "APA Detect",
    }
    return mapping.get(name, name)


def _metric_label(name: str) -> str:
    mapping = {
        "f1": "F1",
        "precision": "Precision",
        "recall": "Recall",
        "cor_pas": "Correlation",
        "mape_pas": "MAPE (barcode)",
        "mape_pas_ct": "MAPE (group)",
    }
    return mapping.get(name, name)


def _metric_order_for_dataset(dataset: str) -> list[str]:
    if dataset == "pas_quantification":
        return ["mape_pas_ct", "mape_pas", "cor_pas"]
    return ["f1", "precision", "recall"]


def _overall_combined_layout() -> list[dict[str, object]]:
    # Fixed order requested by user:
    # APA Detect -> PAS Identification -> PAS Quantification.
    # Metric order within first two groups: f1 -> recall -> precision.
    # Metric order in PAS Quantification: mape(group) -> mape(barcode) -> correlation.
    layout: list[dict[str, object]] = []
    groups = [
        ("apa_detect", "APA Detect", ["f1", "recall", "precision"]),
        ("pas_detect", "PAS Identification", ["f1", "recall", "precision"]),
        ("pas_quantification", "PAS Quantification", ["mape_pas_ct", "mape_pas", "cor_pas"]),
    ]
    for gi, (dataset, group_label, metrics) in enumerate(groups):
        for metric in metrics:
            layout.append(
                {
                    "is_gap": False,
                    "dataset": dataset,
                    "metric": metric,
                    "group_label": group_label,
                }
            )
        if gi < len(groups) - 1:
            layout.append({"is_gap": True})
    return layout


def _overall_combined_metric_label(name: str) -> str:
    mapping = {
        "mape_pas_ct": "MAPE\n(group)",
        "mape_pas": "MAPE\n(barcode)",
    }
    return mapping.get(name, _metric_label(name))


def _filter_comb_label(text: str) -> str:
    return text.replace("|", "\n")


def _short_filter_label(text: str) -> str:
    parts = [p.strip() for p in text.split("|")]
    if len(parts) != 2:
        return text

    def _term_label(term: str) -> str:
        if "_" in term:
            key, val = term.rsplit("_", 1)
            return f"{key}={val}"
        return term

    left = _term_label(parts[0])
    right = _term_label(parts[1])
    return f"{left} + {right}"


def _filter_color_map(filter_combs: list[str]) -> dict[str, str]:
    colors = ["#4c78a8", "#f58518", "#54a24b", "#b279a2", "#e45756", "#72b7b2"]
    return {name: colors[i % len(colors)] for i, name in enumerate(filter_combs)}


def _load_tool_order(meta_path: Path, df: pd.DataFrame) -> list[str]:
    if meta_path.exists():
        loaded = json.loads(meta_path.read_text(encoding="utf-8"))
        order = loaded.get("tool_order", [])
        if isinstance(order, list) and all(isinstance(x, str) for x in order):
            return order
    observed = set(df["tool_a"].dropna().astype(str)) | set(df["tool_b"].dropna().astype(str))
    return sorted(observed)


def _panel_tool_order(panel_df: pd.DataFrame, global_order: list[str]) -> list[str]:
    present = set(panel_df["tool_a"].astype(str)) | set(panel_df["tool_b"].astype(str))
    ordered = [name for name in global_order if name in present]
    leftovers = sorted(name for name in present if name not in ordered)
    return ordered + leftovers


def _build_matrix(panel_df: pd.DataFrame, tool_order: list[str]) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.full((len(tool_order), len(tool_order)), np.nan, dtype=float)
    n_matrix = np.full((len(tool_order), len(tool_order)), np.nan, dtype=float)
    idx = {name: i for i, name in enumerate(tool_order)}

    for item in panel_df.itertuples(index=False):
        if item.tool_a not in idx or item.tool_b not in idx:
            continue
        i = idx[item.tool_a]
        j = idx[item.tool_b]
        matrix[i, j] = float(item.win_rate)
        n_matrix[i, j] = float(item.n_units)

    for i in range(len(tool_order)):
        matrix[i, i] = np.nan
        n_matrix[i, i] = np.nan
    return matrix, n_matrix


def _annotate_cells(
    ax,
    matrix: np.ndarray,
    n_matrix: np.ndarray,
    *,
    as_percent: bool,
    show_n: bool,
    fontsize: float,
) -> None:
    n_tools = matrix.shape[0]
    for i in range(n_tools):
        for j in range(n_tools):
            if not np.isfinite(matrix[i, j]):
                continue
            value = float(matrix[i, j])
            text = f"{value*100:.0f}%" if as_percent else f"{value:.2f}"
            if show_n and np.isfinite(n_matrix[i, j]):
                text = f"{text}\n(n={int(n_matrix[i, j])})"
            ax.text(j, i, text, ha="center", va="center", fontsize=fontsize, gid="cell_annot")


def _draw_panel(
    ax,
    panel_df: pd.DataFrame,
    tool_order: list[str],
    norm: TwoSlopeNorm,
    cmap: str,
    *,
    as_percent: bool,
    show_n: bool,
    annot_fontsize: float,
    show_xlabels: bool,
    show_ylabels: bool,
    title: str | None = None,
    x_tick_fontsize: float = 7.0,
    y_tick_fontsize: float = 7.0,
    title_fontsize: float = 8.0,
    aspect: str = "equal",
    annotate_cells: bool = False,
) -> None:
    matrix, n_matrix = _build_matrix(panel_df, tool_order)
    cmap_obj = plt.get_cmap(cmap)
    n_tools = len(tool_order)
    for i in range(n_tools):
        for j in range(n_tools):
            value = matrix[i, j]
            facecolor = "#f2f2f2" if not np.isfinite(value) else cmap_obj(norm(float(value)))
            ax.add_patch(
                Rectangle(
                    (j - 0.5, i - 0.5),
                    1.0,
                    1.0,
                    facecolor=facecolor,
                    edgecolor="none",
                    linewidth=0,
                )
            )
    ax.set_xlim(-0.5, n_tools - 0.5)
    ax.set_ylim(n_tools - 0.5, -0.5)
    ax.set_aspect(aspect)
    ax.set_xticks(range(len(tool_order)))
    ax.set_yticks(range(len(tool_order)))
    ax.set_xticklabels(tool_order if show_xlabels else [], rotation=90, fontsize=x_tick_fontsize)
    ax.set_yticklabels(tool_order if show_ylabels else [], fontsize=y_tick_fontsize)
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0)
    if annotate_cells:
        _annotate_cells(
            ax,
            matrix,
            n_matrix,
            as_percent=as_percent,
            show_n=show_n,
            fontsize=annot_fontsize,
        )
    if title:
        ax.set_title(title, fontsize=title_fontsize, pad=4)


def _savefig(fig: plt.Figure, out_path: Path, *, tight: bool = True) -> None:
    _resolve_text_overlaps(fig, out_path)
    bbox = "tight" if tight else None
    fig.savefig(out_path, bbox_inches=bbox, dpi=plt.rcParams.get("savefig.dpi", 300))
    plt.close(fig)
    print(f"Wrote figure: {out_path}")


def main() -> None:
    apply_reference_style()
    apply_runtime_rcparams()

    in_dir = topic_intermediate_dir(TOPIC)
    table_path = in_dir / INPUT_TABLE
    meta_path = in_dir / INPUT_META
    if not table_path.exists():
        raise FileNotFoundError(f"Missing pairwise table: {table_path}. Run 04_prepare_pairwise_winrate.py first.")

    df = read_table(table_path)
    required = {
        "dataset",
        "metric",
        "analysis_type",
        "analysis_value",
        "tool_a",
        "tool_b",
        "win_rate",
        "n_units",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise KeyError(f"Missing required columns in pairwise table: {missing}")
    if df.empty:
        raise ValueError("Pairwise table is empty.")

    # Backward compatibility with earlier outputs.
    if "analysis_protocol" not in df.columns:
        df["analysis_protocol"] = np.where(df["analysis_type"] == "protocol", df["analysis_value"], None)
    if "analysis_filter_comb" not in df.columns:
        df["analysis_filter_comb"] = None

    global_tool_order = _load_tool_order(meta_path, df)
    out_dir = topic_figure_dir(TOPIC) / "heatmap"
    out_dir.mkdir(parents=True, exist_ok=True)
    # Remove legacy split by-filter outputs; they are now merged into stacked panels.
    for stale in out_dir.glob("winrate_*_by_protocol_by_filter.pdf"):
        stale.unlink(missing_ok=True)

    cmap = str(json.loads((SCRIPT_DIR.parent / "config" / "plot_params.json").read_text(encoding="utf-8")).get("pairwise_heatmap", {}).get("cmap", "RdBu_r"))
    width_mm = cfg_float("pairwise_heatmap.figsize_mm.0", 120.0)
    height_mm = cfg_float("pairwise_heatmap.figsize_mm.1", 105.0)
    as_percent = cfg_bool("pairwise_heatmap.annot_as_percent", False)
    show_n = cfg_bool("pairwise_heatmap.annot_show_n", False)
    annot_fontsize = cfg_float("pairwise_heatmap.annot_fontsize", 6.0)
    norm = TwoSlopeNorm(vmin=0.0, vcenter=0.5, vmax=1.0)

    for dataset in sorted(df["dataset"].dropna().astype(str).unique()):
        ds = df[df["dataset"].astype(str) == dataset].copy()
        if ds.empty:
            continue
        metric_order = _metric_order_for_dataset(dataset)
        metric_present = sorted(ds["metric"].dropna().astype(str).unique())
        metrics = [m for m in metric_order if m in metric_present]
        metrics += [m for m in metric_present if m not in metrics]

        # Overall: keep one file per metric.
        for metric in metrics:
            sub = ds[ds["metric"].astype(str) == metric].copy()
            overall = sub[sub["analysis_type"] == "overall"].copy()
            if overall.empty:
                continue
            order = _panel_tool_order(overall, global_tool_order)
            fig, ax = plt.subplots(figsize=(width_mm * MM, height_mm * MM))
            _draw_panel(
                ax,
                overall,
                order,
                norm,
                cmap,
                as_percent=as_percent,
                show_n=show_n,
                annot_fontsize=annot_fontsize,
                show_xlabels=True,
                show_ylabels=True,
                title=f"{_dataset_label(dataset)} | {_metric_label(metric)} | Overall",
            )
            fig.subplots_adjust(left=0.1, right=0.86, bottom=0.1, top=0.95)
            cax = fig.add_axes([0.88, 0.2, 0.02, 0.62])
            cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
            cbar.set_label("Win rate: Row > Column", fontsize=8)
            out_path = out_dir / f"winrate_{_slug(dataset)}_{_slug(metric)}_overall.pdf"
            _savefig(fig, out_path)

        # By protocol: merge metrics into 3 columns, or 9 columns if filter-split exists.
        by_protocol_all = ds[ds["analysis_type"] == "protocol"].copy()
        # Filter-split rows are stored under dataset='apa_detect_filter_comb';
        # map them back when rendering dataset='apa_detect' stacked panel.
        if dataset == "apa_detect":
            by_protocol_filter_all = df[
                (df["dataset"].astype(str) == "apa_detect_filter_comb")
                & (df["analysis_type"].astype(str) == "protocol_filter_comb")
            ].copy()
        else:
            by_protocol_filter_all = ds[ds["analysis_type"] == "protocol_filter_comb"].copy()
        if not by_protocol_all.empty and metrics:
            protocols = available_protocol_order(
                by_protocol_all.rename(columns={"analysis_protocol": "protocol"})
            )
            if protocols:
                metric_cols = [m for m in metrics if not by_protocol_all[by_protocol_all["metric"] == m].empty]
                if metric_cols:
                    use_filter_split = (dataset == "apa_detect") and (not by_protocol_filter_all.empty)
                    plot_source = by_protocol_filter_all if use_filter_split else by_protocol_all
                    order = _panel_tool_order(plot_source, global_tool_order)

                    if use_filter_split:
                        filter_combs = sorted(
                            by_protocol_filter_all["analysis_filter_comb"].dropna().astype(str).unique()
                        )
                        # Keep at most three filter combinations in display.
                        filter_combs = filter_combs[:3]
                        col_specs: list[dict[str, object]] = []
                        for metric_idx, metric in enumerate(metric_cols):
                            for fcomb_idx, fcomb in enumerate(filter_combs):
                                col_specs.append(
                                    {
                                        "metric": metric,
                                        "filter_comb": fcomb,
                                        "is_gap": False,
                                        "metric_idx": metric_idx,
                                        "fcomb_idx": fcomb_idx,
                                    }
                                )
                            if metric_idx < len(metric_cols) - 1:
                                col_specs.append(
                                    {
                                        "metric": metric,
                                        "filter_comb": None,
                                        "is_gap": True,
                                        "metric_idx": metric_idx,
                                        "fcomb_idx": -1,
                                    }
                                )
                    else:
                        fig_w = 182.0 * MM
                        fig_h_mm = min(128.0, max(96.0, 22.0 * len(metric_cols) + 34.0))
                        fig, axes = plt.subplots(
                            len(metric_cols),
                            len(protocols),
                            figsize=(fig_w, fig_h_mm * MM),
                            squeeze=False,
                        )
                        for r, metric in enumerate(metric_cols):
                            for c, protocol in enumerate(protocols):
                                psub = by_protocol_all[
                                    (by_protocol_all["analysis_protocol"].astype(str) == protocol)
                                    & (by_protocol_all["metric"].astype(str) == metric)
                                ]
                                _draw_panel(
                                    axes[r][c],
                                    psub,
                                    order,
                                    norm,
                                    cmap,
                                    as_percent=as_percent,
                                    show_n=show_n,
                                    annot_fontsize=max(4.8, min(6.0, annot_fontsize)),
                                    show_xlabels=(r == len(metric_cols) - 1),
                                    show_ylabels=(c == 0),
                                    title=protocol if r == 0 else None,
                                    x_tick_fontsize=7.2,
                                    y_tick_fontsize=7.2,
                                    title_fontsize=7.8,
                                    aspect="equal",
                                    annotate_cells=False,
                                )
                        fig.suptitle(f"{_dataset_label(dataset)} | By Protocol", fontsize=9, y=0.985)
                        fig.subplots_adjust(
                            left=0.135,
                            right=0.895,
                            bottom=0.14,
                            top=0.90,
                            wspace=0.045,
                            hspace=0.10,
                        )
                        for r, metric in enumerate(metric_cols):
                            bbox = axes[r][0].get_position()
                            fig.text(
                                0.025,
                                (bbox.y0 + bbox.y1) / 2.0,
                                _metric_label(metric),
                                rotation=90,
                                ha="center",
                                va="center",
                                fontsize=8.0,
                            )
                        cax = fig.add_axes([0.91, 0.22, 0.012, 0.56])
                        cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
                        cbar.set_label("Win rate: Row > Column", fontsize=8)
                        out_path = out_dir / f"winrate_{_slug(dataset)}_by_protocol_stacked.pdf"
                        _savefig(fig, out_path, tight=False)
                        continue

                    fig_w = 185.0 * MM
                    # Keep size within limits after accounting for labels/colorbar.
                    fig_w = 182.0 * MM
                    fig_h_mm = min(198.0, max(158.0, 20.0 * len(protocols) + 14.0))
                    fig_h = fig_h_mm * MM
                    width_ratios = [
                        0.22 if bool(spec.get("is_gap", False)) else 1.0
                        for spec in col_specs
                    ]
                    fig, axes = plt.subplots(
                        len(protocols),
                        len(col_specs),
                        figsize=(fig_w, fig_h),
                        squeeze=False,
                        gridspec_kw={"width_ratios": width_ratios},
                    )
                    for r, protocol in enumerate(protocols):
                        filter_colors = _filter_color_map(filter_combs) if use_filter_split else {}
                        group_size = len(filter_combs) if use_filter_split else 0
                        for c, spec in enumerate(col_specs):
                            if bool(spec.get("is_gap", False)):
                                axes[r][c].axis("off")
                                continue
                            metric = str(spec["metric"])
                            filter_comb = spec["filter_comb"]
                            if filter_comb is None:
                                psub = by_protocol_all[
                                    (by_protocol_all["analysis_protocol"].astype(str) == protocol)
                                    & (by_protocol_all["metric"].astype(str) == metric)
                                ]
                            else:
                                psub = by_protocol_filter_all[
                                    (by_protocol_filter_all["analysis_protocol"].astype(str) == protocol)
                                    & (by_protocol_filter_all["metric"].astype(str) == metric)
                                    & (by_protocol_filter_all["analysis_filter_comb"].astype(str) == filter_comb)
                                ]
                            annot_fs = max(4.8, min(6.0, annot_fontsize))
                            if r == 0:
                                if filter_comb is None:
                                    col_title = _metric_label(metric)
                                else:
                                    if int(spec.get("fcomb_idx", -1)) == (group_size // 2):
                                        col_title = _metric_label(metric)
                                    else:
                                        col_title = None
                            else:
                                col_title = None
                            _draw_panel(
                                axes[r][c],
                                psub,
                                order,
                                norm,
                                cmap,
                                as_percent=as_percent,
                                show_n=show_n,
                                annot_fontsize=annot_fs,
                                show_xlabels=(r == len(protocols) - 1),
                                show_ylabels=(c == 0),
                                title=col_title,
                                x_tick_fontsize=8.0,
                                y_tick_fontsize=8.0,
                                title_fontsize=8.0,
                                aspect="equal",
                                annotate_cells=False,
                            )
                            if r == 0 and filter_comb is not None:
                                axes[r][c].add_patch(
                                    Rectangle(
                                        (0.02, 1.01),
                                        0.96,
                                        0.06,
                                        transform=axes[r][c].transAxes,
                                        facecolor=filter_colors.get(filter_comb, "#999999"),
                                        edgecolor="none",
                                        clip_on=False,
                                    )
                                )
                            if c == 0:
                                # Explicit per-row protocol tag at left margin.
                                axes[r][c].set_ylabel(protocol, fontsize=8.0, labelpad=10.0)
                    if use_filter_split and filter_combs:
                        handles = [
                            Patch(
                                facecolor=filter_colors[name],
                                edgecolor="none",
                                label=_short_filter_label(name),
                            )
                            for name in filter_combs
                        ]
                        fig.legend(
                            handles=handles,
                            loc="upper center",
                            bbox_to_anchor=(0.5, 0.962),
                            ncol=len(handles),
                            frameon=False,
                            fontsize=7,
                            handlelength=1.1,
                            columnspacing=0.9,
                        )
                    fig.suptitle(f"{_dataset_label(dataset)} | By Protocol", fontsize=9, y=0.985)
                    if use_filter_split:
                        fig.subplots_adjust(left=0.115, right=0.88, bottom=0.11, top=0.91, wspace=0.008, hspace=0.018)
                    else:
                        fig.subplots_adjust(left=0.12, right=0.88, bottom=0.07, top=0.92, wspace=0.02, hspace=0.02)
                    cax = fig.add_axes([0.895, 0.2, 0.012, 0.62])
                    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
                    cbar.set_label("Win rate: Row > Column", fontsize=8)
                    out_path = out_dir / f"winrate_{_slug(dataset)}_by_protocol_stacked.pdf"
                    _savefig(fig, out_path, tight=False)

    # Combined overall: merge 9 overall panels into one horizontal row with shared colorbar.
    overall_specs = _overall_combined_layout()
    width_ratios = [0.10 if bool(spec.get("is_gap", False)) else 1.0 for spec in overall_specs]
    fig_w_mm = 178.0
    fig_h_mm = 68.0
    fig, axes = plt.subplots(
        1,
        len(overall_specs),
        figsize=(fig_w_mm * MM, fig_h_mm * MM),
        squeeze=False,
        gridspec_kw={"width_ratios": width_ratios},
    )
    showed_ylabels = False
    first_panel_in_group = {"apa_detect": True, "pas_detect": True, "pas_quantification": True}
    for c, spec in enumerate(overall_specs):
        ax = axes[0][c]
        if bool(spec.get("is_gap", False)):
            ax.axis("off")
            continue
        dataset = str(spec["dataset"])
        metric = str(spec["metric"])
        group_label = str(spec["group_label"])
        sub = df[
            (df["dataset"].astype(str) == dataset)
            & (df["metric"].astype(str) == metric)
            & (df["analysis_type"].astype(str) == "overall")
        ].copy()
        if sub.empty:
            ax.axis("off")
            continue
        order = _panel_tool_order(sub, global_tool_order)
        metric_title = _overall_combined_metric_label(metric)
        title = f"{group_label}\n{metric_title}" if first_panel_in_group.get(dataset, False) else metric_title
        _draw_panel(
            ax,
            sub,
            order,
            norm,
            cmap,
            as_percent=as_percent,
            show_n=show_n,
            annot_fontsize=annot_fontsize,
            show_xlabels=True,
            show_ylabels=(not showed_ylabels),
            title=title,
            x_tick_fontsize=7.0,
            y_tick_fontsize=7.0,
            title_fontsize=7.5,
            aspect="equal",
            annotate_cells=False,
        )
        if not showed_ylabels:
            showed_ylabels = True
        first_panel_in_group[dataset] = False
    fig.subplots_adjust(left=0.08, right=0.90, bottom=0.29, top=0.80, wspace=0.012, hspace=0.0)
    cax = fig.add_axes([0.912, 0.29, 0.010, 0.48])
    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label("Win rate: Row > Column", fontsize=7, labelpad=2)
    out_path = out_dir / "winrate_overall_combined_9panel.pdf"
    _savefig(fig, out_path, tight=False)


if __name__ == "__main__":
    main()
