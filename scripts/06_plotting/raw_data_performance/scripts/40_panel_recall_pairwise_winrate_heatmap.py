#!/usr/bin/env python3
"""Render pairwise win-rate heatmap for recall using Marsilea."""

from __future__ import annotations

import argparse
import json
import re
import tempfile
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
import marsilea as ma
import marsilea.plotter as mp

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import TOOL_ORDER_CP
from _shared.io import read_table, write_tsv
from _shared.paths import topic_figure_dir, topic_intermediate_dir, topic_root
from _shared.style import apply_reference_style, figsize_from_mm

TOPIC = "raw_data_performance"
PLOT_CONFIG_PATH = topic_root(TOPIC) / "config" / "plot_panels.json"
DEFAULT_PANEL_KEY = "heatmap_main_recall"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render pairwise recall win-rate heatmap.")
    parser.add_argument(
        "--plot-config",
        default=str(PLOT_CONFIG_PATH),
        help="Path to plot_panels.json.",
    )
    parser.add_argument(
        "--panel-key",
        default=DEFAULT_PANEL_KEY,
        help="Config key in plot_panels.json (e.g. heatmap_main_recall).",
    )
    return parser.parse_args()


def _with_unit_id(df: pd.DataFrame) -> pd.DataFrame:
    unit_cols = [
        col
        for col in ["sample_id", "pair_label", "filter_type_1", "filter_type_2"]
        if col in df.columns
    ]
    if not unit_cols:
        raise KeyError("Cannot build unit_id: none of sample_id/pair_label/filter_type_1/filter_type_2 found.")
    out = df.copy()
    out["unit_id"] = out[unit_cols].astype(str).agg("|".join, axis=1)
    return out


def _species_filter(df: pd.DataFrame, analysis_set: str) -> pd.DataFrame:
    mode = str(analysis_set).strip().lower()
    if mode in {"", "all"}:
        return df
    if "sample" not in df.columns:
        return df.iloc[0:0].copy()
    mask = df["sample"].astype(str).str.lower().str.contains(mode, na=False)
    return df[mask].copy()


def _pairwise_rate(a: pd.Series, b: pd.Series, epsilon: float) -> tuple[float, int]:
    valid = a.notna() & b.notna()
    n = int(valid.sum())
    if n == 0:
        return np.nan, 0
    delta = a[valid] - b[valid]
    wins = np.where(delta > epsilon, 1.0, np.where(np.abs(delta) <= epsilon, 0.5, 0.0))
    return float(np.mean(wins)), n


def _pairwise_rate_balanced(
    a: pd.Series,
    b: pd.Series,
    stratum: pd.Series,
    epsilon: float,
) -> tuple[float, int]:
    valid = a.notna() & b.notna() & stratum.notna()
    n_total = int(valid.sum())
    if n_total == 0:
        return np.nan, 0

    rates: list[float] = []
    frame = pd.DataFrame({"a": a[valid], "b": b[valid], "stratum": stratum[valid]})
    for _, grp in frame.groupby("stratum", dropna=True):
        delta = grp["a"] - grp["b"]
        wins = np.where(delta > epsilon, 1.0, np.where(np.abs(delta) <= epsilon, 0.5, 0.0))
        if len(wins) > 0:
            rates.append(float(np.mean(wins)))

    if not rates:
        return np.nan, 0
    return float(np.mean(rates)), n_total


def _slugify(text: str) -> str:
    s = re.sub(r"[^0-9A-Za-z]+", "_", str(text).strip())
    s = s.strip("_")
    return s or "group"


def _with_suffix(file_name: str, suffix: str) -> str:
    p = Path(file_name)
    return f"{p.stem}_{suffix}{p.suffix}"


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


def _render_marsilea_panel(
    win: np.ndarray,
    tool_order: list[str],
    *,
    cmap: str,
    center: float,
    vmin: float,
    vmax: float,
    annot_fmt: str,
    lower_triangle_only: bool,
    width: float,
    height: float,
    title: str,
    show_n: bool,
    nmat: np.ndarray,
    output_path: Path,
    show_legend: bool,
    font_size: float,
) -> None:
    n_tools = len(tool_order)
    mask = np.zeros_like(win, dtype=bool)
    if lower_triangle_only:
        mask[np.triu_indices(n_tools, k=0)] = True
    else:
        np.fill_diagonal(mask, True)

    heat = ma.Heatmap(
        win,
        cmap=cmap,
        center=center,
        vmin=vmin,
        vmax=vmax,
        mask=mask,
        annot=win,
        fmt=annot_fmt,
        annot_kws={"fontsize": font_size},
        linewidth=0.5,
        linecolor="white",
        width=width,
        height=height,
    )
    heat.add_left(mp.Labels(tool_order, fontsize=font_size))
    heat.add_bottom(mp.Labels(tool_order, fontsize=font_size))
    if show_legend:
        heat.add_legends()
    heat.render()
    ax = heat.get_main_ax()
    ax.set_title(title, loc="left", fontsize=max(font_size, 7))
    if show_n:
        for iy in range(n_tools):
            for ix in range(n_tools):
                if mask[iy, ix] or np.isnan(win[iy, ix]):
                    continue
                ax.text(ix + 0.5, iy + 0.72, f"n={nmat[iy, ix]}", ha="center", va="center", fontsize=font_size)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    heat.save(output_path)


def main() -> None:
    args = parse_args()
    apply_reference_style()

    cfg = json.loads(Path(args.plot_config).read_text(encoding="utf-8"))
    panel = cfg.get(args.panel_key, {})
    if not panel:
        raise KeyError(f"Missing panel config: {args.panel_key}")

    input_name = str(panel.get("input_name", "raw_data_recall_gt_filtered.parquet"))
    metric_col = str(panel.get("metric_col", "recall"))
    output_figure = str(panel.get("output_figure", "raw_recall_pairwise_winrate_heatmap.pdf"))
    output_points = str(panel.get("output_points", "raw_recall_pairwise_winrate.tsv"))
    figsize_mm = panel.get("figsize_mm", [98, 92])
    combined_figsize_mm = panel.get("combined_figsize_mm", [178, 62])
    epsilon = float(panel.get("epsilon", 0.0))
    lower_triangle_only = bool(panel.get("lower_triangle_only", True))
    show_n = bool(panel.get("show_n", False))
    balanced_by_stratum = bool(panel.get("balanced_by_stratum", False))
    stratum_col = str(panel.get("stratum_col", "filter_pair"))
    annot_fmt = str(panel.get("annot_fmt", ".2f"))
    font_size = float(panel.get("font_size", 7))
    cmap = str(panel.get("cmap", "RdBu_r"))
    center = float(panel.get("center", 0.5))
    vmin = float(panel.get("vmin", 0.0))
    vmax = float(panel.get("vmax", 1.0))
    title = str(panel.get("title", "Pairwise Win Rate"))
    analysis_set = str(panel.get("analysis_set", "all"))
    split_by_filter_pair = bool(panel.get("split_by_filter_pair", True))
    filter_group_col = str(panel.get("filter_group_col", "filter_pair"))
    combine_groups_into_one_figure = bool(panel.get("combine_groups_into_one_figure", args.panel_key == "heatmap_main_recall"))

    tool_order_cfg = panel.get("tool_order", TOOL_ORDER_CP)
    tool_order_global = [str(t) for t in tool_order_cfg]

    in_path = topic_intermediate_dir(TOPIC) / input_name
    if not in_path.exists():
        raise FileNotFoundError(f"Missing input table: {in_path}")

    df = read_table(in_path)
    required = {"tool", metric_col}
    missing = sorted(required - set(df.columns))
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    df = _with_unit_id(df)
    df = _species_filter(df, analysis_set)
    df[metric_col] = pd.to_numeric(df[metric_col], errors="coerce")
    df = df[df[metric_col].notna()].copy()
    if df.empty:
        raise ValueError("No valid rows after filtering input for heatmap.")

    if split_by_filter_pair and filter_group_col in df.columns:
        groups = [(str(k), g.copy()) for k, g in df.groupby(filter_group_col, dropna=False)]
    else:
        groups = [("all", df.copy())]

    if isinstance(figsize_mm, list) and len(figsize_mm) == 2:
        width, height = figsize_from_mm(float(figsize_mm[0]), float(figsize_mm[1]))
    else:
        width, height = figsize_from_mm(98, 92)

    rows: list[dict] = []
    panel_payloads: list[dict] = []
    for group_name, gdf in groups:
        agg = gdf.groupby(["unit_id", "tool"], as_index=False)[metric_col].mean()
        pivot = agg.pivot(index="unit_id", columns="tool", values=metric_col)

        present = set(pivot.columns.astype(str).tolist())
        tool_order = [t for t in tool_order_global if t in present]
        tool_order.extend(sorted(t for t in present if t not in tool_order))
        if len(tool_order) < 2:
            print(f"Skip group {group_name}: fewer than 2 tools present.")
            continue

        if balanced_by_stratum and stratum_col in gdf.columns:
            unit_stratum = (
                gdf[["unit_id", stratum_col]]
                .dropna()
                .drop_duplicates(subset=["unit_id"])
                .set_index("unit_id")[stratum_col]
                .reindex(pivot.index)
            )
        else:
            unit_stratum = pd.Series(index=pivot.index, dtype="object")

        n_tools = len(tool_order)
        win = np.full((n_tools, n_tools), np.nan, dtype=float)
        nmat = np.zeros((n_tools, n_tools), dtype=int)

        # Cell(y, x): win-rate that y-axis tool beats x-axis tool.
        for iy, ty in enumerate(tool_order):
            col_y = pivot[ty] if ty in pivot.columns else pd.Series(index=pivot.index, dtype=float)
            for ix, tx in enumerate(tool_order):
                if iy == ix:
                    continue
                col_x = pivot[tx] if tx in pivot.columns else pd.Series(index=pivot.index, dtype=float)
                # Cell(y, x): win-rate after axis-wise comparison inversion.
                if balanced_by_stratum and stratum_col in gdf.columns:
                    rate, n = _pairwise_rate_balanced(col_y, col_x, unit_stratum, epsilon)
                else:
                    rate, n = _pairwise_rate(col_y, col_x, epsilon)
                win[iy, ix] = rate
                nmat[iy, ix] = n
                rows.append(
                    {
                        "filter_group": group_name,
                        "tool_y": ty,
                        "tool_x": tx,
                        "win_rate": rate,
                        "n": n,
                        "epsilon": epsilon,
                        "balanced_by_stratum": balanced_by_stratum,
                        "analysis_set": analysis_set,
                        "panel_key": args.panel_key,
                    }
                )
        panel_payloads.append(
            {
                "group_name": group_name,
                "win": win,
                "nmat": nmat,
                "tool_order": tool_order,
            }
        )

    if combine_groups_into_one_figure and len(panel_payloads) > 1:
        with tempfile.TemporaryDirectory(prefix="marsilea_heatmap_") as td:
            td_path = Path(td)
            panel_pngs: list[tuple[str, Path]] = []
            for p in panel_payloads:
                png_path = td_path / f"{_slugify(str(p['group_name']))}.png"
                _render_marsilea_panel(
                    p["win"],
                    p["tool_order"],
                    cmap=cmap,
                    center=center,
                    vmin=vmin,
                    vmax=vmax,
                    annot_fmt=annot_fmt,
                    lower_triangle_only=lower_triangle_only,
                    width=width,
                    height=height,
                    title=_format_filter_pair_label(str(p["group_name"])),
                    show_n=show_n,
                    nmat=p["nmat"],
                    output_path=png_path,
                    show_legend=False,
                    font_size=font_size,
                )
                panel_pngs.append((str(p["group_name"]), png_path))

            n = len(panel_pngs)
            if isinstance(combined_figsize_mm, list) and len(combined_figsize_mm) == 2:
                combined_size = figsize_from_mm(float(combined_figsize_mm[0]), float(combined_figsize_mm[1]))
            else:
                combined_size = figsize_from_mm(178, 62)
            fig = plt.figure(figsize=combined_size)
            gs = fig.add_gridspec(1, n + 1, width_ratios=[1] * n + [0.045], wspace=0.03)
            axes = [fig.add_subplot(gs[0, i]) for i in range(n)]
            cax_host = fig.add_subplot(gs[0, n])
            cax_host.axis("off")

            for ax, (gname, path) in zip(axes, panel_pngs):
                img = mpimg.imread(path)
                ax.imshow(img)
                ax.set_title(_format_filter_pair_label(gname), fontsize=max(font_size, 7))
                ax.axis("off")

            sm = plt.cm.ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=plt.get_cmap(cmap))
            sm.set_array([])
            cax = cax_host.inset_axes([0.18, 0.26, 0.64, 0.48])
            cbar = fig.colorbar(sm, cax=cax)
            cbar.ax.tick_params(labelsize=font_size)
            cbar.set_label("Win rate", fontsize=font_size)
            fig.suptitle(title, fontsize=max(font_size + 1, 8), y=0.98)
            fig.subplots_adjust(left=0.02, right=0.98, top=0.9, bottom=0.05)

            out_fig = topic_figure_dir(TOPIC) / output_figure
            out_fig.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out_fig, dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"Wrote combined figure: {out_fig}")
    else:
        for p in panel_payloads:
            group_name = str(p["group_name"])
            if len(panel_payloads) == 1:
                out_fig = topic_figure_dir(TOPIC) / output_figure
            else:
                fig_name = _with_suffix(output_figure, _slugify(group_name))
                out_fig = topic_figure_dir(TOPIC) / fig_name
            _render_marsilea_panel(
                p["win"],
                p["tool_order"],
                cmap=cmap,
                center=center,
                vmin=vmin,
                vmax=vmax,
                annot_fmt=annot_fmt,
                lower_triangle_only=lower_triangle_only,
                width=width,
                height=height,
                title=f"{title} | {_format_filter_pair_label(group_name)}",
                show_n=show_n,
                nmat=p["nmat"],
                output_path=out_fig,
                show_legend=True,
                font_size=font_size,
            )
            print(f"Wrote figure: {out_fig}")

    out_points = topic_intermediate_dir(TOPIC) / output_points
    write_tsv(pd.DataFrame(rows), out_points)
    print(f"Wrote pairwise table: {out_points}")


if __name__ == "__main__":
    main()
