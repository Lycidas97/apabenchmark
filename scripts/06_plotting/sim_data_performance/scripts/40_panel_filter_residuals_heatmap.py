#!/usr/bin/env python3
"""Render APA filter residual heatmap for sim data performance (reference-style marsilea)."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import marsilea as ma
import marsilea.plotter as mp

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.paths import topic_figure_dir
from _shared.style import apply_reference_style

from _plot_helpers import TOPIC, load_apa_residual_by_filter_table
from _plot_config import apply_runtime_rcparams, cfg, figsize_for_figure

OUTPUT_NAME = "filter_residuals_heatmap.pdf"
MEAN_TRACK_COLOR = "#386b98"
MM = 1 / 25.4

FILTER_TYPE_1_MAP = {
    "dexseq_0.05": "DexSeq",
    "fisher_0.05": "Fisher",
    "wilcox_PDUI_0.05": "Wilcox PDUI",
    "wilcox_PPUI_0.05": "Wilcox PPUI",
    "wilcox_RWUI_0.05": "Wilcox RWUI",
    "wilcox_DWUI_0.05": "Wilcox DWUI",
}

FILTER_TYPE_1_ORDER = [
    "dexseq_0.05",
    "fisher_0.05",
    "wilcox_PDUI_0.05",
    "wilcox_PPUI_0.05",
    "wilcox_RWUI_0.05",
    "wilcox_DWUI_0.05",
]

FILTER_TYPE_2_ORDER = [
    "PDUI_0.1",
    "PDUI_0.2",
    "PDUI_0.3",
    "PDUI_0.4",
    "PDUI_0.5",
    "PPUI_0.1",
    "PPUI_0.2",
    "PPUI_0.3",
    "PPUI_0.4",
    "PPUI_0.5",
    "RWUI_0.1",
    "RWUI_0.2",
    "RWUI_0.3",
    "RWUI_0.4",
    "RWUI_0.5",
    "DWUI_0.1",
    "DWUI_0.2",
    "DWUI_0.3",
    "DWUI_0.4",
    "DWUI_0.5",
    "MPRO_0.1",
    "MPRO_0.2",
    "MPRO_0.3",
    "MPRO_0.4",
    "MPRO_0.5",
    "dexseq_log2fc_0.25",
    "dexseq_log2fc_0.5",
    "dexseq_log2fc_0.75",
    "dexseq_log2fc_1",
    "dexseq_log2fc_1.25",
]


def format_filter_type_2(name: str) -> str:
    if name.startswith("dexseq_log2fc"):
        return name.replace("dexseq_log2fc_", "DexSeq log2FC ")
    if name.startswith("MPRO"):
        return name.replace("_", " > ")
    return name.replace("_", " diff > ")


def main() -> None:
    apply_reference_style()
    apply_runtime_rcparams()
    # Match reference notebook heatmap tick settings.
    plt.rcParams["xtick.direction"] = "out"
    plt.rcParams["ytick.direction"] = "out"
    plt.rcParams["xtick.bottom"] = False
    plt.rcParams["xtick.minor.bottom"] = True
    plt.rcParams["ytick.left"] = False

    df = load_apa_residual_by_filter_table()

    required = {"filter_type_1", "filter_type_2", "f1_residuals"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise KeyError(f"Missing required columns for heatmap: {missing}")
    if df.empty:
        raise ValueError("Residual-by-filter table is empty.")

    df = df.copy()
    df["filter_type_1"] = (
        df["filter_type_1"].astype("string").str.strip().replace({"": pd.NA, "nan": pd.NA, "None": pd.NA, "<NA>": pd.NA})
    )
    df["filter_type_2"] = (
        df["filter_type_2"].astype("string").str.strip().replace({"": pd.NA, "nan": pd.NA, "None": pd.NA, "<NA>": pd.NA})
    )
    df = df.dropna(subset=["filter_type_1", "filter_type_2", "f1_residuals"]).copy()
    if df.empty:
        raise ValueError("Residual-by-filter table has no valid filter pairs after NaN cleanup.")

    observed_f1 = set(df["filter_type_1"].astype(str))
    observed_f2 = set(df["filter_type_2"].astype(str))
    f1_order = [name for name in FILTER_TYPE_1_ORDER if name in observed_f1]
    f2_order = [name for name in FILTER_TYPE_2_ORDER if name in observed_f2]
    if not f1_order:
        f1_order = sorted(df["filter_type_1"].dropna().astype(str).unique())
    if not f2_order:
        f2_order = sorted(df["filter_type_2"].dropna().astype(str).unique())

    df["filter_type_1"] = pd.Categorical(df["filter_type_1"], categories=f1_order, ordered=True)
    df["filter_type_2"] = pd.Categorical(df["filter_type_2"], categories=f2_order, ordered=True)
    df = df.sort_values(["filter_type_1", "filter_type_2"])

    pivot = df.pivot(index="filter_type_1", columns="filter_type_2", values="f1_residuals")
    pivot = pivot.dropna(axis=0, how="all").dropna(axis=1, how="all")
    if pivot.empty:
        raise ValueError("Heatmap pivot is empty after removing all-NaN rows/columns.")

    ordered_cols = [name for name in FILTER_TYPE_2_ORDER if name in set(pivot.columns.astype(str))]
    ordered_cols += [name for name in pivot.columns.astype(str) if name not in ordered_cols]
    pivot = pivot.reindex(columns=ordered_cols)

    row_mean = pivot.mean(axis=1).fillna(0.0)
    col_mean = pivot.mean(axis=0).fillna(0.0)

    row_labels = [FILTER_TYPE_1_MAP.get(name, name) for name in pivot.index.astype(str)]
    col_labels = [format_filter_type_2(name) for name in ordered_cols]

    # Keep backward compatibility with legacy key while preferring per-figure key.
    width_in, height_in = figsize_for_figure(OUTPUT_NAME, "heatmap_wh_mm")
    width_mm = width_in / MM
    height_mm = height_in / MM
    legacy_heatmap_size = cfg("size_presets_mm.heatmap_wh_mm", None)
    if isinstance(legacy_heatmap_size, (list, tuple)) and len(legacy_heatmap_size) == 2:
        try:
            width_mm = float(legacy_heatmap_size[0])
            height_mm = float(legacy_heatmap_size[1])
        except (TypeError, ValueError):
            pass

    h = ma.Heatmap(
        pivot,
        cmap="RdBu",
        label="Mean Residual",
        width=width_mm * MM,
        height=height_mm * MM,
        cbar_kws={"width": 1, "height": 10},
    )

    h.add_right(
        mp.Numbers(
            row_mean,
            color=MEAN_TRACK_COLOR,
            label="Mean Residual",
            show_value=False,
            props={"fontsize": 7},
        ),
        size=8 * MM,
        pad=0.5 * MM,
    )
    h.add_top(
        mp.Numbers(
            col_mean,
            color=MEAN_TRACK_COLOR,
            label="Mean DE APA F1",
            show_value=False,
            props={"fontsize": 7},
        ),
        size=10 * MM,
        pad=1 * MM,
        name="top",
    )

    var_num = 5
    h.cut_cols(list(range(var_num, len(ordered_cols), var_num)), spacing=0.01)
    h.add_left(mp.Labels(row_labels, rotation=0, fontsize=8), size=7 * MM, pad=0.5 * MM)
    h.add_bottom(mp.Labels(col_labels, rotation=90, fontsize=8), size=10 * MM, pad=0.5 * MM)
    h.add_legends("right")
    h.render()

    out_path = topic_figure_dir(TOPIC) / OUTPUT_NAME
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, bbox_inches="tight", dpi=plt.rcParams.get("savefig.dpi", 300))
    plt.close("all")
    print(f"Wrote figure: {out_path}")


if __name__ == "__main__":
    main()
