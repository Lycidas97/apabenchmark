#!/usr/bin/env python3
"""Render F1 vs stringency linear-regression demo panel."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.paths import topic_figure_dir
from _shared.style import apply_reference_style, save_figure

from _plot_helpers import TOPIC, load_apa_residual_table
from _plot_config import apply_runtime_rcparams, cfg, cfg_float, cfg_int, figsize_for_figure

OUTPUT_NAME = "f1_vs_stringency_linear_regression.pdf"


def _fit_line(x: pd.Series, y: pd.Series) -> tuple[float, float, float]:
    x_num = pd.to_numeric(x, errors="coerce")
    y_num = pd.to_numeric(y, errors="coerce")
    valid = x_num.notna() & y_num.notna()
    if int(valid.sum()) < 2:
        return np.nan, np.nan, np.nan

    xv = x_num.loc[valid].to_numpy(dtype=float)
    yv = y_num.loc[valid].to_numpy(dtype=float)
    design = np.column_stack((np.ones(len(xv)), xv))
    coef = np.linalg.lstsq(design, yv, rcond=None)[0]
    y_hat = design @ coef
    ss_res = float(np.sum((yv - y_hat) ** 2))
    ss_tot = float(np.sum((yv - np.mean(yv)) ** 2))
    r2 = np.nan if ss_tot == 0 else 1.0 - (ss_res / ss_tot)
    intercept = float(coef[0])
    slope = float(coef[1])
    return intercept, slope, r2


def main() -> None:
    apply_reference_style()
    apply_runtime_rcparams()
    df = load_apa_residual_table()

    required = {"filter_type_comb", "gt_recall", "f1"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise KeyError(f"Missing required columns for regression panel: {missing}")

    plot_df = df[list(required)].copy()
    plot_df["gt_recall"] = pd.to_numeric(plot_df["gt_recall"], errors="coerce")
    plot_df["f1"] = pd.to_numeric(plot_df["f1"], errors="coerce")
    plot_df["filter_type_comb"] = plot_df["filter_type_comb"].astype("string")
    plot_df = plot_df.dropna(subset=["filter_type_comb", "gt_recall", "f1"]).copy()
    if plot_df.empty:
        raise ValueError("Residual table has no valid rows for F1/stringency regression.")

    # Match residual scatter density: one averaged point per filter combination.
    agg_df = (
        plot_df.groupby("filter_type_comb", as_index=False)
        .agg(
            gt_recall=("gt_recall", "mean"),
            f1=("f1", "mean"),
            n=("f1", "size"),
        )
        .sort_values("filter_type_comb")
    )
    if agg_df.empty:
        raise ValueError("No valid averaged points after filter-combination aggregation.")

    fig, ax = plt.subplots(
        figsize=figsize_for_figure(OUTPUT_NAME, "f1_stringency_regression"),
        dpi=cfg_int("render.per_figure_dpi", 300),
    )
    sns.scatterplot(
        data=agg_df,
        x="gt_recall",
        y="f1",
        s=14,
        alpha=0.55,
        color="#454545",
        edgecolor=None,
        ax=ax,
        legend=False,
    )
    sns.regplot(
        data=agg_df,
        x="gt_recall",
        y="f1",
        scatter=False,
        ci=None,
        line_kws={"color": "#d34123", "linewidth": 1.0},
        ax=ax,
    )
    intercept, slope, r2 = _fit_line(agg_df["gt_recall"], agg_df["f1"])
    eq = (
        f"F1 = {intercept:.3f} + {slope:.3f} * stringency\n"
        f"$R^2$ = {r2:.3f}, n = {len(agg_df)} filter combinations"
    )
    ax.text(0.03, 0.97, eq, transform=ax.transAxes, va="top", ha="left", fontsize=8)
    ax.set_xlabel("Stringency (GT recall)")
    ax.set_ylabel("Mean F1")
    ax.grid(
        True,
        linestyle=str(cfg("grid.by_protocol_linestyle", ":")),
        linewidth=cfg_float("grid.by_protocol_linewidth", 0.3),
        alpha=cfg_float("grid.overall_alpha", 0.6),
    )

    out_path = topic_figure_dir(TOPIC) / OUTPUT_NAME
    save_figure(fig, out_path)
    plt.close(fig)
    print(f"Wrote figure: {out_path}")


if __name__ == "__main__":
    main()
