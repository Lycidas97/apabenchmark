#!/usr/bin/env python3
"""Render protocol overlay of normalized peak profiles with Gaussian smoothing."""

from __future__ import annotations

import argparse
from statistics import NormalDist
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import PROTOCOL_ORDER
from _shared.io import read_table, write_tsv
from _shared.paths import topic_figure_dir, topic_intermediate_dir
from _shared.style import (
    apply_reference_style,
    draw_placeholder,
    figsize_from_mm,
    figsize_from_preset,
    palette,
    save_figure,
)

TOPIC = "peak_model_params"
INPUT_NAME = "protocol_normalized_peak_profiles.tsv"
OUTPUT_NAME = "protocol_normalized_peak_profiles_overlay.pdf"
SUMMARY_NAME = "protocol_normalized_peak_profiles_overlay_summary.tsv"
SMOOTHED_NAME = "protocol_normalized_peak_profiles_smoothed.tsv"
REQUIRED_COLUMNS = {"sample_full_id", "protocol", "position", "coverage_norm"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot protocol-wise overlay of normalized peak profile curves"
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=topic_intermediate_dir(TOPIC) / INPUT_NAME,
        help="Input long-format normalized profile TSV",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=topic_figure_dir(TOPIC) / OUTPUT_NAME,
        help="Output figure path",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=topic_intermediate_dir(TOPIC) / SUMMARY_NAME,
        help="Output TSV for protocol-level mean/CI curves",
    )
    parser.add_argument(
        "--smoothed-output",
        type=Path,
        default=topic_intermediate_dir(TOPIC) / SMOOTHED_NAME,
        help="Output TSV for smoothed sample-level curves",
    )
    parser.add_argument(
        "--sigma",
        type=float,
        default=10.0,
        help="Gaussian smoothing sigma in nt (set 0 to disable smoothing)",
    )
    parser.add_argument(
        "--ci-level",
        type=float,
        default=0.95,
        help="Confidence interval level for shaded band",
    )
    parser.add_argument(
        "--fig-preset",
        type=str,
        default="full_wide",
        help="Figure size preset from _shared.style.SIZE_PRESETS_MM",
    )
    parser.add_argument("--width-mm", type=float, default=None, help="Manual figure width in mm")
    parser.add_argument("--height-mm", type=float, default=None, help="Manual figure height in mm")
    parser.add_argument("--legend-ncol", type=int, default=4, help="Legend column count")
    parser.add_argument("--line-width", type=float, default=1.0, help="Line width")
    parser.add_argument("--band-alpha", type=float, default=0.18, help="Shaded CI alpha")
    return parser.parse_args()


def resolve_figsize(args: argparse.Namespace) -> tuple[float, float]:
    if (args.width_mm is None) ^ (args.height_mm is None):
        raise ValueError("width-mm and height-mm must be provided together")
    if args.width_mm is not None and args.height_mm is not None:
        return figsize_from_mm(args.width_mm, args.height_mm)
    return figsize_from_preset(args.fig_preset)


def gaussian_smooth(values: np.ndarray, sigma: float) -> np.ndarray:
    if sigma <= 0:
        smoothed = values.astype(np.float64, copy=True)
    else:
        radius = max(1, int(round(4.0 * sigma)))
        x = np.arange(-radius, radius + 1, dtype=np.float64)
        kernel = np.exp(-0.5 * (x / sigma) ** 2)
        kernel /= kernel.sum()
        smoothed = np.convolve(values, kernel, mode="same")

    total = float(np.sum(smoothed))
    if total > 0:
        smoothed = smoothed / total
    return smoothed


def smooth_per_sample(df: pd.DataFrame, sigma: float) -> pd.DataFrame:
    positions = np.sort(df["position"].unique())
    parts: list[pd.DataFrame] = []

    for (sample_full_id, protocol), group in df.groupby(["sample_full_id", "protocol"], sort=False):
        series = (
            group[["position", "coverage_norm"]]
            .copy()
            .drop_duplicates(subset=["position"], keep="last")
            .set_index("position")["coverage_norm"]
            .astype(np.float64)
            .reindex(positions, fill_value=0.0)
        )
        smoothed = gaussian_smooth(series.to_numpy(), sigma=sigma)
        parts.append(
            pd.DataFrame(
                {
                    "sample_full_id": sample_full_id,
                    "protocol": protocol,
                    "position": positions,
                    "coverage_norm_smoothed": smoothed,
                }
            )
        )

    if not parts:
        return pd.DataFrame(
            columns=["sample_full_id", "protocol", "position", "coverage_norm_smoothed"]
        )
    return pd.concat(parts, ignore_index=True)


def summarize_protocol_curves(smoothed_df: pd.DataFrame, ci_level: float) -> pd.DataFrame:
    if smoothed_df.empty:
        return pd.DataFrame(
            columns=["protocol", "position", "mean", "std", "n_samples", "sem", "ci_low", "ci_high"]
        )

    summary = (
        smoothed_df.groupby(["protocol", "position"], as_index=False)
        .agg(
            mean=("coverage_norm_smoothed", "mean"),
            std=("coverage_norm_smoothed", "std"),
            n_samples=("sample_full_id", "nunique"),
        )
        .sort_values(["protocol", "position"], kind="stable")
        .reset_index(drop=True)
    )

    summary["std"] = summary["std"].fillna(0.0)
    summary["sem"] = summary["std"] / np.sqrt(summary["n_samples"].clip(lower=1))

    ci_level = float(np.clip(ci_level, 1e-9, 1 - 1e-9))
    z = NormalDist().inv_cdf(0.5 + ci_level / 2.0)
    summary["ci_low"] = (summary["mean"] - z * summary["sem"]).clip(lower=0.0)
    summary["ci_high"] = (summary["mean"] + z * summary["sem"]).clip(lower=0.0)
    return summary


def ordered_protocols(values: pd.Series) -> list[str]:
    present = pd.Index(values.dropna().astype(str).unique())
    ordered = [p for p in PROTOCOL_ORDER if p in present]
    extra = sorted([p for p in present if p not in PROTOCOL_ORDER])
    return ordered + extra


def main() -> None:
    args = parse_args()
    apply_reference_style()

    fig, ax = plt.subplots(figsize=resolve_figsize(args))

    try:
        df = read_table(args.input)
    except Exception as exc:
        draw_placeholder(ax, "Protocol Normalized Peak Profiles", f"Could not read input: {exc}")
        save_figure(fig, args.output)
        plt.close(fig)
        print(f"Wrote placeholder figure: {args.output}")
        return

    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        draw_placeholder(
            ax,
            "Protocol Normalized Peak Profiles",
            f"Missing required columns: {sorted(missing)}",
        )
        save_figure(fig, args.output)
        plt.close(fig)
        print(f"Wrote placeholder figure: {args.output}")
        return

    work = df.copy()
    work["position"] = pd.to_numeric(work["position"], errors="coerce")
    work["coverage_norm"] = pd.to_numeric(work["coverage_norm"], errors="coerce")
    work = work.dropna(subset=["sample_full_id", "protocol", "position", "coverage_norm"])
    work["position"] = work["position"].astype(np.int64)

    smoothed = smooth_per_sample(work, sigma=float(args.sigma))
    summary = summarize_protocol_curves(smoothed, ci_level=float(args.ci_level))

    write_tsv(smoothed, args.smoothed_output)
    write_tsv(summary, args.summary_output)

    protocols = ordered_protocols(summary["protocol"]) if not summary.empty else []
    protocol_n = (
        smoothed[["protocol", "sample_full_id"]]
        .drop_duplicates()
        .groupby("protocol")
        .size()
        .to_dict()
    )

    if not protocols:
        draw_placeholder(
            ax,
            "Protocol Normalized Peak Profiles",
            "No valid curves available after filtering.",
        )
    else:
        colors = palette(len(protocols))
        for color, protocol in zip(colors, protocols):
            cur = summary[summary["protocol"] == protocol].sort_values("position")
            if cur.empty:
                continue
            n_samples = int(protocol_n.get(protocol, 0))
            ax.plot(
                cur["position"],
                cur["mean"],
                color=color,
                linewidth=args.line_width,
                label=f"{protocol} (n={n_samples})",
            )
            ax.fill_between(
                cur["position"],
                cur["ci_low"],
                cur["ci_high"],
                color=color,
                alpha=args.band_alpha,
                linewidth=0,
            )

        ax.axvline(0, color="gray", linestyle="--", linewidth=0.6)
        ax.set_xlabel("Distance to PAS (nt)")
        ax.set_ylabel("Normalized coverage")
        ax.set_title(
            f"Protocol-normalized peak profiles (Gaussian sigma={args.sigma:g})",
            loc="left",
            fontsize=7,
            pad=4,
        )
        legend_ncol = min(max(1, int(args.legend_ncol)), max(1, len(protocols)))
        ax.legend(
            frameon=False,
            ncol=legend_ncol,
            loc="upper center",
            bbox_to_anchor=(0.5, 1.24),
            fontsize=6,
            title="Protocol",
            title_fontsize=6,
        )
        fig.subplots_adjust(left=0.09, right=0.98, bottom=0.16, top=0.78)
        ax.grid(True, linestyle=":", linewidth=0.3, alpha=0.5)

    save_figure(fig, args.output)
    plt.close(fig)

    print(f"Wrote figure: {args.output}")
    print(f"Wrote smoothed sample curves: {args.smoothed_output}")
    print(f"Wrote protocol summary curves: {args.summary_output}")


if __name__ == "__main__":
    main()
