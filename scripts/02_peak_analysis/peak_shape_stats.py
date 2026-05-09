"""
peak_shape_stats.py - 基于 read-level peak 记录计算 coverage 形态统计

输入：
- feather 文件（每条 read 一行），至少包含列：
  - pas
  - distance_to_pas
  - reference_length

输出（均为 feather）：
1) <prefix>_peak_shape_by_pas.feather
   每个 peak（pas）一行，包含偏度、峰度、平均距离、分位数指标等
2) <prefix>_peak_shape_summary.feather
   单行汇总，包含样本级统计（每指标 mean/std）和 pooled coverage 指标
"""

import argparse
import logging
import os
from typing import Dict, List

import numpy as np
import pandas as pd


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute peak shape statistics from read-level peaks")
    parser.add_argument("--input", required=True, help="Input feather file from peak extraction")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--prefix", required=True, help="Output prefix")
    parser.add_argument("--window_start", type=int, default=-400, help="Window start (relative to PAS)")
    parser.add_argument("--window_end", type=int, default=20, help="Window end (relative to PAS)")
    parser.add_argument("--min_reads_per_pas", type=int, default=50, help="Minimum reads per PAS")
    return parser.parse_args()


def weighted_quantiles(values: np.ndarray, weights: np.ndarray, probs: List[float]) -> np.ndarray:
    positive = weights > 0
    if not np.any(positive):
        return np.full(len(probs), np.nan, dtype=np.float64)
    values = values[positive]
    weights = weights[positive]
    cdf = np.cumsum(weights)
    return np.interp(probs, cdf, values)


def reads_to_coverage(
    distance_to_pas: np.ndarray,
    reference_length: np.ndarray,
    window_start: int,
    window_end: int,
) -> np.ndarray:
    """
    将 read-level 记录转换为窗口内 coverage。

    coverage 定义：
    对每条 read，覆盖区间为 [distance_to_pas - reference_length, distance_to_pas)（半开区间）。
    """
    n_bins = window_end - window_start + 1
    diff = np.zeros(n_bins + 1, dtype=np.int64)

    start_pos = distance_to_pas - reference_length
    end_pos = distance_to_pas

    start_idx = np.clip(start_pos - window_start, 0, n_bins).astype(np.int64)
    end_idx = np.clip(end_pos - window_start, 0, n_bins).astype(np.int64)
    valid = start_idx < end_idx
    if np.any(valid):
        np.add.at(diff, start_idx[valid], 1)
        np.add.at(diff, end_idx[valid], -1)
    return np.cumsum(diff[:-1])


def compute_distribution_metrics(
    coverage: np.ndarray,
    positions: np.ndarray,
) -> Dict[str, float]:
    total = float(np.sum(coverage))
    if total <= 0:
        return {}

    weights = coverage / total
    avg_dist = float(np.sum(positions * weights))
    avg_abs_dist = float(np.sum(np.abs(positions) * weights))

    centered = positions - avg_dist
    variance = float(np.sum(weights * centered**2))
    std = float(np.sqrt(max(variance, 0.0)))

    if std > 0:
        skewness = float(np.sum(weights * centered**3) / (std**3))
        kurtosis = float(np.sum(weights * centered**4) / (variance**2))
    else:
        skewness = np.nan
        kurtosis = np.nan

    q10, q50, q90 = weighted_quantiles(
        positions.astype(np.float64),
        weights.astype(np.float64),
        probs=[0.1, 0.5, 0.9],
    )
    effective_width = float(q90 - q10)
    left = float(q50 - q10)
    right = float(q90 - q50)
    eps = 1e-12
    asym_ratio = (left + eps) / (right + eps)
    asym_log2 = float(np.log2(asym_ratio))

    return {
        "average_distance": avg_dist,
        "average_absolute_distance": avg_abs_dist,
        "skewness": skewness,
        "kurtosis": kurtosis,
        "excess_kurtosis": float(kurtosis - 3.0) if np.isfinite(kurtosis) else np.nan,
        "std_distance": std,
        "variance_distance": variance,
        "q10": float(q10),
        "q50": float(q50),
        "q90": float(q90),
        "center_offset_q50": float(q50),
        "effective_width_q10_q90": effective_width,
        "upstream_asymmetry_ratio": float(asym_ratio),
        "upstream_asymmetry_log2": asym_log2,
        "coverage_sum": total,
    }


def summarize_metrics(df: pd.DataFrame, metric_cols: List[str]) -> Dict[str, Dict[str, float]]:
    summary = {}
    for col in metric_cols:
        values = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=np.float64)
        valid = np.isfinite(values)
        if not np.any(valid):
            summary[col] = {"mean": np.nan, "std": np.nan}
            continue
        std = np.std(values[valid], ddof=1) if np.sum(valid) > 1 else np.nan
        summary[col] = {
            "mean": float(np.mean(values[valid])),
            "std": float(std) if np.isfinite(std) else np.nan,
        }
    return summary


def flatten_summary(
    input_file: str,
    n_peaks_raw: int,
    n_peaks_after_filter: int,
    n_peaks_with_metrics: int,
    window_start: int,
    window_end: int,
    per_peak_summary: Dict[str, Dict[str, float]],
    pooled_metrics: Dict[str, float],
) -> pd.DataFrame:
    row = {
        "input_file": input_file,
        "n_peaks_raw": int(n_peaks_raw),
        "n_peaks_after_filter": int(n_peaks_after_filter),
        "n_peaks_with_metrics": int(n_peaks_with_metrics),
        "window_start": int(window_start),
        "window_end": int(window_end),
        "q10_prob": 0.1,
        "q50_prob": 0.5,
        "q90_prob": 0.9,
    }

    for metric, stat_values in per_peak_summary.items():
        row[f"per_peak_{metric}_mean"] = stat_values.get("mean", np.nan)
        row[f"per_peak_{metric}_std"] = stat_values.get("std", np.nan)

    for metric, value in pooled_metrics.items():
        row[f"pooled_{metric}"] = value

    return pd.DataFrame([row])


def main() -> None:
    args = parse_args()
    if args.window_end < args.window_start:
        raise ValueError("window_end must be greater than or equal to window_start")
    os.makedirs(args.output, exist_ok=True)

    logger.info("Loading read-level peaks: %s", args.input)
    df = pd.read_feather(args.input)
    required_cols = {"pas", "distance_to_pas", "reference_length"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    raw_peak_count = int(df["pas"].nunique())
    pas_sizes = df.groupby("pas").size()
    valid_pas = pas_sizes[pas_sizes > args.min_reads_per_pas].index
    df = df[df["pas"].isin(valid_pas)].copy()
    kept_peak_count = int(df["pas"].nunique())

    if kept_peak_count == 0:
        raise ValueError(
            f"No PAS left after filtering min_reads_per_pas>{args.min_reads_per_pas}."
        )

    positions = np.arange(args.window_start, args.window_end + 1, dtype=np.float64)
    pooled_coverage = np.zeros_like(positions, dtype=np.float64)

    per_peak_rows = []
    logger.info("Converting reads to coverage and computing per-peak metrics...")
    for pas, peak_df in df.groupby("pas", sort=False):
        numeric_df = peak_df[["distance_to_pas", "reference_length"]].apply(
            pd.to_numeric, errors="coerce"
        )
        numeric_df = numeric_df.dropna()
        if numeric_df.empty:
            continue

        distance = numeric_df["distance_to_pas"].to_numpy(dtype=np.int64)
        ref_len = numeric_df["reference_length"].to_numpy(dtype=np.int64)

        coverage = reads_to_coverage(
            distance_to_pas=distance,
            reference_length=ref_len,
            window_start=args.window_start,
            window_end=args.window_end,
        ).astype(np.float64)
        metrics = compute_distribution_metrics(coverage, positions)
        if not metrics:
            continue

        pooled_coverage += coverage
        metrics["pas"] = pas
        metrics["read_count"] = int(len(numeric_df))
        per_peak_rows.append(metrics)

    if not per_peak_rows:
        raise ValueError("No valid peak metrics were computed from the filtered dataset.")

    by_pas_df = pd.DataFrame(per_peak_rows)
    metric_cols = [
        "average_distance",
        "average_absolute_distance",
        "skewness",
        "kurtosis",
        "excess_kurtosis",
        "std_distance",
        "variance_distance",
        "q10",
        "q50",
        "q90",
        "center_offset_q50",
        "effective_width_q10_q90",
        "upstream_asymmetry_ratio",
        "upstream_asymmetry_log2",
        "coverage_sum",
    ]

    per_peak_summary = summarize_metrics(by_pas_df, metric_cols)
    pooled_metrics = compute_distribution_metrics(pooled_coverage, positions)
    summary_df = flatten_summary(
        input_file=args.input,
        n_peaks_raw=raw_peak_count,
        n_peaks_after_filter=kept_peak_count,
        n_peaks_with_metrics=len(by_pas_df),
        window_start=args.window_start,
        window_end=args.window_end,
        per_peak_summary=per_peak_summary,
        pooled_metrics=pooled_metrics,
    )

    by_pas_path = os.path.join(args.output, f"{args.prefix}_peak_shape_by_pas.feather")
    summary_path = os.path.join(args.output, f"{args.prefix}_peak_shape_summary.feather")
    by_pas_df.to_feather(by_pas_path)
    summary_df.to_feather(summary_path)

    logger.info("Saved per-peak metrics: %s", by_pas_path)
    logger.info("Saved summary metrics: %s", summary_path)


if __name__ == "__main__":
    main()
