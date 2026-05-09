#!/usr/bin/env python3
"""Fit workflow-by-peak-morphology fixed-effect models."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import f as f_dist
from scipy.stats import t as t_dist

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_SCRIPT_DIR = SCRIPT_DIR.parent
PLOT_ROOT = SIM_SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.io import ensure_dir, write_json, write_tsv  # noqa: E402

OUTPUT_DIR = SCRIPT_DIR / "output"
MAIN_INPUT = OUTPUT_DIR / "de_apa_descriptor_performance.parquet"
SINGLE_FILTER_INPUT = OUTPUT_DIR / "de_apa_single_filter_descriptor_performance.parquet"

REFERENCE_WORKFLOW = "scAPAtrap"
METRICS = ["precision", "recall", "f1"]
TRANSFORMS = ["logit_clip", "identity"]
PRIMARY_TRANSFORM = "logit_clip"
LOGIT_EPS = 1e-4

DESCRIPTORS = [
    "shape_per_peak_upstream_asymmetry_log2_mean",
    "shape_per_peak_effective_width_q10_q90_mean",
    "shape_per_peak_average_distance_mean",
]
DESCRIPTOR_Z = [f"{name}_z" for name in DESCRIPTORS]
DESCRIPTOR_LABELS = {
    "shape_per_peak_upstream_asymmetry_log2_mean_z": "upstream_asymmetry",
    "shape_per_peak_effective_width_q10_q90_mean_z": "effective_width",
    "shape_per_peak_average_distance_mean_z": "mean_distance_to_pas",
}

REGIMES = {
    "main_text_three_criteria": {
        "path": MAIN_INPUT,
        "instance_cols": ["sample", "match_type"],
        "note": "PAS-resolved workflows under selected shared criteria, criterion-averaged.",
    },
    "dapars2_scmapa_compatible_single_filter": {
        "path": SINGLE_FILTER_INPUT,
        "instance_cols": ["sample", "match_type_harmonized", "filter_type_1", "filter_type_2"],
        "note": "Harmonized endpoint comparison including gene-level endpoint workflows.",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--main-input", type=Path, default=MAIN_INPUT)
    parser.add_argument("--single-filter-input", type=Path, default=SINGLE_FILTER_INPUT)
    parser.add_argument("--logit-eps", type=float, default=LOGIT_EPS)
    return parser.parse_args()


def bh_fdr(p_values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p_values.notna() & np.isfinite(p_values)
    if not valid.any():
        return out
    p = p_values.loc[valid].astype(float)
    order = np.argsort(p.to_numpy())
    ranked = p.to_numpy()[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    result = np.empty(n, dtype=float)
    result[order] = adjusted
    out.loc[p.index] = result
    return out


def safe_logit(values: pd.Series, eps: float) -> tuple[pd.Series, dict[str, Any]]:
    raw = pd.to_numeric(values, errors="coerce")
    clipped = raw.clip(eps, 1.0 - eps)
    transformed = np.log(clipped / (1.0 - clipped))
    non_na = raw.dropna()
    diagnostics = {
        "n_nonmissing": int(non_na.shape[0]),
        "n_zero": int((non_na == 0).sum()),
        "n_one": int((non_na == 1).sum()),
        "n_clipped_low": int((non_na <= eps).sum()),
        "n_clipped_high": int((non_na >= 1.0 - eps).sum()),
    }
    denom = max(1, diagnostics["n_nonmissing"])
    diagnostics.update(
        {
            "frac_zero": diagnostics["n_zero"] / denom,
            "frac_one": diagnostics["n_one"] / denom,
            "frac_clipped_low": diagnostics["n_clipped_low"] / denom,
            "frac_clipped_high": diagnostics["n_clipped_high"] / denom,
        }
    )
    return transformed, diagnostics


def make_instance_id(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    missing = [col for col in cols if col not in df.columns]
    if missing:
        raise KeyError(f"Missing instance columns: {missing}")
    return df[cols].astype("string").fillna("<NA>").agg("||".join, axis=1)


def unique_sample_descriptors(df: pd.DataFrame) -> pd.DataFrame:
    cols = ["empirical_sample_id", "protocol", *DESCRIPTOR_Z]
    missing = [col for col in cols if col not in df.columns]
    if missing:
        raise KeyError(f"Missing descriptor columns: {missing}")
    return df[cols].drop_duplicates(subset=["empirical_sample_id"]).dropna(subset=DESCRIPTOR_Z).copy()


def add_protocol_centered_descriptors(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    sample_df = unique_sample_descriptors(out)
    centered_cols = []
    for col in DESCRIPTOR_Z:
        centered = sample_df[col] - sample_df.groupby("protocol")[col].transform("mean")
        sd = centered.std(ddof=0)
        centered_col = f"{col}_protocol_centered_z"
        sample_df[centered_col] = (centered - centered.mean()) / sd if pd.notna(sd) and sd > 0 else np.nan
        centered_cols.append(centered_col)
    return out.merge(sample_df[["empirical_sample_id", *centered_cols]], on="empirical_sample_id", how="left")


def add_pca_scores(df: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    out = df.copy()
    sample_df = unique_sample_descriptors(out)
    x = sample_df[DESCRIPTOR_Z].to_numpy(dtype=float)
    x = x - x.mean(axis=0, keepdims=True)
    _, s, vt = np.linalg.svd(x, full_matrices=False)
    scores = x @ vt.T
    var = (s**2) / max(1, x.shape[0] - 1)
    explained = var / var.sum() if var.sum() > 0 else np.full_like(var, np.nan)
    pca_rows: list[dict[str, Any]] = []
    for pc_idx in range(min(2, scores.shape[1])):
        col = f"morphology_pc{pc_idx + 1}"
        sample_df[col] = scores[:, pc_idx]
        pca_rows.append(
            {
                "diagnostic_type": "pca_explained_variance",
                "component": f"PC{pc_idx + 1}",
                "value": float(explained[pc_idx]),
            }
        )
        for desc_col, loading in zip(DESCRIPTOR_Z, vt[pc_idx, :]):
            pca_rows.append(
                {
                    "diagnostic_type": "pca_loading",
                    "component": f"PC{pc_idx + 1}",
                    "descriptor": DESCRIPTOR_LABELS[desc_col],
                    "value": float(loading),
                }
            )
    return out.merge(sample_df[["empirical_sample_id", "morphology_pc1", "morphology_pc2"]], on="empirical_sample_id", how="left"), pca_rows


def descriptor_diagnostics(regime: str, df: pd.DataFrame) -> list[dict[str, Any]]:
    sample_df = unique_sample_descriptors(df)
    rows: list[dict[str, Any]] = []
    corr = sample_df[DESCRIPTOR_Z].corr()
    for a in DESCRIPTOR_Z:
        for b in DESCRIPTOR_Z:
            rows.append(
                {
                    "regime": regime,
                    "diagnostic_type": "descriptor_correlation",
                    "descriptor": DESCRIPTOR_LABELS[a],
                    "descriptor_2": DESCRIPTOR_LABELS[b],
                    "value": float(corr.loc[a, b]),
                    "n_samples": int(len(sample_df)),
                }
            )
    x_all = sample_df[DESCRIPTOR_Z].to_numpy(dtype=float)
    for idx, col in enumerate(DESCRIPTOR_Z):
        y = x_all[:, idx]
        x = np.delete(x_all, idx, axis=1)
        x = np.column_stack([np.ones(x.shape[0]), x])
        beta = np.linalg.lstsq(x, y, rcond=None)[0]
        resid = y - x @ beta
        rss = float(np.sum(resid**2))
        tss = float(np.sum((y - y.mean()) ** 2))
        r2 = 1.0 - rss / tss if tss > 0 else np.nan
        vif = 1.0 / (1.0 - r2) if pd.notna(r2) and r2 < 1 else np.inf
        rows.append(
            {
                "regime": regime,
                "diagnostic_type": "vif",
                "descriptor": DESCRIPTOR_LABELS[col],
                "value": float(vif),
                "r2": float(r2),
                "n_samples": int(len(sample_df)),
            }
        )
    protocol_counts = sample_df.groupby("protocol").size()
    rows.extend(
        [
            {
                "regime": regime,
                "diagnostic_type": "protocol_sample_count",
                "metric_name": "n_protocols",
                "value": int(protocol_counts.shape[0]),
            },
            {
                "regime": regime,
                "diagnostic_type": "protocol_sample_count",
                "metric_name": "min_samples_per_protocol",
                "value": int(protocol_counts.min()) if not protocol_counts.empty else 0,
            },
            {
                "regime": regime,
                "diagnostic_type": "protocol_sample_count",
                "metric_name": "median_samples_per_protocol",
                "value": float(protocol_counts.median()) if not protocol_counts.empty else np.nan,
            },
            {
                "regime": regime,
                "diagnostic_type": "protocol_sample_count",
                "metric_name": "n_protocols_single_sample",
                "value": int((protocol_counts == 1).sum()),
            },
        ]
    )
    return rows


def demean_by_instance(y: np.ndarray, x: np.ndarray, instance: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    frame = pd.DataFrame({"_y": y, "_instance": instance.to_numpy()})
    y_dm = y - frame.groupby("_instance")["_y"].transform("mean").to_numpy()
    if x.shape[1] == 0:
        return y_dm, x
    x_df = pd.DataFrame(x)
    x_df["_instance"] = instance.to_numpy()
    means = x_df.groupby("_instance").transform("mean").drop(columns=[], errors="ignore")
    x_dm = x - means.iloc[:, : x.shape[1]].to_numpy()
    return y_dm, x_dm


def fit_ols(y: np.ndarray, x: np.ndarray, instance: pd.Series) -> dict[str, Any]:
    y_dm, x_dm = demean_by_instance(y, x, instance)
    n = len(y_dm)
    g = int(instance.nunique())
    xtx = x_dm.T @ x_dm if x_dm.shape[1] else np.zeros((0, 0))
    rank = int(np.linalg.matrix_rank(xtx)) if x_dm.shape[1] else 0
    if x_dm.shape[1] == 0 or rank == 0:
        fitted = np.zeros_like(y_dm)
        beta = np.zeros(x_dm.shape[1])
        xtx_inv = np.zeros((x_dm.shape[1], x_dm.shape[1]))
    else:
        # The model has many rows but few columns; solving through X'X avoids
        # repeated large SVDs while preserving rank-deficient pseudoinverse behavior.
        xtx_inv = np.linalg.pinv(xtx)
        beta = xtx_inv @ (x_dm.T @ y_dm)
        fitted = x_dm @ beta
    resid = y_dm - fitted
    rss = float(np.sum(resid**2))
    df_resid = int(n - g - rank)
    sigma2 = rss / df_resid if df_resid > 0 else np.nan
    cov = sigma2 * xtx_inv if x_dm.shape[1] else np.zeros((0, 0))
    return {
        "beta": beta,
        "cov": cov,
        "rss": rss,
        "rank": rank,
        "df_resid": df_resid,
        "n_rows": n,
        "n_instances": g,
    }


def workflow_order(df: pd.DataFrame) -> list[str]:
    tools = sorted(df["tool"].dropna().astype(str).unique())
    if REFERENCE_WORKFLOW not in tools:
        raise ValueError(f"Reference workflow {REFERENCE_WORKFLOW!r} is missing.")
    return [REFERENCE_WORKFLOW] + [tool for tool in tools if tool != REFERENCE_WORKFLOW]


def workflow_design(df: pd.DataFrame, tools: list[str]) -> tuple[np.ndarray, list[dict[str, str]]]:
    nonref = tools[1:]
    cols = [(df["tool"].astype(str) == tool).astype(float).to_numpy() for tool in nonref]
    meta = [{"term_type": "workflow_main", "workflow": tool, "term": f"workflow[{tool}]"} for tool in nonref]
    return (np.column_stack(cols) if cols else np.empty((len(df), 0))), meta


def interaction_design(
    df: pd.DataFrame,
    tools: list[str],
    variables: list[str],
    *,
    term_type: str,
    labels: dict[str, str] | None = None,
) -> tuple[np.ndarray, list[dict[str, str]]]:
    nonref = tools[1:]
    cols = []
    meta = []
    for tool in nonref:
        tool_mask = (df["tool"].astype(str) == tool).astype(float).to_numpy()
        for var in variables:
            cols.append(tool_mask * pd.to_numeric(df[var], errors="coerce").to_numpy(dtype=float))
            label = labels.get(var, var) if labels else var
            meta.append(
                {
                    "term_type": term_type,
                    "workflow": tool,
                    "variable": var,
                    "descriptor": label,
                    "term": f"workflow[{tool}]:{label}",
                }
            )
    return (np.column_stack(cols) if cols else np.empty((len(df), 0))), meta


def protocol_interaction_design(df: pd.DataFrame, tools: list[str]) -> tuple[np.ndarray, list[dict[str, str]]]:
    protocols = sorted(df["protocol"].dropna().astype(str).unique())
    if not protocols:
        return np.empty((len(df), 0)), []
    ref_protocol = protocols[0]
    nonref_protocols = [p for p in protocols if p != ref_protocol]
    cols = []
    meta = []
    for tool in tools[1:]:
        tool_mask = (df["tool"].astype(str) == tool).astype(float).to_numpy()
        for protocol in nonref_protocols:
            cols.append(tool_mask * (df["protocol"].astype(str) == protocol).astype(float).to_numpy())
            meta.append(
                {
                    "term_type": "workflow_x_protocol",
                    "workflow": tool,
                    "protocol": protocol,
                    "reference_protocol": ref_protocol,
                    "term": f"workflow[{tool}]:protocol[{protocol}]",
                }
            )
    return (np.column_stack(cols) if cols else np.empty((len(df), 0))), meta


def combine_designs(*parts: tuple[np.ndarray, list[dict[str, str]]]) -> tuple[np.ndarray, list[dict[str, str]]]:
    arrays = [arr for arr, _ in parts if arr.shape[1] > 0]
    metas: list[dict[str, str]] = []
    for _, meta in parts:
        metas.extend(meta)
    if not arrays:
        n = parts[0][0].shape[0] if parts else 0
        return np.empty((n, 0)), metas
    return np.column_stack(arrays), metas


def analysis_frame(
    df: pd.DataFrame,
    *,
    metric: str,
    transform: str,
    eps: float,
    instance_cols: list[str],
    required_columns: list[str],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    work = df.copy()
    n_start = len(work)
    work["_response_raw"] = pd.to_numeric(work[metric], errors="coerce")
    if transform == "logit_clip":
        work["_response"], clip_diag = safe_logit(work["_response_raw"], eps)
    elif transform == "identity":
        work["_response"] = work["_response_raw"]
        non_na = work["_response_raw"].dropna()
        clip_diag = {
            "n_nonmissing": int(non_na.shape[0]),
            "n_zero": int((non_na == 0).sum()),
            "n_one": int((non_na == 1).sum()),
            "n_clipped_low": 0,
            "n_clipped_high": 0,
            "frac_zero": float((non_na == 0).mean()) if len(non_na) else np.nan,
            "frac_one": float((non_na == 1).mean()) if len(non_na) else np.nan,
            "frac_clipped_low": 0.0,
            "frac_clipped_high": 0.0,
        }
    else:
        raise ValueError(f"Unknown transform: {transform}")
    needed = ["_response", "tool", *instance_cols, *required_columns]
    if "protocol" in work.columns:
        needed.append("protocol")
    needed = list(dict.fromkeys(needed))
    work = work.dropna(subset=[col for col in needed if col in work.columns]).copy()
    if "_instance_id" not in work.columns:
        work["_instance_id"] = make_instance_id(work, instance_cols)
    tool_counts = work.groupby("_instance_id")["tool"].nunique()
    keep_instances = set(tool_counts[tool_counts >= 2].index)
    n_lt2 = int((tool_counts < 2).sum())
    work = work[work["_instance_id"].isin(keep_instances)].copy()
    diagnostics = {
        **clip_diag,
        "n_rows_start": int(n_start),
        "n_rows_dropped_missing_metric": int(df[metric].isna().sum()) if metric in df.columns else n_start,
        "n_rows_after_missing_filters": int(len(work)),
        "n_instances_dropped_lt2_tools": n_lt2,
    }
    return work, diagnostics


def nested_test(reduced: dict[str, Any], full: dict[str, Any]) -> dict[str, Any]:
    df_num = int(full["rank"] - reduced["rank"])
    df_den = int(full["n_rows"] - full["n_instances"] - full["rank"])
    flag = "ok"
    f_stat = np.nan
    p_value = np.nan
    partial_r2 = np.nan
    if df_num <= 0:
        flag = "non_positive_df_num"
    elif df_den <= 0:
        flag = "non_positive_df_den"
    elif reduced["rss"] <= 0:
        flag = "rss_reduced_nonpositive"
    elif full["rss"] > reduced["rss"] + 1e-8:
        flag = "rss_full_greater_than_reduced"
    else:
        partial_r2 = (reduced["rss"] - full["rss"]) / reduced["rss"]
        f_stat = ((reduced["rss"] - full["rss"]) / df_num) / (full["rss"] / df_den)
        p_value = float(f_dist.sf(f_stat, df_num, df_den)) if np.isfinite(f_stat) else np.nan
        if full["rank"] < len(full.get("term_meta", [])):
            flag = "rank_deficient_but_testable"
    return {
        "rank_reduced": int(reduced["rank"]),
        "rank_full": int(full["rank"]),
        "df_num": df_num,
        "df_den": df_den,
        "rss_reduced": float(reduced["rss"]),
        "rss_full": float(full["rss"]),
        "partial_r2": float(partial_r2) if np.isfinite(partial_r2) else np.nan,
        "f_stat": float(f_stat) if np.isfinite(f_stat) else np.nan,
        "p_value": p_value,
        "diagnostic_flag": flag,
    }


def coefficient_rows(
    *,
    regime: str,
    metric: str,
    transform: str,
    model_name: str,
    model_family: str,
    tested_effect: str,
    fit: dict[str, Any],
    term_meta: list[dict[str, str]],
) -> list[dict[str, Any]]:
    rows = []
    beta = fit["beta"]
    cov = fit["cov"]
    df_resid = fit["df_resid"]
    for idx, meta in enumerate(term_meta):
        if meta.get("term_type") != "workflow_x_univariate_descriptor":
            continue
        se = float(np.sqrt(cov[idx, idx])) if cov.size and cov[idx, idx] >= 0 else np.nan
        estimate = float(beta[idx])
        t_stat = estimate / se if pd.notna(se) and se > 0 else np.nan
        p_value = float(2 * t_dist.sf(abs(t_stat), df_resid)) if np.isfinite(t_stat) and df_resid > 0 else np.nan
        rows.append(
            {
                "regime": regime,
                "metric": metric,
                "response_transform": transform,
                "model": model_name,
                "model_family": model_family,
                "tested_effect": tested_effect,
                "workflow": meta.get("workflow"),
                "descriptor": meta.get("descriptor", meta.get("variable")),
                "coefficient": estimate,
                "standard_error": se,
                "t_stat": t_stat,
                "p_value": p_value,
                "reference_workflow": REFERENCE_WORKFLOW,
                "coefficient_interpretation": (
                    "standardized descriptor slope difference relative to scAPAtrap "
                    "in a single-descriptor benchmark-instance-demeaned model"
                ),
            }
        )
    return rows


def posthoc_rows(coef_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for row in coef_rows:
        if row.get("model_family") != "univariate_descriptor":
            continue
        rows.append(
            {
                **row,
                "contrast_type": "workflow_descriptor_slope_difference",
                "contrast": f"{row['workflow']} - {REFERENCE_WORKFLOW}",
                "estimate": row["coefficient"],
            }
        )
    return rows


def fit_internal_descriptor_model(y: np.ndarray, x: np.ndarray) -> dict[str, Any]:
    n = len(y)
    if n < 3:
        return {
            "intercept": np.nan,
            "slope": np.nan,
            "standard_error": np.nan,
            "t_stat": np.nan,
            "p_value": np.nan,
            "rss_reduced": np.nan,
            "rss_full": np.nan,
            "r2": np.nan,
            "partial_r2": np.nan,
            "df_num": 1,
            "df_den": n - 2,
            "diagnostic_flag": "insufficient_rows",
        }
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        return {
            "intercept": np.nan,
            "slope": np.nan,
            "standard_error": np.nan,
            "t_stat": np.nan,
            "p_value": np.nan,
            "rss_reduced": np.nan,
            "rss_full": np.nan,
            "r2": np.nan,
            "partial_r2": np.nan,
            "df_num": 1,
            "df_den": n - 2,
            "diagnostic_flag": "nonfinite_values",
        }
    if float(np.var(x)) <= 0:
        y_mean = float(np.mean(y))
        rss_reduced = float(np.sum((y - y_mean) ** 2))
        return {
            "intercept": y_mean,
            "slope": np.nan,
            "standard_error": np.nan,
            "t_stat": np.nan,
            "p_value": np.nan,
            "rss_reduced": rss_reduced,
            "rss_full": np.nan,
            "r2": np.nan,
            "partial_r2": np.nan,
            "df_num": 1,
            "df_den": n - 2,
            "diagnostic_flag": "zero_descriptor_variance",
        }

    design = np.column_stack([np.ones(n), x])
    xtx_inv = np.linalg.pinv(design.T @ design)
    beta = xtx_inv @ (design.T @ y)
    fitted = design @ beta
    residuals = y - fitted
    rss_full = float(np.sum(residuals**2))
    y_mean = float(np.mean(y))
    rss_reduced = float(np.sum((y - y_mean) ** 2))
    df_den = n - 2
    if df_den <= 0:
        flag = "non_positive_df_den"
        sigma2 = np.nan
    else:
        flag = "ok"
        sigma2 = rss_full / df_den
    se = float(np.sqrt(sigma2 * xtx_inv[1, 1])) if np.isfinite(sigma2) and xtx_inv[1, 1] >= 0 else np.nan
    slope = float(beta[1])
    t_stat = slope / se if pd.notna(se) and se > 0 else np.nan
    p_value = float(2 * t_dist.sf(abs(t_stat), df_den)) if np.isfinite(t_stat) and df_den > 0 else np.nan
    r2 = 1.0 - rss_full / rss_reduced if rss_reduced > 0 else np.nan
    partial_r2 = (rss_reduced - rss_full) / rss_reduced if rss_reduced > 0 else np.nan
    return {
        "intercept": float(beta[0]),
        "slope": slope,
        "standard_error": se,
        "t_stat": t_stat,
        "p_value": p_value,
        "rss_reduced": rss_reduced,
        "rss_full": rss_full,
        "r2": float(r2) if np.isfinite(r2) else np.nan,
        "partial_r2": float(partial_r2) if np.isfinite(partial_r2) else np.nan,
        "df_num": 1,
        "df_den": df_den,
        "diagnostic_flag": flag,
    }


def workflow_internal_variance_rows(
    *,
    regime: str,
    df: pd.DataFrame,
    metric: str,
    transform: str,
    eps: float,
    instance_cols: list[str],
) -> list[dict[str, Any]]:
    work, _ = analysis_frame(
        df,
        metric=metric,
        transform=transform,
        eps=eps,
        instance_cols=instance_cols,
        required_columns=DESCRIPTOR_Z,
    )
    validate_analysis_frame(work, regime=regime, metric=metric, transform=transform, model_name="workflow_internal_variance")
    grouped = work.groupby("_instance_id")["_response"]
    instance_sum = grouped.transform("sum")
    instance_count = grouped.transform("count")
    work = work.copy()
    work["_n_baseline_workflows"] = instance_count - 1
    work["_loo_instance_baseline"] = (instance_sum - work["_response"]) / work["_n_baseline_workflows"]
    work["_internal_response"] = work["_response"] - work["_loo_instance_baseline"]
    work = work[work["_n_baseline_workflows"] >= 1].copy()

    rows: list[dict[str, Any]] = []
    for workflow in sorted(work["tool"].dropna().astype(str).unique()):
        tool_df = work[work["tool"].astype(str) == workflow].copy()
        baseline_counts = tool_df["_n_baseline_workflows"]
        for desc_col in DESCRIPTOR_Z:
            desc_label = DESCRIPTOR_LABELS[desc_col]
            model_df = tool_df[["_instance_id", "_internal_response", desc_col, "_n_baseline_workflows"]].dropna().copy()
            n_rows = int(len(model_df))
            n_instances = int(model_df["_instance_id"].nunique()) if n_rows else 0
            if n_rows:
                fit = fit_internal_descriptor_model(
                    model_df["_internal_response"].to_numpy(dtype=float),
                    model_df[desc_col].to_numpy(dtype=float),
                )
                n_base = model_df["_n_baseline_workflows"]
                n_base_min = int(n_base.min())
                n_base_median = float(n_base.median())
                n_base_max = int(n_base.max())
            else:
                fit = {
                    "intercept": np.nan,
                    "slope": np.nan,
                    "standard_error": np.nan,
                    "t_stat": np.nan,
                    "p_value": np.nan,
                    "rss_reduced": np.nan,
                    "rss_full": np.nan,
                    "r2": np.nan,
                    "partial_r2": np.nan,
                    "df_num": 1,
                    "df_den": np.nan,
                    "diagnostic_flag": "insufficient_rows",
                }
                n_base_min = np.nan
                n_base_median = np.nan
                n_base_max = np.nan
            fdr_group = "|".join([regime, metric, transform, "workflow_internal"])
            rows.append(
                {
                    "regime": regime,
                    "metric": metric,
                    "response_transform": transform,
                    "is_primary": bool(transform == PRIMARY_TRANSFORM),
                    "workflow": workflow,
                    "descriptor": desc_label,
                    "n_rows": n_rows,
                    "n_instances": n_instances,
                    "n_baseline_workflows_min": n_base_min,
                    "n_baseline_workflows_median": n_base_median,
                    "n_baseline_workflows_max": n_base_max,
                    **fit,
                    "fdr_group": fdr_group,
                }
            )
    return rows


def add_internal_fdr(internal: pd.DataFrame) -> pd.DataFrame:
    if internal.empty:
        return internal
    internal["fdr_bh"] = np.nan
    for _, idx in internal.groupby("fdr_group").groups.items():
        internal.loc[idx, "fdr_bh"] = bh_fdr(internal.loc[idx, "p_value"])
    return internal


def validate_analysis_frame(work: pd.DataFrame, *, regime: str, metric: str, transform: str, model_name: str) -> None:
    duplicated = work.duplicated(subset=["_instance_id", "tool"], keep=False)
    if duplicated.any():
        examples = (
            work.loc[duplicated, ["_instance_id", "tool"]]
            .drop_duplicates()
            .head(10)
            .to_dict("records")
        )
        raise ValueError(
            "Duplicate instance_id x workflow rows entering model "
            f"{model_name} for {regime}/{metric}/{transform}: {examples}"
        )


def comparison_row(
    *,
    regime: str,
    metric: str,
    transform: str,
    model_family: str,
    tested_effect: str,
    descriptor: str,
    model_reduced: str,
    model_full: str,
    work: pd.DataFrame,
    reduced: dict[str, Any],
    full: dict[str, Any],
    diag: dict[str, Any],
) -> dict[str, Any]:
    fdr_group = "|".join([regime, transform, metric, model_family])
    return {
        "regime": regime,
        "metric": metric,
        "response_transform": transform,
        "is_primary": bool(transform == PRIMARY_TRANSFORM and model_family == "univariate_descriptor"),
        "model_family": model_family,
        "tested_effect": tested_effect,
        "descriptor": descriptor,
        "model_reduced": model_reduced,
        "model_full": model_full,
        "n_rows": int(len(work)),
        "n_instances": int(work["_instance_id"].nunique()),
        "n_tools": int(work["tool"].nunique()),
        **nested_test(reduced, full),
        "fdr_group": fdr_group,
        **diag,
    }


def run_regime_models(regime: str, df: pd.DataFrame, instance_cols: list[str], eps: float) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    comparison_rows: list[dict[str, Any]] = []
    coefficient_out: list[dict[str, Any]] = []
    posthoc_out: list[dict[str, Any]] = []
    instance_rows: list[dict[str, Any]] = []

    missing_instance = [col for col in instance_cols if col not in df.columns]
    if missing_instance:
        raise KeyError(f"Missing instance columns for {regime}: {missing_instance}")
    null_instance = df[instance_cols].isna().any(axis=1)
    if null_instance.any():
        examples = df.loc[null_instance, instance_cols].head(10).to_dict("records")
        raise ValueError(f"Missing instance identifiers for {regime}: {examples}")

    df = df.copy()
    df["_instance_id"] = make_instance_id(df, instance_cols)
    df = add_protocol_centered_descriptors(df)
    df, pca_diag = add_pca_scores(df)
    centered_cols = [f"{col}_protocol_centered_z" for col in DESCRIPTOR_Z]
    centered_labels = {
        centered_col: f"{DESCRIPTOR_LABELS[base_col]}_protocol_centered"
        for base_col, centered_col in zip(DESCRIPTOR_Z, centered_cols)
    }
    pca_cols = ["morphology_pc1", "morphology_pc2"]
    omnibus_descriptor = "+".join(DESCRIPTOR_LABELS[col] for col in DESCRIPTOR_Z)

    for metric in METRICS:
        for transform in TRANSFORMS:
            morph_required = DESCRIPTOR_Z
            protocol_required = ["protocol"]
            centered_required = ["protocol", *centered_cols]
            pca_required = pca_cols

            work, diag = analysis_frame(
                df,
                metric=metric,
                transform=transform,
                eps=eps,
                instance_cols=instance_cols,
                required_columns=morph_required,
            )
            validate_analysis_frame(work, regime=regime, metric=metric, transform=transform, model_name="descriptor_models")
            tools = workflow_order(work)
            y = work["_response"].to_numpy(dtype=float)
            inst = work["_instance_id"]
            wf_x, wf_meta = workflow_design(work, tools)
            fe_x = np.empty((len(work), 0))
            fit_fe = fit_ols(y, fe_x, inst)
            fit_baseline = fit_ols(y, wf_x, inst)
            fit_baseline["term_meta"] = wf_meta

            comparison_rows.append(
                comparison_row(
                    regime=regime,
                    metric=metric,
                    transform=transform,
                    model_family="workflow_main",
                    tested_effect="workflow_main",
                    descriptor="",
                    model_reduced="instance_fe_only",
                    model_full="baseline_workflow",
                    work=work,
                    reduced=fit_fe,
                    full=fit_baseline,
                    diag=diag,
                )
            )

            # Primary standalone descriptor tests.
            primary_fits: dict[str, dict[str, Any]] = {}
            primary_term_meta: dict[str, list[dict[str, str]]] = {}
            for desc_col in DESCRIPTOR_Z:
                desc_label = DESCRIPTOR_LABELS[desc_col]
                tested = f"workflow_x_{desc_label}"
                desc_x, desc_meta = interaction_design(
                    work,
                    tools,
                    [desc_col],
                    term_type="workflow_x_univariate_descriptor",
                    labels=DESCRIPTOR_LABELS,
                )
                full_x, full_meta = combine_designs((wf_x, wf_meta), (desc_x, desc_meta))
                fit_desc = fit_ols(y, full_x, inst)
                fit_desc["term_meta"] = full_meta
                primary_fits[desc_col] = fit_desc
                primary_term_meta[desc_col] = full_meta
                comparison_rows.append(
                    comparison_row(
                        regime=regime,
                        metric=metric,
                        transform=transform,
                        model_family="univariate_descriptor",
                        tested_effect=tested,
                        descriptor=desc_label,
                        model_reduced="baseline_workflow",
                        model_full=tested,
                        work=work,
                        reduced=fit_baseline,
                        full=fit_desc,
                        diag=diag,
                    )
                )
                coef_rows = coefficient_rows(
                    regime=regime,
                    metric=metric,
                    transform=transform,
                    model_name=tested,
                    model_family="univariate_descriptor",
                    tested_effect=tested,
                    fit=fit_desc,
                    term_meta=full_meta,
                )
                coefficient_out.extend(coef_rows)
                posthoc_out.extend(posthoc_rows(coef_rows))

            # Secondary omnibus morphology block.
            morph_x, morph_meta = interaction_design(
                work,
                tools,
                DESCRIPTOR_Z,
                term_type="workflow_x_morphology_omnibus",
                labels=DESCRIPTOR_LABELS,
            )
            morph_full_x, morph_full_meta = combine_designs((wf_x, wf_meta), (morph_x, morph_meta))
            fit_morph = fit_ols(y, morph_full_x, inst)
            fit_morph["term_meta"] = morph_full_meta
            comparison_rows.append(
                comparison_row(
                    regime=regime,
                    metric=metric,
                    transform=transform,
                    model_family="omnibus_morphology_block",
                    tested_effect="workflow_x_morphology_omnibus",
                    descriptor=omnibus_descriptor,
                    model_reduced="baseline_workflow",
                    model_full="workflow_x_morphology_omnibus",
                    work=work,
                    reduced=fit_baseline,
                    full=fit_morph,
                    diag=diag,
                )
            )

            # Conditional descriptor sensitivity: add one descriptor after the other two.
            for desc_col in DESCRIPTOR_Z:
                desc_label = DESCRIPTOR_LABELS[desc_col]
                other_cols = [col for col in DESCRIPTOR_Z if col != desc_col]
                other_labels = [DESCRIPTOR_LABELS[col] for col in other_cols]
                reduced_x_part, reduced_meta_part = interaction_design(
                    work,
                    tools,
                    other_cols,
                    term_type="workflow_x_conditional_descriptor_reduced",
                    labels=DESCRIPTOR_LABELS,
                )
                reduced_x, reduced_meta = combine_designs((wf_x, wf_meta), (reduced_x_part, reduced_meta_part))
                fit_reduced = fit_ols(y, reduced_x, inst)
                fit_reduced["term_meta"] = reduced_meta
                tested = f"workflow_x_{desc_label}_given_{'_'.join(other_labels)}"
                comparison_rows.append(
                    comparison_row(
                        regime=regime,
                        metric=metric,
                        transform=transform,
                        model_family="conditional_descriptor",
                        tested_effect=tested,
                        descriptor=desc_label,
                        model_reduced=f"baseline_workflow_plus_{'_'.join(other_labels)}",
                        model_full="workflow_x_morphology_omnibus",
                        work=work,
                        reduced=fit_reduced,
                        full=fit_morph,
                        diag=diag,
                    )
                )

            tool_counts = work.groupby("_instance_id")["tool"].nunique()
            instance_rows.append(
                {
                    "regime": regime,
                    "metric": metric,
                    "response_transform": transform,
                    "n_instances": int(tool_counts.shape[0]),
                    "median_tools_per_instance": float(tool_counts.median()) if not tool_counts.empty else np.nan,
                    "min_tools_per_instance": int(tool_counts.min()) if not tool_counts.empty else 0,
                    "max_tools_per_instance": int(tool_counts.max()) if not tool_counts.empty else 0,
                    "n_instances_dropped_lt2_tools": int(diag["n_instances_dropped_lt2_tools"]),
                }
            )

            # Protocol comparison.
            work_p, diag_p = analysis_frame(
                df,
                metric=metric,
                transform=transform,
                eps=eps,
                instance_cols=instance_cols,
                required_columns=protocol_required,
            )
            validate_analysis_frame(work_p, regime=regime, metric=metric, transform=transform, model_name="workflow_x_protocol")
            tools_p = workflow_order(work_p)
            y_p = work_p["_response"].to_numpy(dtype=float)
            inst_p = work_p["_instance_id"]
            wf_p, wf_meta_p = workflow_design(work_p, tools_p)
            proto_x, proto_meta = protocol_interaction_design(work_p, tools_p)
            proto_full_x, proto_full_meta = combine_designs((wf_p, wf_meta_p), (proto_x, proto_meta))
            fit_wf_p = fit_ols(y_p, wf_p, inst_p)
            fit_proto = fit_ols(y_p, proto_full_x, inst_p)
            fit_proto["term_meta"] = proto_full_meta
            comparison_rows.append(
                comparison_row(
                    regime=regime,
                    metric=metric,
                    transform=transform,
                    model_family="protocol_interaction",
                    tested_effect="workflow_x_protocol",
                    descriptor="protocol",
                    model_reduced="baseline_workflow",
                    model_full="workflow_x_protocol",
                    work=work_p,
                    reduced=fit_wf_p,
                    full=fit_proto,
                    diag=diag_p,
                )
            )

            # Protocol-centered descriptor sensitivity.
            work_c, diag_c = analysis_frame(
                df,
                metric=metric,
                transform=transform,
                eps=eps,
                instance_cols=instance_cols,
                required_columns=centered_required,
            )
            validate_analysis_frame(work_c, regime=regime, metric=metric, transform=transform, model_name="protocol_centered_descriptor")
            tools_c = workflow_order(work_c)
            y_c = work_c["_response"].to_numpy(dtype=float)
            inst_c = work_c["_instance_id"]
            wf_c, wf_meta_c = workflow_design(work_c, tools_c)
            proto_c, proto_meta_c = protocol_interaction_design(work_c, tools_c)
            reduced_c_x, reduced_c_meta = combine_designs((wf_c, wf_meta_c), (proto_c, proto_meta_c))
            fit_reduced_c = fit_ols(y_c, reduced_c_x, inst_c)
            fit_reduced_c["term_meta"] = reduced_c_meta
            for centered_col in centered_cols:
                centered_label = centered_labels[centered_col]
                centered_x, centered_meta = interaction_design(
                    work_c,
                    tools_c,
                    [centered_col],
                    term_type="workflow_x_protocol_centered_descriptor",
                    labels=centered_labels,
                )
                full_c_x, full_c_meta = combine_designs((wf_c, wf_meta_c), (proto_c, proto_meta_c), (centered_x, centered_meta))
                fit_full_c = fit_ols(y_c, full_c_x, inst_c)
                fit_full_c["term_meta"] = full_c_meta
                base_label = centered_label.replace("_protocol_centered", "")
                comparison_rows.append(
                    comparison_row(
                        regime=regime,
                        metric=metric,
                        transform=transform,
                        model_family="protocol_centered_descriptor",
                        tested_effect=f"workflow_x_protocol_centered_{base_label}",
                        descriptor=centered_label,
                        model_reduced="workflow_x_protocol",
                        model_full=f"workflow_x_protocol_centered_{base_label}",
                        work=work_c,
                        reduced=fit_reduced_c,
                        full=fit_full_c,
                        diag=diag_c,
                    )
                )

            # PCA sensitivity.
            work_pca, diag_pca = analysis_frame(
                df,
                metric=metric,
                transform=transform,
                eps=eps,
                instance_cols=instance_cols,
                required_columns=pca_required,
            )
            validate_analysis_frame(work_pca, regime=regime, metric=metric, transform=transform, model_name="workflow_x_morphology_pca")
            tools_pca = workflow_order(work_pca)
            y_pca = work_pca["_response"].to_numpy(dtype=float)
            inst_pca = work_pca["_instance_id"]
            wf_pca, wf_meta_pca = workflow_design(work_pca, tools_pca)
            pc_x, pc_meta = interaction_design(
                work_pca,
                tools_pca,
                pca_cols,
                term_type="workflow_x_morphology_pca",
                labels={"morphology_pc1": "PC1", "morphology_pc2": "PC2"},
            )
            full_pca_x, full_pca_meta = combine_designs((wf_pca, wf_meta_pca), (pc_x, pc_meta))
            fit_wf_pca = fit_ols(y_pca, wf_pca, inst_pca)
            fit_pca = fit_ols(y_pca, full_pca_x, inst_pca)
            fit_pca["term_meta"] = full_pca_meta
            comparison_rows.append(
                comparison_row(
                    regime=regime,
                    metric=metric,
                    transform=transform,
                    model_family="morphology_pca_sensitivity",
                    tested_effect="workflow_x_morphology_pca",
                    descriptor="PC1+PC2",
                    model_reduced="baseline_workflow",
                    model_full="workflow_x_morphology_pca",
                    work=work_pca,
                    reduced=fit_wf_pca,
                    full=fit_pca,
                    diag=diag_pca,
                )
            )
    diag_rows = descriptor_diagnostics(regime, df)
    for row in pca_diag:
        row["regime"] = regime
        diag_rows.append(row)
    return comparison_rows, coefficient_out, posthoc_out, diag_rows + instance_rows




def fit_lm(y: np.ndarray, x: np.ndarray) -> dict[str, Any]:
    n = len(y)
    if x.ndim == 1:
        x = x.reshape(-1, 1)
    design = np.column_stack([np.ones(n), x]) if x.shape[1] else np.ones((n, 1))
    rank = int(np.linalg.matrix_rank(design.T @ design))
    df_resid = int(n - rank)
    if n == 0 or rank == 0:
        return {"beta": np.zeros(design.shape[1]), "cov": np.zeros((design.shape[1], design.shape[1])), "rss": np.nan, "rank": rank, "df_resid": df_resid, "n_rows": n, "diagnostic_flag": "insufficient_rows"}
    xtx_inv = np.linalg.pinv(design.T @ design)
    beta = xtx_inv @ (design.T @ y)
    resid = y - design @ beta
    rss = float(np.sum(resid**2))
    sigma2 = rss / df_resid if df_resid > 0 else np.nan
    cov = sigma2 * xtx_inv if np.isfinite(sigma2) else np.full((design.shape[1], design.shape[1]), np.nan)
    flag = "ok" if df_resid > 0 else "non_positive_df_den"
    return {"beta": beta, "cov": cov, "rss": rss, "rank": rank, "df_resid": df_resid, "n_rows": n, "diagnostic_flag": flag}


def nested_lm_test(reduced: dict[str, Any], full: dict[str, Any]) -> dict[str, Any]:
    df_num = int(full["rank"] - reduced["rank"])
    df_den = int(full["df_resid"])
    out = {
        "rank_reduced": int(reduced["rank"]),
        "rank_full": int(full["rank"]),
        "df_num": df_num,
        "df_den": df_den,
        "rss_reduced": float(reduced["rss"]) if np.isfinite(reduced["rss"]) else np.nan,
        "rss_full": float(full["rss"]) if np.isfinite(full["rss"]) else np.nan,
        "r2": np.nan,
        "partial_r2": np.nan,
        "f_stat": np.nan,
        "p_value": np.nan,
        "diagnostic_flag": "ok",
    }
    if df_num <= 0:
        out["diagnostic_flag"] = "non_positive_df_num"
    elif df_den <= 0:
        out["diagnostic_flag"] = "non_positive_df_den"
    elif not np.isfinite(reduced["rss"]) or reduced["rss"] <= 0:
        out["diagnostic_flag"] = "rss_reduced_nonpositive"
    elif not np.isfinite(full["rss"]):
        out["diagnostic_flag"] = "rss_full_nonfinite"
    elif full["rss"] > reduced["rss"] + 1e-8:
        out["diagnostic_flag"] = "rss_full_greater_than_reduced"
    else:
        out["r2"] = (reduced["rss"] - full["rss"]) / reduced["rss"]
        out["partial_r2"] = out["r2"]
        out["f_stat"] = ((reduced["rss"] - full["rss"]) / df_num) / (full["rss"] / df_den)
        out["p_value"] = float(f_dist.sf(out["f_stat"], df_num, df_den)) if np.isfinite(out["f_stat"]) else np.nan
    return out


def protocol_design(df: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    protocols = sorted(df["protocol"].dropna().astype(str).unique())
    nonref = protocols[1:]
    if not nonref:
        return np.empty((len(df), 0)), protocols
    cols = [(df["protocol"].astype(str) == protocol).astype(float).to_numpy() for protocol in nonref]
    return np.column_stack(cols), protocols


def prepare_relative_frame(
    df: pd.DataFrame,
    *,
    regime: str,
    metric: str,
    transform: str,
    eps: float,
    instance_cols: list[str],
    required_columns: list[str] | None = None,
) -> pd.DataFrame:
    required = ["protocol", *DESCRIPTOR_Z]
    if required_columns:
        required.extend(required_columns)
    work, _ = analysis_frame(
        df,
        metric=metric,
        transform=transform,
        eps=eps,
        instance_cols=instance_cols,
        required_columns=list(dict.fromkeys(required)),
    )
    validate_analysis_frame(work, regime=regime, metric=metric, transform=transform, model_name="relative_response")
    grouped = work.groupby("_instance_id")["_response"]
    instance_sum = grouped.transform("sum")
    instance_count = grouped.transform("count")
    work = work.copy()
    work["_n_baseline_workflows"] = instance_count - 1
    work = work[work["_n_baseline_workflows"] >= 1].copy()
    work["_loo_instance_baseline"] = (instance_sum.loc[work.index] - work["_response"]) / work["_n_baseline_workflows"]
    work["_relative_response"] = work["_response"] - work["_loo_instance_baseline"]
    return work


def prepare_complete_case_relative_frame(
    df: pd.DataFrame,
    *,
    regime: str,
    metric: str,
    transform: str,
    eps: float,
    instance_cols: list[str],
    expected_workflows: list[str],
    required_columns: list[str] | None = None,
) -> pd.DataFrame:
    required = ["protocol", *DESCRIPTOR_Z]
    if required_columns:
        required.extend(required_columns)
    work, _ = analysis_frame(
        df,
        metric=metric,
        transform=transform,
        eps=eps,
        instance_cols=instance_cols,
        required_columns=list(dict.fromkeys(required)),
    )
    validate_analysis_frame(work, regime=regime, metric=metric, transform=transform, model_name="complete_case_relative_response")
    expected = set(expected_workflows)
    instance_tools = work.groupby("_instance_id")["tool"].agg(lambda values: set(values.astype(str)))
    keep_instances = instance_tools[instance_tools == expected].index
    work = work[work["_instance_id"].isin(keep_instances)].copy()
    if work.empty:
        return work
    grouped = work.groupby("_instance_id")["_response"]
    instance_sum = grouped.transform("sum")
    instance_count = grouped.transform("count")
    work["_n_baseline_workflows"] = instance_count - 1
    work = work[work["_n_baseline_workflows"] >= 1].copy()
    work["_loo_instance_baseline"] = (instance_sum.loc[work.index] - work["_response"]) / work["_n_baseline_workflows"]
    work["_relative_response"] = work["_response"] - work["_loo_instance_baseline"]
    return work


def sample_collapsed_relative_frame(rel: pd.DataFrame) -> pd.DataFrame:
    source_study_col = "tissue"
    needed = ["empirical_sample_id", "tool", "protocol", source_study_col, "_relative_response", *DESCRIPTOR_Z]
    missing = [col for col in needed if col not in rel.columns]
    if missing:
        raise KeyError(f"Missing columns for sample-collapsed sensitivity: {missing}")
    group_cols = ["empirical_sample_id", "tool"]
    agg_spec: dict[str, Any] = {
        "_relative_response": "mean",
        "protocol": "first",
        source_study_col: "first",
    }
    for col in DESCRIPTOR_Z:
        agg_spec[col] = "first"
    collapsed = rel.dropna(subset=needed).groupby(group_cols, as_index=False).agg(agg_spec)
    collapsed = collapsed.rename(columns={source_study_col: "study"})
    duplicated = collapsed.duplicated(subset=["empirical_sample_id", "tool"], keep=False)
    if duplicated.any():
        examples = collapsed.loc[duplicated, ["empirical_sample_id", "tool"]].head(10).to_dict("records")
        raise ValueError(f"Duplicate empirical_sample_id x workflow rows after sample collapse: {examples}")
    collapsed["_instance_id"] = collapsed["empirical_sample_id"].astype("string")
    return collapsed


def coverage_diagnostic_rows(
    *,
    regime: str,
    rel: pd.DataFrame,
    metric: str,
    transform: str,
    expected_workflows: list[str],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    expected = set(expected_workflows)
    instance_tools = rel.groupby("_instance_id")["tool"].agg(lambda values: set(values.astype(str)))
    k = instance_tools.map(len)
    complete_case = instance_tools.map(lambda tools: tools == expected)
    duplicated = rel.duplicated(subset=["_instance_id", "tool"], keep=False)
    sample_protocol = rel[["empirical_sample_id", "protocol"]].drop_duplicates()
    rows.append(
        {
            "regime": regime,
            "metric": metric,
            "response_transform": transform,
            "diagnostic_type": "overall_coverage",
            "n_rows": int(len(rel)),
            "n_instances": int(instance_tools.shape[0]),
            "n_empirical_samples": int(rel["empirical_sample_id"].nunique()),
            "n_protocols": int(rel["protocol"].nunique()),
            "n_expected_workflows": int(len(expected_workflows)),
            "k_min": int(k.min()) if not k.empty else 0,
            "k_median": float(k.median()) if not k.empty else np.nan,
            "k_max": int(k.max()) if not k.empty else 0,
            "n_complete_case_instances": int(complete_case.sum()) if not complete_case.empty else 0,
            "frac_complete_case_instances": float(complete_case.mean()) if not complete_case.empty else np.nan,
            "n_duplicate_instance_workflow_rows": int(duplicated.sum()),
        }
    )
    for protocol, sub in sample_protocol.groupby("protocol"):
        rows.append(
            {
                "regime": regime,
                "metric": metric,
                "response_transform": transform,
                "diagnostic_type": "samples_per_protocol",
                "protocol": protocol,
                "n_empirical_samples": int(sub["empirical_sample_id"].nunique()),
            }
        )
    for workflow in expected_workflows:
        has_workflow = instance_tools.map(lambda tools, workflow=workflow: workflow in tools)
        rows.append(
            {
                "regime": regime,
                "metric": metric,
                "response_transform": transform,
                "diagnostic_type": "workflow_missing_instances",
                "workflow": workflow,
                "n_instances_missing_workflow": int((~has_workflow).sum()) if not has_workflow.empty else 0,
                "frac_instances_missing_workflow": float((~has_workflow).mean()) if not has_workflow.empty else np.nan,
            }
        )
    return rows


def protocol_response_by_workflow_rows(regime: str, rel: pd.DataFrame, metric: str, transform: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    estimates: list[dict[str, Any]] = []
    for workflow in sorted(rel["tool"].dropna().astype(str).unique()):
        tool_df = rel[rel["tool"].astype(str) == workflow].dropna(subset=["_relative_response", "protocol"]).copy()
        n_rows = int(len(tool_df))
        n_instances = int(tool_df["_instance_id"].nunique()) if n_rows else 0
        n_protocols = int(tool_df["protocol"].dropna().astype(str).nunique()) if n_rows else 0
        y = tool_df["_relative_response"].to_numpy(dtype=float) if n_rows else np.array([])
        if n_rows < 3:
            test = {"rss_reduced": np.nan, "rss_full": np.nan, "r2": np.nan, "partial_r2": np.nan, "df_num": np.nan, "df_den": np.nan, "f_stat": np.nan, "p_value": np.nan, "diagnostic_flag": "insufficient_rows"}
        elif n_protocols < 2:
            reduced = fit_lm(y, np.empty((n_rows, 0)))
            test = {"rss_reduced": reduced["rss"], "rss_full": np.nan, "r2": np.nan, "partial_r2": np.nan, "df_num": 0, "df_den": reduced["df_resid"], "f_stat": np.nan, "p_value": np.nan, "diagnostic_flag": "fewer_than_2_protocols"}
        else:
            proto_x, _ = protocol_design(tool_df)
            reduced = fit_lm(y, np.empty((n_rows, 0)))
            full = fit_lm(y, proto_x)
            test = nested_lm_test(reduced, full)
        rows.append({
            "regime": regime,
            "metric": metric,
            "response_transform": transform,
            "is_primary": bool(transform == PRIMARY_TRANSFORM),
            "workflow": workflow,
            "model_family": "workflow_protocol_response",
            "n_rows": n_rows,
            "n_instances": n_instances,
            "n_protocols": n_protocols,
            **test,
            "fdr_group": "|".join([regime, metric, transform, "workflow_protocol_response"]),
        })
        for protocol, sub in tool_df.groupby("protocol"):
            values = sub["_relative_response"].to_numpy(dtype=float)
            se = float(np.std(values, ddof=1) / np.sqrt(len(values))) if len(values) > 1 else np.nan
            estimates.append({
                "regime": regime,
                "metric": metric,
                "response_transform": transform,
                "is_primary": bool(transform == PRIMARY_TRANSFORM),
                "workflow": workflow,
                "protocol": protocol,
                "mean_relative_response": float(np.mean(values)) if len(values) else np.nan,
                "standard_error": se,
                "n_instances": int(sub["_instance_id"].nunique()),
            })
    return rows, estimates


def descriptor_response_by_workflow_rows(regime: str, rel: pd.DataFrame, metric: str, transform: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for workflow in sorted(rel["tool"].dropna().astype(str).unique()):
        tool_df = rel[rel["tool"].astype(str) == workflow].copy()
        for desc_col in DESCRIPTOR_Z:
            desc_label = DESCRIPTOR_LABELS[desc_col]
            model_df = tool_df[["_instance_id", "_relative_response", desc_col]].dropna().copy()
            n_rows = int(len(model_df))
            n_instances = int(model_df["_instance_id"].nunique()) if n_rows else 0
            if n_rows < 3:
                fit = {"intercept": np.nan, "slope": np.nan, "standard_error": np.nan, "t_stat": np.nan, "p_value": np.nan, "rss_reduced": np.nan, "rss_full": np.nan, "r2": np.nan, "partial_r2": np.nan, "df_num": 1, "df_den": np.nan, "diagnostic_flag": "insufficient_rows"}
            else:
                fit = fit_internal_descriptor_model(
                    model_df["_relative_response"].to_numpy(dtype=float),
                    model_df[desc_col].to_numpy(dtype=float),
                )
            rows.append({
                "regime": regime,
                "metric": metric,
                "response_transform": transform,
                "is_primary": bool(transform == PRIMARY_TRANSFORM),
                "workflow": workflow,
                "descriptor": desc_label,
                "model_family": "workflow_descriptor_response",
                "n_rows": n_rows,
                "n_instances": n_instances,
                **fit,
                "r2_descriptor": fit.get("r2", np.nan),
                "fdr_group": "|".join([regime, metric, transform, "workflow_descriptor_response"]),
            })
    return rows


def study_loo_protocol_rows(
    regime: str,
    sample_rel: pd.DataFrame,
    metric: str,
    transform: str,
    sample_protocol: pd.DataFrame,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    idx = ["regime", "metric", "response_transform", "workflow"]
    source = sample_protocol.set_index(idx)
    studies = sorted(sample_rel["study"].dropna().astype(str).unique())
    for workflow in sorted(sample_rel["tool"].dropna().astype(str).unique()):
        base_key = (regime, metric, transform, workflow)
        base = source.loc[base_key].to_dict() if base_key in source.index else {}
        for study in studies:
            loo = sample_rel[(sample_rel["tool"].astype(str) == workflow) & (sample_rel["study"].astype(str) != study)].copy()
            protocol_rows, _ = protocol_response_by_workflow_rows(regime, loo, metric, transform)
            row = protocol_rows[0] if protocol_rows else {
                "n_rows": 0,
                "n_instances": 0,
                "n_protocols": 0,
                "r2": np.nan,
                "partial_r2": np.nan,
                "f_stat": np.nan,
                "p_value": np.nan,
                "diagnostic_flag": "insufficient_rows",
            }
            rows.append(
                {
                    "regime": regime,
                    "metric": metric,
                    "response_transform": transform,
                    "is_primary": bool(transform == PRIMARY_TRANSFORM),
                    "workflow": workflow,
                    "left_out_study": study,
                    "n_studies_total": int(len(studies)),
                    "n_rows": row.get("n_rows", 0),
                    "n_instances": row.get("n_instances", 0),
                    "n_protocols": row.get("n_protocols", 0),
                    "sample_collapsed_R2": base.get("r2", np.nan),
                    "sample_collapsed_FDR": base.get("fdr_bh", np.nan),
                    "loo_R2": row.get("r2", np.nan),
                    "loo_partial_R2": row.get("partial_r2", np.nan),
                    "loo_F_stat": row.get("f_stat", np.nan),
                    "loo_p_value": row.get("p_value", np.nan),
                    "loo_diagnostic_flag": row.get("diagnostic_flag", "insufficient_rows"),
                }
            )
    return rows, study_loo_protocol_summary(rows)


def study_loo_protocol_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not rows:
        return []
    df = pd.DataFrame(rows)
    out: list[dict[str, Any]] = []
    for keys, sub in df.groupby(["regime", "metric", "response_transform", "workflow"], dropna=False):
        ok = sub["loo_diagnostic_flag"].eq("ok") & sub["loo_R2"].notna()
        base_r2 = sub["sample_collapsed_R2"].dropna()
        base_fdr = sub["sample_collapsed_FDR"].dropna()
        r2 = sub.loc[ok, "loo_R2"]
        deltas = (r2 - base_r2.iloc[0]).abs() if not base_r2.empty else pd.Series(dtype=float)
        out.append(
            {
                "regime": keys[0],
                "metric": keys[1],
                "response_transform": keys[2],
                "is_primary": bool(keys[2] == PRIMARY_TRANSFORM),
                "workflow": keys[3],
                "sample_collapsed_R2": float(base_r2.iloc[0]) if not base_r2.empty else np.nan,
                "sample_collapsed_FDR": float(base_fdr.iloc[0]) if not base_fdr.empty else np.nan,
                "loo_R2_min": float(r2.min()) if not r2.empty else np.nan,
                "loo_R2_max": float(r2.max()) if not r2.empty else np.nan,
                "loo_R2_max_abs_delta": float(deltas.max()) if not deltas.empty else np.nan,
                "n_studies_tested": int(sub["left_out_study"].nunique()),
                "n_studies_ok": int(ok.sum()),
                "all_loo_significant": bool((sub.loc[ok, "loo_p_value"] < 0.05).all()) if ok.any() else False,
                "diagnostic_flag": "ok" if ok.all() else "some_leave_study_out_failed",
            }
        )
    return out


def study_loo_descriptor_rows(
    regime: str,
    sample_rel: pd.DataFrame,
    metric: str,
    transform: str,
    sample_descriptor: pd.DataFrame,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    idx = ["regime", "metric", "response_transform", "workflow", "descriptor"]
    source = sample_descriptor.set_index(idx)
    studies = sorted(sample_rel["study"].dropna().astype(str).unique())
    for workflow in sorted(sample_rel["tool"].dropna().astype(str).unique()):
        for desc_col in DESCRIPTOR_Z:
            desc_label = DESCRIPTOR_LABELS[desc_col]
            base_key = (regime, metric, transform, workflow, desc_label)
            base = source.loc[base_key].to_dict() if base_key in source.index else {}
            for study in studies:
                loo = sample_rel[
                    (sample_rel["tool"].astype(str) == workflow) & (sample_rel["study"].astype(str) != study)
                ].copy()
                descriptor_rows = descriptor_response_by_workflow_rows(regime, loo, metric, transform)
                row = next((item for item in descriptor_rows if item.get("workflow") == workflow and item.get("descriptor") == desc_label), None)
                if row is None:
                    row = {
                        "n_rows": 0,
                        "n_instances": 0,
                        "slope": np.nan,
                        "standard_error": np.nan,
                        "t_stat": np.nan,
                        "p_value": np.nan,
                        "r2_descriptor": np.nan,
                        "diagnostic_flag": "insufficient_rows",
                    }
                rows.append(
                    {
                        "regime": regime,
                        "metric": metric,
                        "response_transform": transform,
                        "is_primary": bool(transform == PRIMARY_TRANSFORM),
                        "workflow": workflow,
                        "descriptor": desc_label,
                        "left_out_study": study,
                        "n_studies_total": int(len(studies)),
                        "n_rows": row.get("n_rows", 0),
                        "n_instances": row.get("n_instances", 0),
                        "sample_collapsed_slope": base.get("slope", np.nan),
                        "sample_collapsed_R2_descriptor": base.get("r2_descriptor", np.nan),
                        "sample_collapsed_FDR": base.get("fdr_bh", np.nan),
                        "loo_slope": row.get("slope", np.nan),
                        "loo_R2_descriptor": row.get("r2_descriptor", np.nan),
                        "loo_p_value": row.get("p_value", np.nan),
                        "loo_diagnostic_flag": row.get("diagnostic_flag", "insufficient_rows"),
                    }
                )
    return rows, study_loo_descriptor_summary(rows)


def study_loo_descriptor_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not rows:
        return []
    df = pd.DataFrame(rows)
    out: list[dict[str, Any]] = []
    for keys, sub in df.groupby(["regime", "metric", "response_transform", "workflow", "descriptor"], dropna=False):
        ok = sub["loo_diagnostic_flag"].eq("ok") & sub["loo_R2_descriptor"].notna()
        base_slope = sub["sample_collapsed_slope"].dropna()
        base_r2 = sub["sample_collapsed_R2_descriptor"].dropna()
        base_fdr = sub["sample_collapsed_FDR"].dropna()
        slopes = sub.loc[ok, "loo_slope"].dropna()
        r2 = sub.loc[ok, "loo_R2_descriptor"]
        r2_delta = (r2 - base_r2.iloc[0]).abs() if not base_r2.empty else pd.Series(dtype=float)
        if not base_slope.empty and not slopes.empty:
            direction_agreement = float((np.sign(slopes) == np.sign(base_slope.iloc[0])).mean())
        else:
            direction_agreement = np.nan
        out.append(
            {
                "regime": keys[0],
                "metric": keys[1],
                "response_transform": keys[2],
                "is_primary": bool(keys[2] == PRIMARY_TRANSFORM),
                "workflow": keys[3],
                "descriptor": keys[4],
                "sample_collapsed_slope": float(base_slope.iloc[0]) if not base_slope.empty else np.nan,
                "sample_collapsed_R2_descriptor": float(base_r2.iloc[0]) if not base_r2.empty else np.nan,
                "sample_collapsed_FDR": float(base_fdr.iloc[0]) if not base_fdr.empty else np.nan,
                "loo_slope_min": float(slopes.min()) if not slopes.empty else np.nan,
                "loo_slope_max": float(slopes.max()) if not slopes.empty else np.nan,
                "loo_slope_direction_agreement_fraction": direction_agreement,
                "loo_R2_min": float(r2.min()) if not r2.empty else np.nan,
                "loo_R2_max": float(r2.max()) if not r2.empty else np.nan,
                "loo_R2_max_abs_delta": float(r2_delta.max()) if not r2_delta.empty else np.nan,
                "n_studies_tested": int(sub["left_out_study"].nunique()),
                "n_studies_ok": int(ok.sum()),
                "all_loo_significant": bool((sub.loc[ok, "loo_p_value"] < 0.05).all()) if ok.any() else False,
                "diagnostic_flag": "ok" if ok.all() else "some_leave_study_out_failed",
            }
        )
    return out


def descriptor_explains_protocol_rows(regime: str, rel: pd.DataFrame, metric: str, transform: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    centered_cols = [f"{col}_protocol_centered_z" for col in DESCRIPTOR_Z]
    for workflow in sorted(rel["tool"].dropna().astype(str).unique()):
        tool_df = rel[rel["tool"].astype(str) == workflow].dropna(subset=["_relative_response", "protocol", *DESCRIPTOR_Z, *centered_cols]).copy()
        n_rows = int(len(tool_df))
        n_instances = int(tool_df["_instance_id"].nunique()) if n_rows else 0
        n_protocols = int(tool_df["protocol"].dropna().astype(str).nunique()) if n_rows else 0
        base = {
            "regime": regime,
            "metric": metric,
            "response_transform": transform,
            "is_primary": bool(transform == PRIMARY_TRANSFORM),
            "workflow": workflow,
            "n_rows": n_rows,
            "n_instances": n_instances,
            "n_protocols": n_protocols,
        }
        if n_rows < 5 or n_protocols < 2:
            rows.append({**base, "r2_protocol": np.nan, "r2_all_descriptors": np.nan, "partial_r2_protocol_after_descriptors": np.nan, "partial_r2_within_protocol_descriptors": np.nan, "p_protocol": np.nan, "p_all_descriptors": np.nan, "p_protocol_after_descriptors": np.nan, "p_within_protocol_descriptors": np.nan, "diagnostic_flag": "insufficient_rows" if n_rows < 5 else "fewer_than_2_protocols"})
            continue
        y = tool_df["_relative_response"].to_numpy(dtype=float)
        intercept = fit_lm(y, np.empty((n_rows, 0)))
        proto_x, _ = protocol_design(tool_df)
        desc_x = tool_df[DESCRIPTOR_Z].to_numpy(dtype=float)
        centered_x = tool_df[centered_cols].to_numpy(dtype=float)
        proto = fit_lm(y, proto_x)
        desc = fit_lm(y, desc_x)
        desc_plus_proto = fit_lm(y, np.column_stack([desc_x, proto_x]))
        proto_plus_centered = fit_lm(y, np.column_stack([proto_x, centered_x]))
        protocol_test = nested_lm_test(intercept, proto)
        desc_test = nested_lm_test(intercept, desc)
        protocol_after_desc = nested_lm_test(desc, desc_plus_proto)
        centered_after_protocol = nested_lm_test(proto, proto_plus_centered)
        flag = "ok"
        for test in [protocol_test, desc_test, protocol_after_desc, centered_after_protocol]:
            if test["diagnostic_flag"] != "ok":
                flag = test["diagnostic_flag"]
                break
        rows.append({
            **base,
            "r2_protocol": protocol_test["r2"],
            "r2_all_descriptors": desc_test["r2"],
            "partial_r2_protocol_after_descriptors": protocol_after_desc["partial_r2"],
            "partial_r2_within_protocol_descriptors": centered_after_protocol["partial_r2"],
            "p_protocol": protocol_test["p_value"],
            "p_all_descriptors": desc_test["p_value"],
            "p_protocol_after_descriptors": protocol_after_desc["p_value"],
            "p_within_protocol_descriptors": centered_after_protocol["p_value"],
            "df_protocol_num": protocol_test["df_num"],
            "df_protocol_den": protocol_test["df_den"],
            "df_all_descriptors_num": desc_test["df_num"],
            "df_all_descriptors_den": desc_test["df_den"],
            "diagnostic_flag": flag,
        })
    return rows


def add_group_fdr(df: pd.DataFrame, *, p_col: str, out_col: str, group_cols: list[str]) -> pd.DataFrame:
    if df.empty:
        return df
    df[out_col] = np.nan
    for _, idx in df.groupby(group_cols).groups.items():
        df.loc[idx, out_col] = bh_fdr(df.loc[idx, p_col])
    return df


def transform_consistency_rows(protocol: pd.DataFrame, descriptor: pd.DataFrame) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not descriptor.empty:
        wide = descriptor.pivot_table(
            index=["regime", "metric", "workflow", "descriptor"],
            columns="response_transform",
            values=["slope", "r2_descriptor"],
            aggfunc="first",
        )
        if {"logit_clip", "identity"}.issubset(set(wide.columns.get_level_values(1))):
            flat = wide.reset_index()
            for (regime, metric), sub in flat.groupby([("regime", ""), ("metric", "")]):
                slope_logit = sub[("slope", "logit_clip")]
                slope_identity = sub[("slope", "identity")]
                r2_logit = sub[("r2_descriptor", "logit_clip")]
                r2_identity = sub[("r2_descriptor", "identity")]
                valid_slope = slope_logit.notna() & slope_identity.notna()
                valid_r2 = r2_logit.notna() & r2_identity.notna()
                rows.append(
                    {
                        "comparison_type": "descriptor_response",
                        "regime": regime,
                        "metric": metric,
                        "n_terms": int(valid_slope.sum()),
                        "direction_agreement_fraction": (
                            float((np.sign(slope_logit[valid_slope]) == np.sign(slope_identity[valid_slope])).mean())
                            if valid_slope.any()
                            else np.nan
                        ),
                        "slope_spearman": (
                            float(slope_logit[valid_slope].corr(slope_identity[valid_slope], method="spearman"))
                            if valid_slope.sum() >= 2
                            else np.nan
                        ),
                        "abs_slope_spearman": (
                            float(slope_logit[valid_slope].abs().corr(slope_identity[valid_slope].abs(), method="spearman"))
                            if valid_slope.sum() >= 2
                            else np.nan
                        ),
                        "r2_spearman": (
                            float(r2_logit[valid_r2].corr(r2_identity[valid_r2], method="spearman"))
                            if valid_r2.sum() >= 2
                            else np.nan
                        ),
                    }
                )
    if not protocol.empty:
        wide_p = protocol.pivot_table(
            index=["regime", "metric", "workflow"],
            columns="response_transform",
            values="r2",
            aggfunc="first",
        )
        if {"logit_clip", "identity"}.issubset(set(wide_p.columns)):
            flat_p = wide_p.reset_index()
            for (regime, metric), sub in flat_p.groupby(["regime", "metric"]):
                valid = sub["logit_clip"].notna() & sub["identity"].notna()
                rows.append(
                    {
                        "comparison_type": "protocol_response",
                        "regime": regime,
                        "metric": metric,
                        "n_terms": int(valid.sum()),
                        "direction_agreement_fraction": np.nan,
                        "slope_spearman": np.nan,
                        "abs_slope_spearman": np.nan,
                        "r2_spearman": (
                            float(sub.loc[valid, "logit_clip"].corr(sub.loc[valid, "identity"], method="spearman"))
                            if valid.sum() >= 2
                            else np.nan
                        ),
                    }
                )
    return rows


def add_fdr_tables(comparisons: pd.DataFrame, coefficients: pd.DataFrame, posthoc: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if not comparisons.empty:
        comparisons["fdr_bh"] = np.nan
        for _, idx in comparisons.groupby("fdr_group").groups.items():
            comparisons.loc[idx, "fdr_bh"] = bh_fdr(comparisons.loc[idx, "p_value"])
    if not coefficients.empty:
        coefficients["fdr_bh"] = np.nan
        for _, idx in coefficients.groupby(["regime", "response_transform", "metric", "model_family", "tested_effect"]).groups.items():
            coefficients.loc[idx, "fdr_bh"] = bh_fdr(coefficients.loc[idx, "p_value"])
    if not posthoc.empty:
        posthoc["fdr_bh"] = np.nan
        for _, idx in posthoc.groupby(["regime", "response_transform", "metric", "model_family", "tested_effect"]).groups.items():
            posthoc.loc[idx, "fdr_bh"] = bh_fdr(posthoc.loc[idx, "p_value"])
    return comparisons, coefficients, posthoc


def main() -> None:
    args = parse_args()
    out_dir = ensure_dir(args.output_dir)
    regime_paths = {
        "main_text_three_criteria": args.main_input,
        "dapars2_scmapa_compatible_single_filter": args.single_filter_input,
    }
    all_pooled: list[dict[str, Any]] = []
    all_protocol: list[dict[str, Any]] = []
    all_protocol_estimates: list[dict[str, Any]] = []
    all_descriptor: list[dict[str, Any]] = []
    all_explanation: list[dict[str, Any]] = []
    all_sample_collapsed_protocol: list[dict[str, Any]] = []
    all_sample_collapsed_descriptor: list[dict[str, Any]] = []
    all_complete_case_protocol: list[dict[str, Any]] = []
    all_complete_case_descriptor: list[dict[str, Any]] = []
    all_study_loo_protocol: list[dict[str, Any]] = []
    all_study_loo_descriptor: list[dict[str, Any]] = []
    all_study_loo_protocol_summary: list[dict[str, Any]] = []
    all_study_loo_descriptor_summary: list[dict[str, Any]] = []
    all_coverage: list[dict[str, Any]] = []
    all_diagnostics: list[dict[str, Any]] = []
    all_instances: list[dict[str, Any]] = []
    tool_universe: dict[str, list[str]] = {}
    instance_cols_by_regime: dict[str, list[str]] = {}

    for regime, cfg in REGIMES.items():
        path = regime_paths[regime]
        if not path.exists():
            raise FileNotFoundError(f"Missing descriptor performance table: {path}")
        df = pd.read_parquet(path)
        instance_cols = cfg["instance_cols"]
        instance_cols_by_regime[regime] = instance_cols
        tool_universe[regime] = sorted(df["tool"].dropna().astype(str).unique().tolist())

        # Secondary pooled/global heterogeneity tests.
        pooled, _coefs, _posthoc, diag_and_inst = run_regime_models(regime, df, instance_cols, args.logit_eps)
        all_pooled.extend(pooled)
        for row in diag_and_inst:
            if "median_tools_per_instance" in row:
                all_instances.append(row)
            else:
                all_diagnostics.append(row)

        # Primary workflow-specific analyses based on leave-one-out relative response.
        df_model = df.copy()
        df_model["_instance_id"] = make_instance_id(df_model, instance_cols)
        df_model = add_protocol_centered_descriptors(df_model)
        for metric in METRICS:
            for transform in TRANSFORMS:
                expected_workflows = tool_universe[regime]
                rel = prepare_relative_frame(
                    df_model,
                    regime=regime,
                    metric=metric,
                    transform=transform,
                    eps=args.logit_eps,
                    instance_cols=instance_cols,
                    required_columns=[f"{col}_protocol_centered_z" for col in DESCRIPTOR_Z],
                )
                all_coverage.extend(
                    coverage_diagnostic_rows(
                        regime=regime,
                        rel=rel,
                        metric=metric,
                        transform=transform,
                        expected_workflows=expected_workflows,
                    )
                )
                protocol_rows, estimate_rows = protocol_response_by_workflow_rows(regime, rel, metric, transform)
                all_protocol.extend(protocol_rows)
                all_protocol_estimates.extend(estimate_rows)
                all_descriptor.extend(descriptor_response_by_workflow_rows(regime, rel, metric, transform))
                all_explanation.extend(descriptor_explains_protocol_rows(regime, rel, metric, transform))

                sample_rel = sample_collapsed_relative_frame(rel)
                sample_protocol_rows, _ = protocol_response_by_workflow_rows(regime, sample_rel, metric, transform)
                sample_protocol_for_loo = add_group_fdr(
                    pd.DataFrame(sample_protocol_rows),
                    p_col="p_value",
                    out_col="fdr_bh",
                    group_cols=["regime", "metric", "response_transform"],
                )
                sample_descriptor_rows = descriptor_response_by_workflow_rows(regime, sample_rel, metric, transform)
                sample_descriptor_for_loo = add_group_fdr(
                    pd.DataFrame(sample_descriptor_rows),
                    p_col="p_value",
                    out_col="fdr_bh",
                    group_cols=["regime", "metric", "response_transform"],
                )
                all_sample_collapsed_protocol.extend(sample_protocol_rows)
                all_sample_collapsed_descriptor.extend(sample_descriptor_rows)
                study_protocol_rows, study_protocol_summary = study_loo_protocol_rows(
                    regime, sample_rel, metric, transform, sample_protocol_for_loo
                )
                study_descriptor_rows, study_descriptor_summary = study_loo_descriptor_rows(
                    regime, sample_rel, metric, transform, sample_descriptor_for_loo
                )
                all_study_loo_protocol.extend(study_protocol_rows)
                all_study_loo_descriptor.extend(study_descriptor_rows)
                all_study_loo_protocol_summary.extend(study_protocol_summary)
                all_study_loo_descriptor_summary.extend(study_descriptor_summary)

                complete_rel = prepare_complete_case_relative_frame(
                    df_model,
                    regime=regime,
                    metric=metric,
                    transform=transform,
                    eps=args.logit_eps,
                    instance_cols=instance_cols,
                    expected_workflows=expected_workflows,
                    required_columns=[f"{col}_protocol_centered_z" for col in DESCRIPTOR_Z],
                )
                if not complete_rel.empty:
                    complete_protocol_rows, _ = protocol_response_by_workflow_rows(regime, complete_rel, metric, transform)
                    all_complete_case_protocol.extend(complete_protocol_rows)
                    all_complete_case_descriptor.extend(descriptor_response_by_workflow_rows(regime, complete_rel, metric, transform))
                else:
                    for workflow in expected_workflows:
                        all_complete_case_protocol.append(
                            {
                                "regime": regime,
                                "metric": metric,
                                "response_transform": transform,
                                "is_primary": bool(transform == PRIMARY_TRANSFORM),
                                "workflow": workflow,
                                "model_family": "workflow_protocol_response",
                                "n_rows": 0,
                                "n_instances": 0,
                                "n_protocols": 0,
                                "rss_reduced": np.nan,
                                "rss_full": np.nan,
                                "r2": np.nan,
                                "partial_r2": np.nan,
                                "df_num": np.nan,
                                "df_den": np.nan,
                                "f_stat": np.nan,
                                "p_value": np.nan,
                                "diagnostic_flag": "no_complete_case_instances",
                                "fdr_group": "|".join([regime, metric, transform, "workflow_protocol_response"]),
                            }
                        )
                        for desc_col in DESCRIPTOR_Z:
                            all_complete_case_descriptor.append(
                                {
                                    "regime": regime,
                                    "metric": metric,
                                    "response_transform": transform,
                                    "is_primary": bool(transform == PRIMARY_TRANSFORM),
                                    "workflow": workflow,
                                    "descriptor": DESCRIPTOR_LABELS[desc_col],
                                    "model_family": "workflow_descriptor_response",
                                    "n_rows": 0,
                                    "n_instances": 0,
                                    "intercept": np.nan,
                                    "slope": np.nan,
                                    "standard_error": np.nan,
                                    "t_stat": np.nan,
                                    "p_value": np.nan,
                                    "rss_reduced": np.nan,
                                    "rss_full": np.nan,
                                    "r2": np.nan,
                                    "partial_r2": np.nan,
                                    "df_num": 1,
                                    "df_den": np.nan,
                                    "diagnostic_flag": "no_complete_case_instances",
                                    "r2_descriptor": np.nan,
                                    "fdr_group": "|".join([regime, metric, transform, "workflow_descriptor_response"]),
                                }
                            )

    pooled = pd.DataFrame(all_pooled)
    protocol = pd.DataFrame(all_protocol)
    protocol_estimates = pd.DataFrame(all_protocol_estimates)
    descriptor = pd.DataFrame(all_descriptor)
    explanation = pd.DataFrame(all_explanation)
    sample_collapsed_protocol = pd.DataFrame(all_sample_collapsed_protocol)
    sample_collapsed_descriptor = pd.DataFrame(all_sample_collapsed_descriptor)
    complete_case_protocol = pd.DataFrame(all_complete_case_protocol)
    complete_case_descriptor = pd.DataFrame(all_complete_case_descriptor)
    study_loo_protocol = pd.DataFrame(all_study_loo_protocol)
    study_loo_descriptor = pd.DataFrame(all_study_loo_descriptor)
    study_loo_protocol_summary = pd.DataFrame(all_study_loo_protocol_summary)
    study_loo_descriptor_summary = pd.DataFrame(all_study_loo_descriptor_summary)
    coverage = pd.DataFrame(all_coverage)

    pooled, _, _ = add_fdr_tables(pooled, pd.DataFrame(), pd.DataFrame())
    protocol = add_group_fdr(protocol, p_col="p_value", out_col="fdr_bh", group_cols=["regime", "metric", "response_transform"])
    descriptor = add_group_fdr(descriptor, p_col="p_value", out_col="fdr_bh", group_cols=["regime", "metric", "response_transform"])
    sample_collapsed_protocol = add_group_fdr(sample_collapsed_protocol, p_col="p_value", out_col="fdr_bh", group_cols=["regime", "metric", "response_transform"])
    sample_collapsed_descriptor = add_group_fdr(sample_collapsed_descriptor, p_col="p_value", out_col="fdr_bh", group_cols=["regime", "metric", "response_transform"])
    complete_case_protocol = add_group_fdr(complete_case_protocol, p_col="p_value", out_col="fdr_bh", group_cols=["regime", "metric", "response_transform"])
    complete_case_descriptor = add_group_fdr(complete_case_descriptor, p_col="p_value", out_col="fdr_bh", group_cols=["regime", "metric", "response_transform"])
    for p_col, fdr_col in [
        ("p_protocol", "fdr_protocol"),
        ("p_all_descriptors", "fdr_all_descriptors"),
        ("p_protocol_after_descriptors", "fdr_protocol_after_descriptors"),
        ("p_within_protocol_descriptors", "fdr_within_protocol_descriptors"),
    ]:
        explanation = add_group_fdr(explanation, p_col=p_col, out_col=fdr_col, group_cols=["regime", "metric", "response_transform"])
    consistency = pd.DataFrame(transform_consistency_rows(protocol, descriptor))

    write_tsv(protocol, out_dir / "protocol_response_by_workflow.tsv")
    write_tsv(protocol_estimates, out_dir / "protocol_response_estimates.tsv")
    write_tsv(descriptor, out_dir / "descriptor_response_by_workflow.tsv")
    write_tsv(explanation, out_dir / "descriptor_explains_protocol_response.tsv")
    write_tsv(sample_collapsed_protocol, out_dir / "sample_collapsed_protocol_response_by_workflow.tsv")
    write_tsv(sample_collapsed_descriptor, out_dir / "sample_collapsed_descriptor_response_by_workflow.tsv")
    write_tsv(complete_case_protocol, out_dir / "complete_case_protocol_response_by_workflow.tsv")
    write_tsv(complete_case_descriptor, out_dir / "complete_case_descriptor_response_by_workflow.tsv")
    write_tsv(study_loo_protocol, out_dir / "study_loo_protocol_response_by_workflow.tsv")
    write_tsv(study_loo_descriptor, out_dir / "study_loo_descriptor_response_by_workflow.tsv")
    write_tsv(study_loo_protocol_summary, out_dir / "study_loo_protocol_sensitivity_summary.tsv")
    write_tsv(study_loo_descriptor_summary, out_dir / "study_loo_descriptor_sensitivity_summary.tsv")
    write_tsv(coverage, out_dir / "morphology_coverage_diagnostics.tsv")
    write_tsv(consistency, out_dir / "transform_consistency_summary.tsv")
    write_tsv(pooled, out_dir / "pooled_protocol_descriptor_model_comparison.tsv")
    write_tsv(pd.DataFrame(all_diagnostics), out_dir / "morphology_descriptor_diagnostics.tsv")
    write_tsv(pd.DataFrame(all_instances), out_dir / "morphology_instance_summary.tsv")
    write_json(
        {
            "instance_id_cols_by_regime": instance_cols_by_regime,
            "regime_tool_universe": tool_universe,
            "regime_notes": {regime: cfg["note"] for regime, cfg in REGIMES.items()},
            "response_transforms": TRANSFORMS,
            "primary_transform": PRIMARY_TRANSFORM,
            "logit_eps": float(args.logit_eps),
            "descriptor_columns": DESCRIPTORS,
            "descriptor_z_columns": DESCRIPTOR_Z,
            "primary_protocol_output": "protocol_response_by_workflow.tsv",
            "primary_descriptor_output": "descriptor_response_by_workflow.tsv",
            "protocol_estimates_output": "protocol_response_estimates.tsv",
            "descriptor_protocol_explanation_output": "descriptor_explains_protocol_response.tsv",
            "sample_collapsed_protocol_output": "sample_collapsed_protocol_response_by_workflow.tsv",
            "sample_collapsed_descriptor_output": "sample_collapsed_descriptor_response_by_workflow.tsv",
            "complete_case_protocol_output": "complete_case_protocol_response_by_workflow.tsv",
            "complete_case_descriptor_output": "complete_case_descriptor_response_by_workflow.tsv",
            "study_loo_protocol_output": "study_loo_protocol_response_by_workflow.tsv",
            "study_loo_descriptor_output": "study_loo_descriptor_response_by_workflow.tsv",
            "study_loo_protocol_summary_output": "study_loo_protocol_sensitivity_summary.tsv",
            "study_loo_descriptor_summary_output": "study_loo_descriptor_sensitivity_summary.tsv",
            "coverage_diagnostics_output": "morphology_coverage_diagnostics.tsv",
            "transform_consistency_output": "transform_consistency_summary.tsv",
            "secondary_global_output": "pooled_protocol_descriptor_model_comparison.tsv",
            "relative_response_definition": (
                "For benchmark instance i and workflow m, R_i,m = g(Y_i,m) - "
                "mean_{m' != m} g(Y_i,m'), using workflows with valid measurements in the same instance."
            ),
            "fdr_protocol_response": "BH within regime x metric x response_transform across workflows",
            "fdr_descriptor_response": "BH within regime x metric x response_transform across workflow x descriptor tests",
            "fdr_descriptor_explanation": "Separate BH within regime x metric x response_transform for each p-value family",
            "sensitivity_note": (
                "Sample-collapsed sensitivity averages relative response within empirical_sample_id x workflow "
                "before refitting workflow-specific protocol and descriptor models. Complete-case sensitivity "
                "keeps only instances containing every workflow in the regime tool universe and recomputes the "
                "leave-one-out baseline within that subset. Study leave-one-out sensitivity derives study groups "
                "from the available source grouping column and refits sample-collapsed models after removing "
                "one study at a time."
            ),
            "model_note": (
                "Workflow-specific protocol and descriptor responses are fit separately for each workflow "
                "using leave-one-out instance-adjusted relative performance. No workflow-specific output "
                "uses scAPAtrap as a reference. Pooled workflow interaction models are secondary global "
                "heterogeneity tests only."
            ),
            "non_additivity_note": (
                "Protocol and descriptor explanatory fractions are descriptive and should not be interpreted "
                "as additive or causal because protocol labels and peak descriptors may be correlated."
            ),
        },
        out_dir / "morphology_model_metadata.json",
    )
    print(f"Wrote workflow response model outputs to: {out_dir}")


if __name__ == "__main__":
    main()
