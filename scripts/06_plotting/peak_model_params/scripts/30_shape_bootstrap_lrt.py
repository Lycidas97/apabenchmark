#!/usr/bin/env python3
"""Run bootstrap LRT on shape core parameters.

Statistical design follows reference notebook `peak_model_lrt.ipynb`:
- Fixed effects: species, protocol, and species:protocol interaction
- Random effect: tissue (random intercept)
- Bootstrap-based LRT p-values

This script applies the same tests to shape core features:
- all numeric columns prefixed with `shape_per_peak_`
- all numeric columns prefixed with `shape_pooled_`
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from tqdm import tqdm
import warnings

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.io import read_table, write_json, write_tsv
from _shared.paths import topic_intermediate_dir

TOPIC = "peak_model_params"
INPUT_NAME = "peak_model_params_prepared.tsv"
OUTPUT_NAME = "shape_bootstrap_lrt_results.tsv"
METADATA_NAME = "shape_bootstrap_lrt_metadata.json"

SHAPE_PREFIXES = ("shape_per_peak_", "shape_pooled_")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run bootstrap LRT tests on shape core parameters")
    parser.add_argument(
        "--input",
        type=Path,
        default=topic_intermediate_dir(TOPIC) / INPUT_NAME,
        help="Prepared input table",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=topic_intermediate_dir(TOPIC) / OUTPUT_NAME,
        help="Output results TSV",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=topic_intermediate_dir(TOPIC) / METADATA_NAME,
        help="Output metadata JSON",
    )
    parser.add_argument(
        "--b",
        type=int,
        default=5000,
        help="Bootstrap iterations (default: 5000)",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument(
        "--max-metrics",
        type=int,
        default=None,
        help="Optional cap for number of metrics (useful for smoke tests)",
    )
    parser.add_argument(
        "--min-obs",
        type=int,
        default=20,
        help="Minimum non-null observations required per metric",
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=100,
        help="Max retries per bootstrap iteration for mixed model bootstrap",
    )
    return parser.parse_args()


def select_shape_metrics(df: pd.DataFrame) -> list[str]:
    metrics = []
    for col in df.columns:
        if not col.startswith(SHAPE_PREFIXES):
            continue
        if pd.api.types.is_numeric_dtype(df[col]):
            metrics.append(col)
    return sorted(metrics)


def build_analysis_df(df: pd.DataFrame) -> pd.DataFrame:
    required = {"sample_id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in input: {sorted(missing)}")

    work = df.copy()
    parts = work["sample_id"].astype(str).str.split("_", expand=True)
    if parts.shape[1] < 3:
        raise ValueError("sample_id format is unexpected; need protocol/species/tissue in first 3 fields.")

    work["protocol_raw"] = parts.iloc[:, 0]
    work["species"] = parts.iloc[:, 1]
    work["tissue"] = parts.iloc[:, 2]

    # Match reference notebook rule: keep protocols with two species represented.
    protocol_counts = work.groupby("protocol_raw")["species"].nunique()
    valid_protocols = protocol_counts[protocol_counts == 2].index
    filtered = work[work["protocol_raw"].isin(valid_protocols)].copy()
    filtered["protocol"] = filtered["protocol_raw"].astype(str)
    filtered["species"] = filtered["species"].astype(str)
    filtered["tissue"] = filtered["tissue"].astype(str)
    return filtered


def bootstrap_lrt(
    *,
    null_formula: str,
    full_formula: str,
    group_var: str,
    re_formula: str,
    data: pd.DataFrame,
    b: int,
    rng: np.random.Generator,
    max_retries: int,
) -> tuple[float, float, int, bool]:
    null_model = smf.mixedlm(null_formula, data=data, groups=data[group_var], re_formula=re_formula)
    null_result = null_model.fit(reml=False, maxiter=2000)

    full_model = smf.mixedlm(full_formula, data=data, groups=data[group_var], re_formula=re_formula)
    full_result = full_model.fit(reml=False, maxiter=2000)

    observed_converged = bool(
        null_result.converged
        and full_result.converged
        and np.isfinite(null_result.llf)
        and np.isfinite(full_result.llf)
    )
    lrt_observed = max(-2 * (null_result.llf - full_result.llf), 0.0)

    x_null = null_model.exog
    fixef = null_result.fe_params
    ranef_var = (
        float(null_result.cov_re.values[0][0])
        if getattr(null_result, "cov_re", None) is not None
        else 0.0
    )
    resid_var = float(null_result.scale)
    ranef_var = max(ranef_var, 0.0)
    resid_var = max(resid_var, 1e-12)

    lrt_bootstrap: list[float] = []
    groups = data[group_var]
    unique_groups = groups.unique().tolist()

    for _ in tqdm(range(b), desc="Bootstrapping mixed LRT", leave=False):
        success = False
        retries = 0
        while not success and retries < max_retries:
            try:
                group_effects = {
                    grp: rng.normal(0.0, np.sqrt(ranef_var)) for grp in unique_groups
                }
                re = groups.map(group_effects).to_numpy()
                y = x_null.dot(fixef) + re + rng.normal(0.0, np.sqrt(resid_var), len(data))

                boot_df = data.copy()
                boot_df["value"] = y

                null_boot = smf.mixedlm(
                    null_formula,
                    data=boot_df,
                    groups=boot_df[group_var],
                    re_formula=re_formula,
                )
                null_boot_result = null_boot.fit(reml=False, maxiter=10000, method="powell")

                full_boot = smf.mixedlm(
                    full_formula,
                    data=boot_df,
                    groups=boot_df[group_var],
                    re_formula=re_formula,
                )
                full_boot_result = full_boot.fit(reml=False, maxiter=10000, method="powell")

                lrt = max(-2 * (null_boot_result.llf - full_boot_result.llf), 0.0)
                if null_boot_result.converged and full_boot_result.converged and lrt >= 0:
                    lrt_bootstrap.append(lrt)
                    success = True
                else:
                    retries += 1
            except Exception:
                retries += 1

    if not lrt_bootstrap:
        return lrt_observed, float("nan"), 0, observed_converged

    arr = np.asarray(lrt_bootstrap)
    pval = (
        float("nan")
        if not np.isfinite(lrt_observed)
        else (np.sum(arr >= lrt_observed) + 1) / (len(arr) + 1)
    )
    return lrt_observed, float(pval), int(len(arr)), observed_converged


def bootstrap_lrt_random(
    *,
    data: pd.DataFrame,
    group_var: str,
    b: int,
    rng: np.random.Generator,
    max_retries: int,
) -> tuple[float, float, int, bool]:
    null_formula = "value ~ species * protocol"
    null_model = smf.ols(null_formula, data=data)
    null_result = null_model.fit()

    full_formula = "value ~ species * protocol"
    full_model = smf.mixedlm(
        full_formula,
        groups=data[group_var],
        re_formula="1",
        data=data,
    )
    full_result = full_model.fit(reml=False, method="lbfgs", maxiter=3000)

    observed_converged = bool(
        full_result.converged and np.isfinite(null_result.llf) and np.isfinite(full_result.llf)
    )
    lrt_observed = max(-2 * (null_result.llf - full_result.llf), 0.0)

    x = null_model.exog
    fixef = null_result.params
    resid_var = max(float(null_result.mse_resid), 1e-6)

    lrt_bootstrap: list[float] = []
    for _ in tqdm(range(b), desc="Bootstrapping tissue RE LRT", leave=False):
        success = False
        retries = 0
        while not success and retries < max_retries:
            try:
                y = x.dot(fixef) + rng.normal(0.0, np.sqrt(resid_var), len(x))
                boot_df = data.copy()
                boot_df["value"] = y

                null_boot = smf.ols(null_formula, data=boot_df)
                null_boot_result = null_boot.fit()

                full_boot = smf.mixedlm(
                    full_formula,
                    groups=boot_df[group_var],
                    re_formula="1",
                    data=boot_df,
                )
                full_boot_result = full_boot.fit(reml=False, method="nm", maxiter=3000)

                lrt = max(-2 * (null_boot_result.llf - full_boot_result.llf), 0.0)
                if full_boot_result.converged and lrt >= 0:
                    lrt_bootstrap.append(lrt)
                    success = True
                else:
                    retries += 1
            except Exception:
                retries += 1

    if not lrt_bootstrap:
        return lrt_observed, float("nan"), 0, observed_converged

    arr = np.asarray(lrt_bootstrap)
    pval = (
        float("nan")
        if not np.isfinite(lrt_observed)
        else (np.sum(arr >= lrt_observed) + 1) / (len(arr) + 1)
    )
    return lrt_observed, float(pval), int(len(arr)), observed_converged


def add_bh_fdr(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for effect in ["Species", "Protocol", "Interaction", "Tissue"]:
        pcol = f"{effect}_p"
        qcol = f"{effect}_q"
        out[qcol] = np.nan
        if pcol not in out.columns:
            continue
        pvals = pd.to_numeric(out[pcol], errors="coerce")
        mask = np.isfinite(pvals)
        if mask.sum() == 0:
            continue
        _, qvals, _, _ = multipletests(pvals[mask], method="fdr_bh")
        out.loc[mask, qcol] = qvals
    return out


def run_metric_tests(
    *,
    metric: str,
    analysis_df: pd.DataFrame,
    b: int,
    min_obs: int,
    max_retries: int,
    rng: np.random.Generator,
) -> dict[str, Any]:
    current = analysis_df[["protocol", "species", "tissue", metric]].copy()
    current["value"] = pd.to_numeric(current[metric], errors="coerce")
    current = current.dropna(subset=["value"])

    row: dict[str, Any] = {
        "Metric": metric,
        "n_obs": int(len(current)),
        "n_protocol": int(current["protocol"].nunique()),
        "n_species": int(current["species"].nunique()),
        "n_tissue": int(current["tissue"].nunique()),
        "B_target": int(b),
        "Observed_converged": None,
        "Error": None,
    }

    if len(current) < min_obs:
        row["Error"] = f"insufficient observations (<{min_obs})"
        return row
    if current["species"].nunique() < 2:
        row["Error"] = "species has <2 levels"
        return row
    if current["protocol"].nunique() < 2:
        row["Error"] = "protocol has <2 levels"
        return row
    if current["tissue"].nunique() < 2:
        row["Error"] = "tissue has <2 levels"
        return row

    try:
        species_lrt, species_p, species_boot, species_converged = bootstrap_lrt(
            null_formula="value ~ protocol",
            full_formula="value ~ species * protocol",
            group_var="tissue",
            re_formula="1",
            data=current,
            b=b,
            rng=rng,
            max_retries=max_retries,
        )
        protocol_lrt, protocol_p, protocol_boot, protocol_converged = bootstrap_lrt(
            null_formula="value ~ species",
            full_formula="value ~ species * protocol",
            group_var="tissue",
            re_formula="1",
            data=current,
            b=b,
            rng=rng,
            max_retries=max_retries,
        )
        interaction_lrt, interaction_p, interaction_boot, interaction_converged = bootstrap_lrt(
            null_formula="value ~ species + protocol",
            full_formula="value ~ species * protocol",
            group_var="tissue",
            re_formula="1",
            data=current,
            b=b,
            rng=rng,
            max_retries=max_retries,
        )
        tissue_lrt, tissue_p, tissue_boot, tissue_converged = bootstrap_lrt_random(
            data=current,
            group_var="tissue",
            b=b,
            rng=rng,
            max_retries=max_retries,
        )

        row.update(
            {
                "Species_LRT": species_lrt,
                "Species_p": species_p,
                "Species_boot_n": species_boot,
                "Protocol_LRT": protocol_lrt,
                "Protocol_p": protocol_p,
                "Protocol_boot_n": protocol_boot,
                "Interaction_LRT": interaction_lrt,
                "Interaction_p": interaction_p,
                "Interaction_boot_n": interaction_boot,
                "Tissue_LRT": tissue_lrt,
                "Tissue_p": tissue_p,
                "Tissue_boot_n": tissue_boot,
                "Observed_converged": bool(
                    species_converged
                    and protocol_converged
                    and interaction_converged
                    and tissue_converged
                ),
            }
        )
    except Exception as exc:
        row["Error"] = str(exc)
    return row


def main() -> None:
    warnings.simplefilter("ignore", ConvergenceWarning)
    args = parse_args()
    start = time.time()

    df = read_table(args.input)
    analysis_df = build_analysis_df(df)

    metrics = select_shape_metrics(analysis_df)
    if args.max_metrics is not None:
        metrics = metrics[: args.max_metrics]
    if not metrics:
        raise ValueError("No shape core metrics selected from input table.")

    rng = np.random.default_rng(args.seed)
    results = []
    for metric in metrics:
        print(f"Running bootstrap LRT for metric: {metric}")
        row = run_metric_tests(
            metric=metric,
            analysis_df=analysis_df,
            b=args.b,
            min_obs=args.min_obs,
            max_retries=args.max_retries,
            rng=rng,
        )
        results.append(row)

    results_df = pd.DataFrame(results)
    results_df = add_bh_fdr(results_df)
    results_df = results_df.sort_values("Metric", kind="stable").reset_index(drop=True)

    write_tsv(results_df, args.output)
    elapsed = time.time() - start
    write_json(
        {
            "topic": TOPIC,
            "input": str(args.input),
            "output": str(args.output),
            "rows": int(len(results_df)),
            "b": int(args.b),
            "seed": int(args.seed),
            "min_obs": int(args.min_obs),
            "max_retries": int(args.max_retries),
            "analysis_rows": int(len(analysis_df)),
            "analysis_protocols": int(analysis_df["protocol"].nunique()),
            "analysis_species": int(analysis_df["species"].nunique()),
            "analysis_tissues": int(analysis_df["tissue"].nunique()),
            "shape_metric_count": int(len(metrics)),
            "error_count": int(results_df["Error"].notna().sum()) if "Error" in results_df else 0,
            "elapsed_seconds": elapsed,
        },
        args.metadata,
    )

    print(f"Wrote bootstrap LRT results: {args.output}")
    print(f"Wrote metadata: {args.metadata}")


if __name__ == "__main__":
    main()
