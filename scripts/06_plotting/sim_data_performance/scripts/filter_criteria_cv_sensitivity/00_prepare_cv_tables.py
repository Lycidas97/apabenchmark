#!/usr/bin/env python3
"""Prepare post-hoc criterion CV sensitivity tables for DE-APA detection."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[2]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))
TOPIC_SCRIPT_DIR = SCRIPT_DIR.parent
if str(TOPIC_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(TOPIC_SCRIPT_DIR))

from _shared.style import apply_reference_style  # noqa: E402
from _plot_helpers import _pick_apa_filter_set  # noqa: E402

TOPIC_DIR = SCRIPT_DIR.parents[1]
INPUT_TABLE = TOPIC_DIR / "data" / "intermediate" / "sim_data_apa_residuals.parquet"
OUTPUT_DIR = SCRIPT_DIR / "output"
ROW_MATRIX_TABLE = OUTPUT_DIR / "apa_criterion_performance_row_matrix.parquet"
ROW_MATRIX_METADATA = OUTPUT_DIR / "apa_criterion_performance_row_matrix.metadata.json"

REQUIRED_COLUMNS = [
    "sample",
    "tool",
    "protocol",
    "match_type",
    "filter_type_1",
    "filter_type_2",
    "filter_type_comb",
    "precision",
    "recall",
    "f1",
    "gt_recall",
    "f1_residuals",
]

SAMPLE_SUFFIX_RE = re.compile(
    r"^(?P<background>.+?)_(?P<genome>hg38|mm10)_pas(?P<pas_set>[0-9]+)_gn(?P<gn>[0-9]+)_rep(?P<rep>[0-9]+)$"
)

TIER_ORDER = ["strict", "moderate", "lenient"]
PRIMARY_METRIC = "f1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT_TABLE)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--n-random-splits", type=int, default=100)
    parser.add_argument("--random-seed", type=int, default=1)
    parser.add_argument("--train-frac", type=float, default=0.5)
    parser.add_argument("--split-unit", choices=["background", "sample"], default="background")
    parser.add_argument("--f1-lower-quantile", type=float, default=0.01)
    parser.add_argument("--f1-upper-quantile", type=float, default=0.99)
    parser.add_argument("--rebuild-row-matrix", action="store_true")
    parser.add_argument("--skip-leave-protocol-out", action="store_true")
    return parser.parse_args()


def empirical_background(sample: str) -> str:
    text = str(sample)
    match = SAMPLE_SUFFIX_RE.match(text)
    if match:
        return match.group("background")
    return text


def criterion_family(filter_type_1: str, filter_type_2: str) -> str:
    f1 = str(filter_type_1)
    f2 = str(filter_type_2)
    if f1.startswith("dexseq"):
        test = "DEXSeq"
    elif f1.startswith("fisher"):
        test = "Fisher"
    elif f1.startswith("wilcox_"):
        parts = f1.split("_")
        test = f"Wilcoxon-{parts[1]}" if len(parts) > 1 else "Wilcoxon"
    else:
        test = f1.split("_")[0]

    if f2.startswith("dexseq_log2fc"):
        effect = "log2FC"
    else:
        effect = f2.split("_")[0]
    return f"{test}-{effect}"


def row_matrix_path(output_dir: Path) -> Path:
    return output_dir / ROW_MATRIX_TABLE.name


def row_matrix_metadata_path(output_dir: Path) -> Path:
    return output_dir / ROW_MATRIX_METADATA.name


def cache_matches(output_dir: Path, *, lower_q: float, upper_q: float) -> bool:
    meta_path = row_matrix_metadata_path(output_dir)
    if not meta_path.exists():
        return False
    try:
        metadata = json.loads(meta_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return False
    return (
        metadata.get("residual_method") == "per_sample_tool_match_linear_f1_on_gt_recall"
        and np.isclose(float(metadata.get("f1_lower_quantile", np.nan)), float(lower_q))
        and np.isclose(float(metadata.get("f1_upper_quantile", np.nan)), float(upper_q))
    )


def add_group_residuals_once(
    df: pd.DataFrame,
    *,
    lower_q: float,
    upper_q: float,
) -> pd.DataFrame:
    """Precompute split-invariant residuals within sample/tool/match groups."""
    out = df.copy()
    group_id = out.groupby(["sample", "tool", "match_type"], dropna=False, sort=False).ngroup()

    f1 = pd.to_numeric(out["f1"], errors="coerce")
    gt_recall = pd.to_numeric(out["gt_recall"], errors="coerce")
    quantiles = f1.groupby(group_id).quantile([lower_q, upper_q]).unstack()
    lo = group_id.map(quantiles[lower_q])
    hi = group_id.map(quantiles[upper_q])
    keep = f1.notna() & gt_recall.notna() & f1.between(lo, hi, inclusive="both")

    work = pd.DataFrame(
        {
            "_row": np.flatnonzero(keep.to_numpy()),
            "_group_id": group_id[keep].to_numpy(dtype=np.int64),
            "_x": gt_recall[keep].to_numpy(dtype=float),
            "_y": f1[keep].to_numpy(dtype=float),
        }
    )
    if work.empty:
        out["f1_residual_cv"] = np.nan
        return out

    work["_xx"] = work["_x"] * work["_x"]
    work["_xy"] = work["_x"] * work["_y"]
    stats = (
        work.groupby("_group_id", sort=False)
        .agg(n=("_y", "size"), sx=("_x", "sum"), sy=("_y", "sum"), sxx=("_xx", "sum"), sxy=("_xy", "sum"))
    )
    denom = stats["n"] * stats["sxx"] - stats["sx"] * stats["sx"]
    slope = (stats["n"] * stats["sxy"] - stats["sx"] * stats["sy"]) / denom.replace(0, np.nan)
    intercept = (stats["sy"] - slope * stats["sx"]) / stats["n"]

    # Rank-deficient groups have constant stringency; residuals reduce to y - mean(y).
    constant = slope.isna()
    slope.loc[constant] = 0.0
    intercept.loc[constant] = stats.loc[constant, "sy"] / stats.loc[constant, "n"]

    coef = pd.DataFrame({"intercept": intercept, "slope": slope})
    work = work.join(coef, on="_group_id")
    residual = work["_y"].to_numpy(dtype=float) - (
        work["intercept"].to_numpy(dtype=float) + work["slope"].to_numpy(dtype=float) * work["_x"].to_numpy(dtype=float)
    )
    values = np.full(len(out), np.nan, dtype=float)
    values[work["_row"].to_numpy(dtype=np.int64)] = residual
    out["f1_residual_cv"] = values
    return out


def load_or_build_matrix(
    input_path: Path,
    output_dir: Path,
    rebuild: bool,
    *,
    lower_q: float,
    upper_q: float,
) -> pd.DataFrame:
    output_dir.mkdir(parents=True, exist_ok=True)
    cached_path = row_matrix_path(output_dir)
    if cached_path.exists() and not rebuild and cache_matches(output_dir, lower_q=lower_q, upper_q=upper_q):
        return pd.read_parquet(cached_path)

    df = pd.read_parquet(input_path, columns=REQUIRED_COLUMNS)
    for col in ["precision", "recall", "f1", "gt_recall", "f1_residuals"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["sample", "tool", "match_type", "filter_type_comb", "f1", "gt_recall"]).copy()
    df["background"] = df["sample"].map(empirical_background)
    df["criterion_family"] = [
        criterion_family(a, b) for a, b in zip(df["filter_type_1"], df["filter_type_2"])
    ]
    df = add_group_residuals_once(df, lower_q=lower_q, upper_q=upper_q)
    df.to_parquet(cached_path, index=False)
    row_matrix_metadata_path(output_dir).write_text(
        json.dumps(
            {
                "input": str(input_path),
                "residual_method": "per_sample_tool_match_linear_f1_on_gt_recall",
                "f1_lower_quantile": float(lower_q),
                "f1_upper_quantile": float(upper_q),
                "n_rows": int(len(df)),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return df


def add_tiers(summary: pd.DataFrame) -> pd.DataFrame:
    out = summary.copy()
    min_s = float(out["stringency"].min())
    max_s = float(out["stringency"].max())
    third_1 = min_s + (max_s - min_s) / 3.0
    third_2 = min_s + 2.0 * (max_s - min_s) / 3.0
    conditions = [
        out["stringency"] <= third_1,
        (out["stringency"] > third_1) & (out["stringency"] <= third_2),
        out["stringency"] > third_2,
    ]
    out["tier"] = np.select(conditions, TIER_ORDER, default="unknown")
    return out


def summarize_criteria_with_residuals(df: pd.DataFrame) -> pd.DataFrame:
    residuals = df.dropna(subset=["f1_residual_cv"]).copy()
    if residuals.empty:
        return pd.DataFrame(
            columns=[
                "filter_type_comb",
                "filter_type_1",
                "filter_type_2",
                "criterion_family",
                "stringency",
                "mean_residual",
                "mean_f1",
                "n_units",
                "tier",
            ]
        )
    summary = (
        residuals.groupby(
            ["filter_type_comb", "filter_type_1", "filter_type_2", "criterion_family"],
            dropna=False,
            as_index=False,
        )
        .agg(
            stringency=("gt_recall", "mean"),
            mean_residual=("f1_residual_cv", "mean"),
            mean_f1=("f1", "mean"),
            n_units=("f1", "size"),
        )
    )
    return add_tiers(summary)


def select_criteria(
    train: pd.DataFrame,
    *,
    lower_q: float,
    upper_q: float,
) -> pd.DataFrame:
    del lower_q, upper_q
    summary = summarize_criteria_with_residuals(train)
    if summary.empty:
        return summary.copy()
    selected = []
    for tier in TIER_ORDER:
        section = summary[summary["tier"] == tier].dropna(subset=["mean_residual"])
        if section.empty:
            continue
        selected.append(section.sort_values(["mean_residual", "mean_f1"], ascending=False).iloc[0])
    return pd.DataFrame(selected).reset_index(drop=True)


def rank_tools(metric_df: pd.DataFrame) -> pd.DataFrame:
    ranked = metric_df.sort_values(["f1", "tool"], ascending=[False, True]).copy()
    ranked["rank"] = np.arange(1, len(ranked) + 1, dtype=int)
    return ranked


def tool_metrics_for_criterion(test: pd.DataFrame, criterion: str) -> pd.DataFrame:
    sub = test[test["filter_type_comb"].astype(str) == str(criterion)].copy()
    if sub.empty:
        return sub
    metrics = (
        sub.groupby("tool", as_index=False)
        .agg(
            precision=("precision", "mean"),
            recall=("recall", "mean"),
            f1=("f1", "mean"),
            n_units=("f1", "size"),
        )
    )
    return rank_tools(metrics)


def reference_criteria_from_full_data(df: pd.DataFrame) -> pd.DataFrame:
    selected = _pick_apa_filter_set(df)
    rows = []
    by_filter = (
        df.groupby(
            ["filter_type_comb", "filter_type_1", "filter_type_2", "criterion_family"],
            dropna=False,
            as_index=False,
        )
        .agg(
            stringency=("gt_recall", "mean"),
            mean_residual=("f1_residuals", "mean"),
            mean_f1=("f1", "mean"),
            n_units=("f1", "size"),
        )
    )
    by_filter = add_tiers(by_filter)
    for tier, criterion in zip(TIER_ORDER, selected):
        match = by_filter[by_filter["filter_type_comb"].astype(str) == str(criterion)]
        if match.empty:
            rows.append({"tier": tier, "filter_type_comb": criterion})
            continue
        current = match.iloc[0].to_dict()
        current["tier"] = tier
        rows.append(current)
    return pd.DataFrame(rows)


def full_data_rankings(df: pd.DataFrame) -> pd.DataFrame:
    selected = reference_criteria_from_full_data(df)
    rows = []
    for _, row in selected.iterrows():
        ranks = tool_metrics_for_criterion(df, str(row["filter_type_comb"]))
        ranks["tier"] = row["tier"]
        ranks["reference_criterion"] = row["filter_type_comb"]
        rows.append(ranks)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True, sort=False)


def concordance(heldout: pd.DataFrame, full_rank: pd.DataFrame, tier: str) -> tuple[float, float]:
    ref = full_rank[full_rank["tier"] == tier][["tool", "rank"]].rename(columns={"rank": "rank_full"})
    cur = heldout[["tool", "rank"]].rename(columns={"rank": "rank_heldout"})
    merged = ref.merge(cur, on="tool", how="inner")
    if len(merged) < 2:
        return np.nan, np.nan
    rho = spearmanr(merged["rank_full"], merged["rank_heldout"]).statistic
    tau = kendalltau(merged["rank_full"], merged["rank_heldout"]).statistic
    return float(rho), float(tau)


def criterion_concordance(
    train_summary: pd.DataFrame,
    heldout_summary: pd.DataFrame,
    tier: str,
) -> dict[str, float | int]:
    train_tier = train_summary.loc[
        train_summary["tier"].astype(str) == str(tier),
        ["filter_type_comb", "mean_residual"],
    ].rename(columns={"mean_residual": "mean_residual_train"})
    heldout_values = heldout_summary[["filter_type_comb", "mean_residual"]].rename(
        columns={"mean_residual": "mean_residual_heldout"}
    )
    merged = train_tier.merge(heldout_values, on="filter_type_comb", how="inner").dropna(
        subset=["mean_residual_train", "mean_residual_heldout"]
    )
    n_candidates = int(len(merged))
    if n_candidates < 2:
        return {"n_candidates": n_candidates, "spearman_rho": np.nan, "kendall_tau": np.nan}
    rho = spearmanr(merged["mean_residual_train"], merged["mean_residual_heldout"]).statistic
    tau = kendalltau(merged["mean_residual_train"], merged["mean_residual_heldout"]).statistic
    return {
        "n_candidates": n_candidates,
        "spearman_rho": float(rho),
        "kendall_tau": float(tau),
    }


def oracle_gap(
    test: pd.DataFrame,
    heldout_summary: pd.DataFrame,
    selected_criterion: str,
    candidate_criteria: list[str],
) -> dict[str, float]:
    selected_key = str(selected_criterion)
    candidates = [str(item) for item in candidate_criteria]
    selected_f1 = pd.to_numeric(
        test.loc[test["filter_type_comb"].astype(str) == selected_key, "f1"],
        errors="coerce",
    ).mean()
    oracle_f1 = (
        test[test["filter_type_comb"].astype(str).isin(candidates)]
        .groupby("filter_type_comb")["f1"]
        .mean()
        .max()
    )

    residual_lookup = heldout_summary.set_index(heldout_summary["filter_type_comb"].astype(str))["mean_residual"]
    selected_residual = residual_lookup.get(selected_key, np.nan)
    oracle_residual = residual_lookup.reindex(candidates).max()

    gap_f1 = np.nan if pd.isna(selected_f1) or pd.isna(oracle_f1) else float(oracle_f1 - selected_f1)
    gap_residual = (
        np.nan
        if pd.isna(selected_residual) or pd.isna(oracle_residual)
        else float(oracle_residual - selected_residual)
    )
    return {
        "gap_f1": gap_f1,
        "gap_residual": gap_residual,
        "selected_heldout_f1": float(selected_f1) if not pd.isna(selected_f1) else np.nan,
        "oracle_heldout_f1": float(oracle_f1) if not pd.isna(oracle_f1) else np.nan,
        "selected_heldout_residual": float(selected_residual) if not pd.isna(selected_residual) else np.nan,
        "oracle_heldout_residual": float(oracle_residual) if not pd.isna(oracle_residual) else np.nan,
    }


def evaluate_split(
    df: pd.DataFrame,
    train_mask: pd.Series,
    *,
    split_id: str,
    split_type: str,
    lower_q: float,
    upper_q: float,
    full_rank: pd.DataFrame,
) -> tuple[list[dict], list[dict], list[dict], list[dict], list[dict]]:
    train = df.loc[train_mask].copy()
    test = df.loc[~train_mask].copy()
    train_summary = summarize_criteria_with_residuals(train)
    selected = select_criteria(train, lower_q=lower_q, upper_q=upper_q)
    # The grouped split never cuts a sample/tool/match regression group, so the
    # cached per-group residuals are identical to residuals refit within held-out rows.
    heldout_summary = summarize_criteria_with_residuals(test)
    selected_rows: list[dict] = []
    metric_rows: list[dict] = []
    regret_rows: list[dict] = []
    concord_rows: list[dict] = []
    criterion_concord_rows: list[dict] = []

    for _, sel in selected.iterrows():
        tier = str(sel["tier"])
        criterion = str(sel["filter_type_comb"])
        selected_rows.append(
            {
                "split_id": split_id,
                "split_type": split_type,
                "tier": tier,
                "filter_type_comb": criterion,
                "filter_type_1": sel["filter_type_1"],
                "filter_type_2": sel["filter_type_2"],
                "criterion_family": sel["criterion_family"],
                "train_stringency": float(sel["stringency"]),
                "train_mean_residual": float(sel["mean_residual"]),
            }
        )
        ranks = tool_metrics_for_criterion(test, criterion)
        if ranks.empty:
            continue
        ranks["split_id"] = split_id
        ranks["split_type"] = split_type
        ranks["tier"] = tier
        ranks["filter_type_comb"] = criterion
        metric_rows.extend(ranks.to_dict("records"))
        rho, tau = concordance(ranks, full_rank, tier)
        concord_rows.append(
            {
                "split_id": split_id,
                "split_type": split_type,
                "tier": tier,
                "spearman_rho": rho,
                "kendall_tau": tau,
            }
        )
        criterion_stats = criterion_concordance(train_summary, heldout_summary, tier)
        criterion_concord_rows.append(
            {
                "split_id": split_id,
                "split_type": split_type,
                "tier": tier,
                **criterion_stats,
            }
        )
        tier_candidates = train_summary.loc[
            train_summary["tier"].astype(str) == tier, "filter_type_comb"
        ].astype(str).tolist()
        gap = oracle_gap(test, heldout_summary, criterion, tier_candidates)
        regret_rows.append(
            {
                "split_id": split_id,
                "split_type": split_type,
                "tier": tier,
                "filter_type_comb": criterion,
                **gap,
            }
        )
    return selected_rows, metric_rows, regret_rows, concord_rows, criterion_concord_rows


def make_random_splits(df: pd.DataFrame, unit_col: str, n: int, train_frac: float, seed: int):
    rng = np.random.default_rng(seed)
    units = np.array(sorted(df[unit_col].dropna().astype(str).unique()))
    n_train = max(1, min(len(units) - 1, int(round(len(units) * train_frac))))
    unit_series = df[unit_col].astype(str)
    for i in range(n):
        train_units = set(rng.choice(units, size=n_train, replace=False))
        yield f"random_{i + 1:03d}", unit_series.isin(train_units)


def make_protocol_splits(df: pd.DataFrame):
    protocols = sorted(df["protocol"].dropna().astype(str).unique())
    for protocol in protocols:
        yield f"leave_protocol_out__{protocol}", df["protocol"].astype(str) != protocol


def main() -> None:
    args = parse_args()
    apply_reference_style()
    matrix = load_or_build_matrix(
        args.input,
        args.output_dir,
        args.rebuild_row_matrix,
        lower_q=args.f1_lower_quantile,
        upper_q=args.f1_upper_quantile,
    )
    unit_col = "background" if args.split_unit == "background" else "sample"
    if unit_col not in matrix.columns:
        raise KeyError(f"Missing split unit column: {unit_col}")

    reference_criteria = reference_criteria_from_full_data(matrix)
    full_rank = full_data_rankings(matrix)
    selected_rows: list[dict] = []
    metric_rows: list[dict] = []
    regret_rows: list[dict] = []
    concord_rows: list[dict] = []
    criterion_concord_rows: list[dict] = []

    for split_id, train_mask in make_random_splits(
        matrix,
        unit_col,
        args.n_random_splits,
        args.train_frac,
        args.random_seed,
    ):
        parts = evaluate_split(
            matrix,
            train_mask,
            split_id=split_id,
            split_type="random_half_split",
            lower_q=args.f1_lower_quantile,
            upper_q=args.f1_upper_quantile,
            full_rank=full_rank,
        )
        for store, part in zip((selected_rows, metric_rows, regret_rows, concord_rows, criterion_concord_rows), parts):
            store.extend(part)

    if not args.skip_leave_protocol_out:
        for split_id, train_mask in make_protocol_splits(matrix):
            parts = evaluate_split(
                matrix,
                train_mask,
                split_id=split_id,
                split_type="leave_protocol_out",
                lower_q=args.f1_lower_quantile,
                upper_q=args.f1_upper_quantile,
                full_rank=full_rank,
            )
            for store, part in zip((selected_rows, metric_rows, regret_rows, concord_rows, criterion_concord_rows), parts):
                store.extend(part)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(selected_rows).to_csv(args.output_dir / "cv_selected_criteria.tsv", sep="\t", index=False)
    tool_metrics = pd.DataFrame(metric_rows)
    tool_metrics.to_csv(args.output_dir / "cv_tool_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(regret_rows).to_csv(args.output_dir / "cv_regret.tsv", sep="\t", index=False)
    pd.DataFrame(concord_rows).to_csv(args.output_dir / "cv_rank_concordance.tsv", sep="\t", index=False)
    pd.DataFrame(criterion_concord_rows).to_csv(
        args.output_dir / "cv_criterion_concordance.tsv", sep="\t", index=False
    )

    if not tool_metrics.empty:
        top = (
            tool_metrics.assign(is_top1=tool_metrics["rank"] == 1, is_top2=tool_metrics["rank"] <= 2)
            .groupby(["split_type", "tier", "tool"], as_index=False)
            .agg(top1_frequency=("is_top1", "mean"), top2_frequency=("is_top2", "mean"))
        )
        top.to_csv(args.output_dir / "cv_top_frequency.tsv", sep="\t", index=False)

    metadata = {
        "input": str(args.input),
        "row_matrix_table": str(row_matrix_path(args.output_dir)),
        "reference_criteria": reference_criteria.to_dict("records"),
        "output_dir": str(args.output_dir),
        "n_rows_matrix": int(len(matrix)),
        "n_random_splits": int(args.n_random_splits),
        "random_seed": int(args.random_seed),
        "train_frac": float(args.train_frac),
        "split_unit": args.split_unit,
        "split_unit_column": unit_col,
        "n_split_units": int(matrix[unit_col].dropna().astype(str).nunique()),
        "n_protocols": int(matrix["protocol"].dropna().astype(str).nunique()),
        "f1_lower_quantile": float(args.f1_lower_quantile),
        "f1_upper_quantile": float(args.f1_upper_quantile),
    }
    (args.output_dir / "cv_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(f"Wrote CV tables to: {args.output_dir}")


if __name__ == "__main__":
    main()
