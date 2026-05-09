#!/usr/bin/env python3
"""Prepare intermediate tables for raw-data performance plotting."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import PROTOCOL_MAP, TOOL_MAP
from _shared.io import read_table, write_json, write_parquet
from _shared.paths import resolve_project_path, topic_intermediate_dir

TOPIC = "raw_data_performance"
DEFAULT_MAIN_PROJECT_ROOT = resolve_project_path("")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate raw benchmark summaries into plotting intermediate files."
    )
    parser.add_argument(
        "--data-root",
        default=str(DEFAULT_MAIN_PROJECT_ROOT),
        help=(
            "Path to benchmark project root (contains data/) or data directory "
            "(contains result/raw)."
        ),
    )
    parser.add_argument(
        "--raw-run-id",
        default="default",
        help="Raw run id under data/result/raw (default: default).",
    )
    parser.add_argument(
        "--run-root",
        default="",
        help="Explicit raw run root path. If set, overrides --data-root/--raw-run-id.",
    )
    return parser.parse_args()


def resolve_project_and_run_root(data_root_arg: str, raw_run_id: str, run_root_arg: str) -> tuple[Path, Path]:
    if run_root_arg.strip():
        run_root = Path(run_root_arg).expanduser().resolve()
        if not run_root.exists():
            raise FileNotFoundError(f"run root does not exist: {run_root}")
        for parent in [run_root.parent.parent.parent, run_root.parent.parent, run_root.parent]:
            if (parent / "data").exists():
                return parent, run_root
        return run_root.parent, run_root

    root = Path(data_root_arg).expanduser().resolve()
    project_data = root / "data"
    if (project_data / "result" / "raw").exists():
        project_root = root
    elif (root / "result" / "raw").exists():
        project_root = root.parent
    else:
        raise FileNotFoundError(
            "Cannot resolve raw output base from --data-root. "
            f"Tried: {project_data / 'result/raw'} and {root / 'result/raw'}"
        )

    run_root = project_root / "data" / "result" / "raw" / raw_run_id
    if not run_root.exists():
        raise FileNotFoundError(f"Raw run root not found: {run_root}")
    return project_root, run_root


def _safe_relative(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _add_identity_columns(df: pd.DataFrame, metric_domain: str, source_file: str, raw_run_id: str) -> pd.DataFrame:
    df = df.copy()
    df["metric_domain"] = metric_domain
    df["source_file"] = source_file
    df["raw_run_id"] = raw_run_id

    if "sample" in df.columns:
        sample_series = df["sample"].astype(str)
        parts = sample_series.str.split("__", n=1, expand=True)
        if parts.shape[1] == 2:
            df["sample_id"] = parts.iloc[:, 0]
            df["pair_label"] = parts.iloc[:, 1]

    if "tool" in df.columns:
        df["tool"] = df["tool"].map(TOOL_MAP).fillna(df["tool"])
    if "protocol" in df.columns:
        df["protocol"] = df["protocol"].map(PROTOCOL_MAP).fillna(df["protocol"])

    return df


def _read_optional_table(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    return read_table(path)


def _read_metric_table(path: Path, *, metric_domain: str, project_root: Path, raw_run_id: str) -> pd.DataFrame | None:
    df = _read_optional_table(path)
    if df is None:
        return None
    return _add_identity_columns(
        df,
        metric_domain=metric_domain,
        source_file=_safe_relative(path, project_root),
        raw_run_id=raw_run_id,
    )


def main() -> None:
    args = parse_args()
    project_root, run_root = resolve_project_and_run_root(
        args.data_root,
        args.raw_run_id,
        args.run_root,
    )
    raw_run_id = run_root.name

    metrics_sources: list[tuple[str, Path]] = [
        (
            "raw_match_performance",
            run_root / "raw_performance" / "summary" / "raw_performance_match_metrics.tsv",
        ),
        (
            "raw_pas_quantify_performance",
            run_root / "raw_performance" / "summary" / "raw_performance_pas_quantify_metrics.tsv",
        ),
        (
            "raw_de_apa_performance",
            run_root / "raw_performance" / "summary" / "raw_performance_de_apa_metrics.tsv",
        ),
    ]

    parallel_dapars_path = (
        run_root / "raw_dapars2_parallel" / "summary" / "raw_dapars2_parallel_metrics.tsv"
    )
    format_dapars_path = run_root / "raw_dapars2_format" / "summary" / "raw_dapars2_metrics.tsv"
    if parallel_dapars_path.exists():
        metrics_sources.append(("raw_dapars2_performance", parallel_dapars_path))
    else:
        metrics_sources.append(("raw_dapars2_performance", format_dapars_path))

    metric_frames: list[pd.DataFrame] = []
    discovered_metric_files: dict[str, str] = {}
    missing_metric_files: dict[str, str] = {}
    for domain, path in metrics_sources:
        df = _read_metric_table(
            path,
            metric_domain=domain,
            project_root=project_root,
            raw_run_id=raw_run_id,
        )
        if df is None:
            missing_metric_files[domain] = _safe_relative(path, project_root)
            continue
        discovered_metric_files[domain] = _safe_relative(path, project_root)
        metric_frames.append(df)

    if not metric_frames:
        expected = [str(path) for _, path in metrics_sources]
        raise FileNotFoundError(f"No metric summary tables found. Expected one of: {expected}")

    merged_metrics = pd.concat(metric_frames, ignore_index=True, sort=False)
    sort_columns = [
        col
        for col in ["metric_domain", "tool", "sample_id", "pair_label", "sample"]
        if col in merged_metrics.columns
    ]
    if sort_columns:
        merged_metrics = merged_metrics.sort_values(sort_columns).reset_index(drop=True)

    pair_manifest_path = run_root / "raw_pair_manifest" / "sample_pair_manifest.tsv"
    pair_manifest_df = _read_optional_table(pair_manifest_path)
    if pair_manifest_df is None:
        pair_manifest_df = pd.DataFrame()

    qc_path = run_root / "raw_qc" / "pas_filter_overview.tsv"
    qc_df = _read_optional_table(qc_path)
    if qc_df is None:
        qc_df = pd.DataFrame()

    failure_paths = [
        run_root / "raw_performance" / "summary" / "raw_performance_missing_failed.tsv",
        run_root / "raw_dapars2_parallel" / "summary" / "raw_dapars2_parallel_missing_failed.tsv",
        run_root / "raw_dapars2_format" / "summary" / "raw_dapars2_missing_failed.tsv",
    ]
    failure_frames: list[pd.DataFrame] = []
    discovered_failure_files: list[str] = []
    for path in failure_paths:
        df = _read_optional_table(path)
        if df is None:
            continue
        df = df.copy()
        df["source_file"] = _safe_relative(path, project_root)
        df["raw_run_id"] = raw_run_id
        failure_frames.append(df)
        discovered_failure_files.append(_safe_relative(path, project_root))

    failures_df = pd.concat(failure_frames, ignore_index=True, sort=False) if failure_frames else pd.DataFrame()

    out_dir = topic_intermediate_dir(TOPIC)
    metrics_out = out_dir / "raw_data_performance_prepared.parquet"
    manifest_out = out_dir / "raw_data_pair_manifest.parquet"
    qc_out = out_dir / "raw_data_qc_overview.parquet"
    failures_out = out_dir / "raw_data_failures.parquet"

    write_parquet(merged_metrics, metrics_out)
    write_parquet(pair_manifest_df, manifest_out)
    write_parquet(qc_df, qc_out)
    write_parquet(failures_df, failures_out)

    metadata = {
        "topic": TOPIC,
        "data_root_argument": str(args.data_root),
        "raw_run_id_argument": str(args.raw_run_id),
        "run_root_argument": str(args.run_root),
        "resolved_project_root": str(project_root),
        "resolved_run_root": str(run_root),
        "plot_project_root": str(resolve_project_path("")),
        "rows": {
            "metrics": int(len(merged_metrics)),
            "pair_manifest": int(len(pair_manifest_df)),
            "qc": int(len(qc_df)),
            "failures": int(len(failures_df)),
        },
        "columns": {
            "metrics": list(merged_metrics.columns),
            "pair_manifest": list(pair_manifest_df.columns),
            "qc": list(qc_df.columns),
            "failures": list(failures_df.columns),
        },
        "discovered_metric_files": discovered_metric_files,
        "missing_metric_files": missing_metric_files,
        "discovered_failure_files": discovered_failure_files,
    }
    write_json(metadata, out_dir / "metadata.json")

    print(f"Wrote metrics intermediate: {metrics_out}")
    print(f"Wrote pair manifest intermediate: {manifest_out}")
    print(f"Wrote QC intermediate: {qc_out}")
    print(f"Wrote failures intermediate: {failures_out}")


if __name__ == "__main__":
    main()
