#!/usr/bin/env python3
"""Prepare sample-level normalized peak profile curves for protocol overlay plotting."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.io import read_table, write_json, write_tsv
from _shared.paths import resolve_project_path, topic_intermediate_dir

TOPIC = "peak_model_params"
INPUT_NAME = "peak_model_params_prepared.tsv"
OUTPUT_NAME = "protocol_normalized_peak_profiles.tsv"
METADATA_NAME = "protocol_normalized_peak_profiles_metadata.json"
REQUIRED_COLUMNS = {"sample_full_id", "protocol"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare normalized per-sample peak profile curves for protocol overlay plotting"
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=topic_intermediate_dir(TOPIC) / INPUT_NAME,
        help="Prepared sample table",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=topic_intermediate_dir(TOPIC) / OUTPUT_NAME,
        help="Output long-format normalized profile TSV",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=topic_intermediate_dir(TOPIC) / METADATA_NAME,
        help="Output metadata JSON",
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=None,
        help="Project root containing data/ or the data directory itself.",
    )
    parser.add_argument(
        "--window-start",
        type=int,
        default=-500,
        help="Coverage window start (relative to PAS)",
    )
    parser.add_argument(
        "--window-end",
        type=int,
        default=200,
        help="Coverage window end (relative to PAS)",
    )
    parser.add_argument(
        "--min-reads-per-pas",
        type=int,
        default=50,
        help="Keep PAS with read count strictly greater than this threshold",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Optional cap for number of samples (smoke tests)",
    )
    return parser.parse_args()


def resolve_data_root(data_root: Path | None) -> Path | None:
    if data_root is None:
        return None
    root = data_root.expanduser().resolve()
    if root.name == "data":
        root = root.parent
    if not (root / "data").is_dir():
        raise FileNotFoundError(f"Missing data directory under project root: {root}")
    return root


def reads_to_coverage(
    distance_to_pas: np.ndarray,
    reference_length: np.ndarray,
    window_start: int,
    window_end: int,
) -> np.ndarray:
    """Convert read-level records to coverage using interval [d-L, d)."""
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

    return np.cumsum(diff[:-1]).astype(np.float64)


def resolve_peak_path(sample_full_id: str, shape_input_file: Any, data_root: Path | None = None) -> Path | None:
    """Resolve peak feather path from the table field plus robust fallbacks."""
    candidates: list[Path] = []

    def add_candidate(path: Path) -> None:
        if path not in candidates:
            candidates.append(path)

    if pd.notna(shape_input_file):
        raw = str(shape_input_file).strip()
        if raw:
            raw_path = Path(raw)
            if raw_path.is_absolute():
                add_candidate(raw_path)
            else:
                add_candidate((Path.cwd() / raw_path).resolve())
                add_candidate(resolve_project_path(raw).resolve())
                add_candidate((resolve_project_path("scripts/02_peak_analysis") / raw_path).resolve())

    sample_peak_name = f"{sample_full_id}_peak.feather"
    root = resolve_data_root(data_root)
    if root is not None:
        add_candidate((root / "data" / "int_data" / "peaks" / sample_peak_name).resolve())
    add_candidate(resolve_project_path(f"data/int_data/peaks/{sample_peak_name}").resolve())

    for path in candidates:
        if path.exists():
            return path
    return None


def build_sample_curve(
    peak_path: Path,
    window_start: int,
    window_end: int,
    min_reads_per_pas: int,
) -> tuple[np.ndarray | None, dict[str, int | str]]:
    """Build normalized raw coverage curve for one sample."""
    try:
        frame = pd.read_feather(peak_path, columns=["pas", "distance_to_pas", "reference_length"])
    except Exception as exc:  # pragma: no cover - defensive path
        return None, {"status": "read_failed", "detail": str(exc)}

    if frame.empty:
        return None, {"status": "empty_peak_file", "detail": "no rows"}

    numeric = frame[["pas", "distance_to_pas", "reference_length"]].copy()
    numeric[["distance_to_pas", "reference_length"]] = numeric[
        ["distance_to_pas", "reference_length"]
    ].apply(pd.to_numeric, errors="coerce")
    numeric = numeric.dropna(subset=["pas", "distance_to_pas", "reference_length"])
    if numeric.empty:
        return None, {"status": "no_numeric_rows", "detail": "all rows invalid"}

    pas_sizes = numeric.groupby("pas", sort=False).size()
    valid_pas = pas_sizes[pas_sizes > min_reads_per_pas].index
    filtered = numeric[numeric["pas"].isin(valid_pas)]
    if filtered.empty:
        return None, {
            "status": "no_pas_after_filter",
            "detail": f"min_reads_per_pas>{min_reads_per_pas}",
        }

    distance = filtered["distance_to_pas"].to_numpy(dtype=np.int64)
    reference = filtered["reference_length"].to_numpy(dtype=np.int64)
    coverage = reads_to_coverage(distance, reference, window_start, window_end)

    total = float(np.sum(coverage))
    if total <= 0:
        return None, {"status": "zero_coverage", "detail": "coverage sum <= 0"}

    return coverage / total, {
        "status": "ok",
        "n_reads_after_filter": int(len(filtered)),
        "n_pas_after_filter": int(len(valid_pas)),
    }


def main() -> None:
    args = parse_args()
    if args.window_end < args.window_start:
        raise ValueError("window_end must be greater than or equal to window_start")

    df = read_table(args.input)
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in input: {sorted(missing)}")

    if "shape_input_file" not in df.columns:
        df = df.copy()
        df["shape_input_file"] = np.nan

    samples = (
        df[["sample_full_id", "protocol", "shape_input_file"]]
        .dropna(subset=["sample_full_id", "protocol"])
        .drop_duplicates(subset=["sample_full_id"], keep="first")
        .sort_values(["protocol", "sample_full_id"], kind="stable")
        .reset_index(drop=True)
    )

    if args.max_samples is not None:
        samples = samples.head(args.max_samples).copy()

    positions = np.arange(args.window_start, args.window_end + 1, dtype=np.int64)

    output_parts: list[pd.DataFrame] = []
    skipped: list[dict[str, str]] = []
    processed = 0

    for row in samples.itertuples(index=False):
        peak_path = resolve_peak_path(row.sample_full_id, row.shape_input_file, args.data_root)
        if peak_path is None:
            skipped.append(
                {
                    "sample_full_id": str(row.sample_full_id),
                    "protocol": str(row.protocol),
                    "status": "missing_peak_file",
                    "detail": str(row.shape_input_file),
                }
            )
            continue

        curve, info = build_sample_curve(
            peak_path=peak_path,
            window_start=args.window_start,
            window_end=args.window_end,
            min_reads_per_pas=args.min_reads_per_pas,
        )
        if curve is None:
            skipped.append(
                {
                    "sample_full_id": str(row.sample_full_id),
                    "protocol": str(row.protocol),
                    "status": str(info.get("status", "unknown")),
                    "detail": str(info.get("detail", "")),
                }
            )
            continue

        output_parts.append(
            pd.DataFrame(
                {
                    "sample_full_id": row.sample_full_id,
                    "protocol": row.protocol,
                    "position": positions,
                    "coverage_norm": curve,
                }
            )
        )
        processed += 1

    if output_parts:
        output_df = pd.concat(output_parts, ignore_index=True)
    else:
        output_df = pd.DataFrame(
            columns=["sample_full_id", "protocol", "position", "coverage_norm"]
        )

    write_tsv(output_df, args.output)

    skipped_df = pd.DataFrame(skipped)
    if skipped_df.empty:
        status_counts: dict[str, int] = {}
    else:
        status_counts = skipped_df["status"].value_counts().to_dict()

    write_json(
        {
            "topic": TOPIC,
            "input": str(args.input),
            "output": str(args.output),
            "window_start": int(args.window_start),
            "window_end": int(args.window_end),
            "min_reads_per_pas": int(args.min_reads_per_pas),
            "max_samples": None if args.max_samples is None else int(args.max_samples),
            "n_samples_in_input": int(df["sample_full_id"].nunique()),
            "n_samples_attempted": int(len(samples)),
            "n_samples_processed": int(processed),
            "n_samples_skipped": int(len(skipped)),
            "skipped_status_counts": status_counts,
            "n_output_rows": int(len(output_df)),
            "n_protocols": int(output_df["protocol"].nunique()) if not output_df.empty else 0,
        },
        args.metadata,
    )

    print(f"Wrote normalized profile table: {args.output}")
    print(f"Wrote metadata: {args.metadata}")


if __name__ == "__main__":
    main()
