#!/usr/bin/env python3
"""Prepare aggregated intermediate table for sim data performance plots."""

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

TOPIC = "sim_data_performance"
OUTPUT_NAME = "sim_data_performance_prepared.parquet"
DEFAULT_MAIN_PROJECT_ROOT = resolve_project_path("")
FILE_SUFFIX_BY_DOMAIN = {
    "match_performance": "_match_performance.tsv",
    "pas_quantify_performance": "_pas_quantify_performance.tsv",
    "de_apa_performance": "_de_apa_performance.tsv",
}


def _safe_relative(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _fill_identity_columns(df: pd.DataFrame, *, tool: str, sample: str) -> pd.DataFrame:
    if "tool" in df.columns:
        df["tool"] = df["tool"].fillna(tool).replace("", tool)
    else:
        df["tool"] = tool

    if "sample" in df.columns:
        df["sample"] = df["sample"].fillna(sample).replace("", sample)
    else:
        df["sample"] = sample

    default_protocol = sample.split("_")[0] if sample else ""
    if "protocol" in df.columns:
        df["protocol"] = df["protocol"].fillna(default_protocol).replace("", default_protocol)
    else:
        df["protocol"] = default_protocol
    return df


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate sim data performance tables into one intermediate file."
    )
    parser.add_argument(
        "--data-root",
        default=str(DEFAULT_MAIN_PROJECT_ROOT),
        help=(
            "Path to benchmark project root (contains data/) or data directory "
            "(contains result/performance). Default: repository checkout root."
        ),
    )
    parser.add_argument(
        "--include-stage-files",
        action="store_true",
        help="Include *_stage10_* and *_stage20_* performance files (default: excluded).",
    )
    parser.add_argument(
        "--max-files-per-domain",
        type=int,
        default=0,
        help="Optional debug cap of source files per metric domain (0 means all).",
    )
    return parser.parse_args()


def resolve_data_paths(data_root_arg: str) -> tuple[Path, Path]:
    root = Path(data_root_arg).expanduser().resolve()
    as_project_data = root / "data"
    as_data_dir = root

    project_perf = as_project_data / "result" / "performance"
    data_perf = as_data_dir / "result" / "performance"

    if project_perf.exists():
        return root, project_perf
    if data_perf.exists():
        return as_data_dir.parent, data_perf

    raise FileNotFoundError(
        "Cannot resolve performance directory from --data-root. "
        f"Tried: {project_perf} and {data_perf}"
    )


def main() -> None:
    args = parse_args()
    data_project_root, performance_dir = resolve_data_paths(args.data_root)
    plot_project_root = resolve_project_path("")
    frames: list[pd.DataFrame] = []
    discovered_files: dict[str, list[str]] = {}

    for domain, suffix in FILE_SUFFIX_BY_DOMAIN.items():
        files = sorted(performance_dir.glob(f"*/*{suffix}"))
        if not args.include_stage_files:
            files = [
                path
                for path in files
                if "_stage10_" not in path.name and "_stage20_" not in path.name
            ]
        if args.max_files_per_domain and args.max_files_per_domain > 0:
            files = files[: args.max_files_per_domain]
        discovered_files[domain] = [_safe_relative(path, data_project_root) for path in files]
        for table_path in files:
            tool_name = table_path.parent.name
            sample_name = table_path.name.removesuffix(suffix)
            df = read_table(table_path)
            df = _fill_identity_columns(df, tool=tool_name, sample=sample_name)
            df["metric_domain"] = domain
            df["source_file"] = _safe_relative(table_path, data_project_root)
            frames.append(df)

    if not frames:
        expected_patterns = ", ".join(
            f"data/result/performance/*/*{suffix}" for suffix in FILE_SUFFIX_BY_DOMAIN.values()
        )
        raise FileNotFoundError(f"No performance tables found. Expected patterns: {expected_patterns}")

    merged = pd.concat(frames, ignore_index=True, sort=False)
    if "tool" in merged.columns:
        merged["tool"] = merged["tool"].map(TOOL_MAP).fillna(merged["tool"])
    if "protocol" in merged.columns:
        merged["protocol"] = merged["protocol"].map(PROTOCOL_MAP).fillna(merged["protocol"])

    sort_columns = [name for name in ("metric_domain", "tool", "sample") if name in merged.columns]
    if sort_columns:
        merged = merged.sort_values(sort_columns).reset_index(drop=True)

    out_dir = topic_intermediate_dir(TOPIC)
    out_path = out_dir / OUTPUT_NAME
    write_parquet(merged, out_path)
    write_json(
        {
            "topic": TOPIC,
            "data_root_argument": str(args.data_root),
            "include_stage_files": bool(args.include_stage_files),
            "max_files_per_domain": int(args.max_files_per_domain),
            "resolved_data_project_root": str(data_project_root),
            "resolved_performance_dir": str(performance_dir),
            "plot_project_root": str(plot_project_root),
            "rows": int(len(merged)),
            "columns": list(merged.columns),
            "file_count_by_domain": {key: len(value) for key, value in discovered_files.items()},
            "source_files_by_domain": discovered_files,
        },
        out_dir / "metadata.json",
    )
    print(f"Wrote intermediate table: {out_path}")


if __name__ == "__main__":
    main()
