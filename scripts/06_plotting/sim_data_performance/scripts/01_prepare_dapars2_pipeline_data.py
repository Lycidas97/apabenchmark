#!/usr/bin/env python3
"""Prepare aggregated manifest for separate DaPars2/DE-APA pipeline outputs."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Callable

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.io import write_json, write_parquet
from _shared.paths import resolve_project_path, topic_intermediate_dir

TOPIC = "sim_data_performance"
DEFAULT_MAIN_PROJECT_ROOT = resolve_project_path("")


def _safe_relative(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate separate DaPars2-format and differential APA pipeline outputs."
    )
    parser.add_argument(
        "--data-root",
        default=str(DEFAULT_MAIN_PROJECT_ROOT),
        help=(
            "Path to benchmark project root (contains data/) or data directory. "
            "Default: repository checkout root."
        ),
    )
    parser.add_argument(
        "--with-row-count",
        action="store_true",
        help="Also count data rows per file (slower for large datasets).",
    )
    return parser.parse_args()


def resolve_data_paths(data_root_arg: str) -> tuple[Path, Path]:
    root = Path(data_root_arg).expanduser().resolve()
    as_project_data = root / "data"
    as_data_dir = root

    if (as_project_data / "result").exists():
        return root, as_project_data
    if (as_data_dir / "result").exists():
        return as_data_dir.parent, as_data_dir

    raise FileNotFoundError(
        "Cannot resolve data directory from --data-root. "
        f"Tried: {as_project_data / 'result'} and {as_data_dir / 'result'}"
    )


def read_header_columns(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        header = fh.readline().strip()
    if not header:
        return []
    return header.split("\t")


def count_data_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        header = fh.readline()
        if not header:
            return 0
        total = sum(1 for _ in fh)
    return total


def parse_dapars2_format_rel(rel: Path) -> tuple[str, str, str]:
    parts = rel.parts
    tool = parts[0] if len(parts) >= 2 else ""
    sample = rel.name.removesuffix("_dapars2_format.txt")
    return tool, sample, ""


def parse_dapars2_format_gt_rel(rel: Path) -> tuple[str, str, str]:
    sample = rel.name.removesuffix("_dapars2_format.txt")
    return "ground_truth", sample, ""


def parse_differential_apa_rel(rel: Path) -> tuple[str, str, str]:
    parts = rel.parts
    tool = parts[0] if len(parts) >= 4 else ""
    sample = parts[1] if len(parts) >= 4 else ""
    comparison = parts[2] if len(parts) >= 4 else ""
    return tool, sample, comparison


def parse_differential_apa_gt_rel(rel: Path) -> tuple[str, str, str]:
    parts = rel.parts
    sample = parts[0] if len(parts) >= 3 else ""
    comparison = parts[1] if len(parts) >= 3 else ""
    return "ground_truth", sample, comparison


def collect_domain_records(
    *,
    domain_name: str,
    source_group: str,
    domain_dir: Path,
    pattern: str,
    parse_rel: Callable[[Path], tuple[str, str, str]],
    data_project_root: Path,
    with_row_count: bool,
) -> tuple[list[dict], list[str]]:
    if not domain_dir.exists():
        return [], []

    files = sorted(domain_dir.glob(pattern))
    records: list[dict] = []
    resolved_files: list[str] = []

    for path in files:
        rel_to_domain = path.relative_to(domain_dir)
        tool, sample, comparison = parse_rel(rel_to_domain)
        header_cols = read_header_columns(path)
        row_count = count_data_rows(path) if with_row_count else pd.NA

        records.append(
            {
                "pipeline_domain": domain_name,
                "source_group": source_group,
                "tool": tool,
                "sample": sample,
                "comparison": comparison,
                "n_columns": len(header_cols),
                "header_columns": "|".join(header_cols),
                "row_count": row_count,
                "file_size_bytes": path.stat().st_size,
                "source_file": _safe_relative(path, data_project_root),
            }
        )
        resolved_files.append(_safe_relative(path, data_project_root))

    return records, resolved_files


def build_summary(manifest: pd.DataFrame, *, with_row_count: bool) -> pd.DataFrame:
    group_cols = ["pipeline_domain", "source_group", "tool", "comparison"]
    agg_spec: dict[str, tuple[str, str]] = {
        "file_count": ("source_file", "size"),
        "total_size_bytes": ("file_size_bytes", "sum"),
    }
    if with_row_count:
        agg_spec["total_row_count"] = ("row_count", "sum")
    summary = manifest.groupby(group_cols, dropna=False).agg(**agg_spec).reset_index()
    return summary.sort_values(group_cols).reset_index(drop=True)


def main() -> None:
    args = parse_args()
    data_project_root, data_dir = resolve_data_paths(args.data_root)
    out_dir = topic_intermediate_dir(TOPIC)

    domain_specs = [
        {
            "domain_name": "dapars2_format",
            "source_group": "predicted",
            "domain_dir": data_dir / "result" / "dapars2_format",
            "pattern": "*/*_dapars2_format.txt",
            "parse_rel": parse_dapars2_format_rel,
        },
        {
            "domain_name": "dapars2_format_ground_truth",
            "source_group": "ground_truth",
            "domain_dir": data_dir / "result" / "dapars2_format_ground_truth",
            "pattern": "*_dapars2_format.txt",
            "parse_rel": parse_dapars2_format_gt_rel,
        },
        {
            "domain_name": "differential_apa",
            "source_group": "predicted",
            "domain_dir": data_dir / "result" / "differential_apa",
            "pattern": "*/*/*/APA_results.txt",
            "parse_rel": parse_differential_apa_rel,
        },
        {
            "domain_name": "differential_apa_ground_truth",
            "source_group": "ground_truth",
            "domain_dir": data_dir / "result" / "differential_apa_ground_truth",
            "pattern": "*/*/APA_results.txt",
            "parse_rel": parse_differential_apa_gt_rel,
        },
    ]

    all_records: list[dict] = []
    source_files_by_domain: dict[str, list[str]] = {}
    for spec in domain_specs:
        records, files = collect_domain_records(
            domain_name=spec["domain_name"],
            source_group=spec["source_group"],
            domain_dir=spec["domain_dir"],
            pattern=spec["pattern"],
            parse_rel=spec["parse_rel"],
            data_project_root=data_project_root,
            with_row_count=args.with_row_count,
        )
        all_records.extend(records)
        source_files_by_domain[spec["domain_name"]] = files

    if not all_records:
        raise FileNotFoundError(
            "No separate pipeline outputs found under: "
            f"{data_dir / 'result' / 'dapars2_format'}, "
            f"{data_dir / 'result' / 'dapars2_format_ground_truth'}, "
            f"{data_dir / 'result' / 'differential_apa'}, "
            f"{data_dir / 'result' / 'differential_apa_ground_truth'}"
        )

    manifest = pd.DataFrame(all_records)
    if "row_count" in manifest.columns:
        manifest["row_count"] = manifest["row_count"].astype("Int64")
    manifest = manifest.sort_values(
        ["pipeline_domain", "tool", "sample", "comparison", "source_file"]
    ).reset_index(drop=True)
    summary = build_summary(manifest, with_row_count=args.with_row_count)

    manifest_path = out_dir / "sim_data_dapars2_pipeline_manifest.parquet"
    summary_path = out_dir / "sim_data_dapars2_pipeline_summary.parquet"
    write_parquet(manifest, manifest_path)
    write_parquet(summary, summary_path)
    write_json(
        {
            "topic": TOPIC,
            "data_root_argument": str(args.data_root),
            "resolved_data_project_root": str(data_project_root),
            "resolved_data_dir": str(data_dir),
            "with_row_count": bool(args.with_row_count),
            "rows_manifest": int(len(manifest)),
            "rows_summary": int(len(summary)),
            "columns_manifest": list(manifest.columns),
            "columns_summary": list(summary.columns),
            "file_count_by_domain": {key: len(value) for key, value in source_files_by_domain.items()},
            "source_files_by_domain": source_files_by_domain,
        },
        out_dir / "sim_data_dapars2_pipeline_metadata.json",
    )
    print(f"Wrote pipeline manifest: {manifest_path}")
    print(f"Wrote pipeline summary: {summary_path}")


if __name__ == "__main__":
    main()
