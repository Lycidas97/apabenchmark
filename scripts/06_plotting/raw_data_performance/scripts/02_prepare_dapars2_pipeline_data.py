#!/usr/bin/env python3
"""Prepare aggregated manifest for raw DaPars2/scMAPA pipeline outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Callable

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.io import write_json, write_parquet
from _shared.paths import resolve_project_path, topic_intermediate_dir, topic_root

TOPIC = "raw_data_performance"
DEFAULT_MAIN_PROJECT_ROOT = resolve_project_path("")
DEFAULT_CONFIG = topic_root(TOPIC) / "config" / "dapars2_pipeline.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate raw DaPars2/scMAPA pipeline outputs into manifest/summary tables."
    )
    parser.add_argument(
        "--data-root",
        default=str(DEFAULT_MAIN_PROJECT_ROOT),
        help="Project root (contains data/) or data directory.",
    )
    parser.add_argument(
        "--raw-run-id",
        default="default",
        help="Raw run id under data/result/raw (ignored when --run-root is set).",
    )
    parser.add_argument(
        "--run-root",
        default="",
        help="Explicit raw run root path (preferred for reproducibility).",
    )
    parser.add_argument(
        "--config",
        default=str(DEFAULT_CONFIG),
        help="JSON config for pipeline-domain and tool filtering.",
    )
    parser.add_argument(
        "--with-row-count",
        action="store_true",
        help="Also count rows for each source file (slower).",
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
    as_project_data = root / "data"
    if (as_project_data / "result" / "raw").exists():
        project_root = root
    elif (root / "result" / "raw").exists():
        project_root = root.parent
    else:
        raise FileNotFoundError(
            "Cannot resolve raw output base from --data-root. "
            f"Tried: {as_project_data / 'result/raw'} and {root / 'result/raw'}"
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
        return sum(1 for _ in fh)


def parse_converted_rel(rel: Path) -> tuple[str, str, str, str]:
    parts = rel.parts
    # {tool}/{sample}/{pair_label}/dapars2_format.txt
    tool = parts[0] if len(parts) >= 4 else ""
    sample = parts[1] if len(parts) >= 4 else ""
    pair_label = parts[2] if len(parts) >= 4 else ""
    comparison = ""
    return tool, sample, pair_label, comparison


def parse_differential_rel(rel: Path) -> tuple[str, str, str, str]:
    parts = rel.parts
    # {tool}/{sample}/{pair_label}/{comparison}/APA_results.txt
    tool = parts[0] if len(parts) >= 5 else ""
    sample = parts[1] if len(parts) >= 5 else ""
    pair_label = parts[2] if len(parts) >= 5 else ""
    comparison = parts[3] if len(parts) >= 5 else ""
    return tool, sample, pair_label, comparison


def parse_performance_rel(rel: Path) -> tuple[str, str, str, str]:
    parts = rel.parts
    # {tool}/{sample}/{pair_label}/dapars2_performance.tsv
    tool = parts[0] if len(parts) >= 4 else ""
    sample = parts[1] if len(parts) >= 4 else ""
    pair_label = parts[2] if len(parts) >= 4 else ""
    comparison = ""
    return tool, sample, pair_label, comparison


def collect_domain_records(
    *,
    domain_name: str,
    source_group: str,
    domain_dir: Path,
    pattern: str,
    parse_rel: Callable[[Path], tuple[str, str, str, str]],
    project_root: Path,
    with_row_count: bool,
    selected_tools: set[str] | None,
    include_ground_truth: bool,
) -> tuple[list[dict], list[str]]:
    if not domain_dir.exists():
        return [], []

    files = sorted(domain_dir.glob(pattern))
    records: list[dict] = []
    resolved_files: list[str] = []

    for path in files:
        rel_to_domain = path.relative_to(domain_dir)
        tool, sample, pair_label, comparison = parse_rel(rel_to_domain)

        if tool == "ground_truth" and not include_ground_truth:
            continue
        if selected_tools is not None and tool != "ground_truth" and tool not in selected_tools:
            continue

        header_cols = read_header_columns(path)
        row_count = count_data_rows(path) if with_row_count else pd.NA
        records.append(
            {
                "pipeline_domain": domain_name,
                "source_group": source_group,
                "tool": tool,
                "sample": sample,
                "pair_label": pair_label,
                "comparison": comparison,
                "n_columns": len(header_cols),
                "header_columns": "|".join(header_cols),
                "row_count": row_count,
                "file_size_bytes": path.stat().st_size,
                "source_file": _safe_relative(path, project_root),
            }
        )
        resolved_files.append(_safe_relative(path, project_root))

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
    config_path = Path(args.config).expanduser().resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Missing config file: {config_path}")

    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    pipeline_root_name = str(cfg.get("pipeline_root", "raw_dapars2_parallel"))
    domains_cfg = cfg.get("domains", {})
    tools_mode = str(cfg.get("tools_mode", "selected")).strip().lower()
    selected_tools_list = [str(v).strip() for v in cfg.get("selected_tools", []) if str(v).strip()]
    include_ground_truth = bool(cfg.get("include_ground_truth", True))

    selected_tools: set[str] | None = None
    if tools_mode == "selected":
        selected_tools = set(selected_tools_list)
    elif tools_mode == "all":
        selected_tools = None
    else:
        raise ValueError(f"Unsupported tools_mode: {tools_mode}. Use 'selected' or 'all'.")

    project_root, run_root = resolve_project_and_run_root(args.data_root, args.raw_run_id, args.run_root)
    pipeline_root = run_root / pipeline_root_name
    if not pipeline_root.exists():
        raise FileNotFoundError(f"Pipeline root not found: {pipeline_root}")

    domain_specs = []
    if bool(domains_cfg.get("converted_inputs", True)):
        domain_specs.append(
            {
                "domain_name": "raw_dapars2_parallel_converted_inputs",
                "source_group": "predicted",
                "domain_dir": pipeline_root / "converted_inputs",
                "pattern": "*/*/*/dapars2_format.txt",
                "parse_rel": parse_converted_rel,
            }
        )
    if bool(domains_cfg.get("differential_apa", True)):
        domain_specs.append(
            {
                "domain_name": "raw_dapars2_parallel_differential_apa",
                "source_group": "mixed",
                "domain_dir": pipeline_root / "differential_apa",
                "pattern": "*/*/*/*/APA_results.txt",
                "parse_rel": parse_differential_rel,
            }
        )
    if bool(domains_cfg.get("performance", True)):
        domain_specs.append(
            {
                "domain_name": "raw_dapars2_parallel_performance",
                "source_group": "predicted",
                "domain_dir": pipeline_root / "performance",
                "pattern": "*/*/*/dapars2_performance.tsv",
                "parse_rel": parse_performance_rel,
            }
        )

    all_records: list[dict] = []
    source_files_by_domain: dict[str, list[str]] = {}
    for spec in domain_specs:
        records, files = collect_domain_records(
            domain_name=spec["domain_name"],
            source_group=spec["source_group"],
            domain_dir=spec["domain_dir"],
            pattern=spec["pattern"],
            parse_rel=spec["parse_rel"],
            project_root=project_root,
            with_row_count=args.with_row_count,
            selected_tools=selected_tools,
            include_ground_truth=include_ground_truth,
        )
        all_records.extend(records)
        source_files_by_domain[spec["domain_name"]] = files

    if not all_records:
        raise FileNotFoundError(
            "No raw dapars2 pipeline outputs found for configured domains/tools under: "
            f"{pipeline_root}"
        )

    manifest = pd.DataFrame(all_records)
    if "row_count" in manifest.columns:
        manifest["row_count"] = manifest["row_count"].astype("Int64")
    manifest = manifest.sort_values(
        ["pipeline_domain", "tool", "sample", "pair_label", "comparison", "source_file"]
    ).reset_index(drop=True)
    summary = build_summary(manifest, with_row_count=args.with_row_count)

    out_dir = topic_intermediate_dir(TOPIC)
    manifest_path = out_dir / "raw_data_dapars2_pipeline_manifest.parquet"
    summary_path = out_dir / "raw_data_dapars2_pipeline_summary.parquet"
    metadata_path = out_dir / "raw_data_dapars2_pipeline_metadata.json"

    write_parquet(manifest, manifest_path)
    write_parquet(summary, summary_path)
    write_json(
        {
            "topic": TOPIC,
            "data_root_argument": str(args.data_root),
            "raw_run_id_argument": str(args.raw_run_id),
            "run_root_argument": str(args.run_root),
            "config": str(config_path),
            "resolved_project_root": str(project_root),
            "resolved_run_root": str(run_root),
            "pipeline_root": str(pipeline_root),
            "with_row_count": bool(args.with_row_count),
            "tools_mode": tools_mode,
            "selected_tools": sorted(selected_tools) if selected_tools is not None else "all",
            "include_ground_truth": include_ground_truth,
            "rows_manifest": int(len(manifest)),
            "rows_summary": int(len(summary)),
            "columns_manifest": list(manifest.columns),
            "columns_summary": list(summary.columns),
            "file_count_by_domain": {k: len(v) for k, v in source_files_by_domain.items()},
            "source_files_by_domain": source_files_by_domain,
        },
        metadata_path,
    )

    print(f"Wrote pipeline manifest: {manifest_path}")
    print(f"Wrote pipeline summary: {summary_path}")


if __name__ == "__main__":
    main()
