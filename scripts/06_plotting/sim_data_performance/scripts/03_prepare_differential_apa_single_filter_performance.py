#!/usr/bin/env python3
"""Prepare single-filter DE-APA performance table from 05 performance outputs."""

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
DEFAULT_MAIN_PROJECT_ROOT = resolve_project_path("")
DE_APA_MAIN_SUFFIX = "_de_apa_performance.tsv"
DE_APA_DAPARS_PARALLEL_SUFFIX = "_de_apa_dapars_parallel_performance.tsv"
DEFAULT_FILTER_TYPE_1 = "scmapa_adPval_0.05"
DEFAULT_FILTER_TYPE_2 = "scmapa_or_0.1"
MAIN_ONLY_TOOLS = {"dapars2", "scmapa"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build single-filter DE-APA performance table from "
            "data/result/performance/*/*_de_apa*_performance.tsv"
        )
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
        "--filter-type-1",
        default=DEFAULT_FILTER_TYPE_1,
        help=(
            "Exact filter_type_1 to keep. Use empty string to disable this filter. "
            f"Default: {DEFAULT_FILTER_TYPE_1}"
        ),
    )
    parser.add_argument(
        "--filter-type-2",
        default=DEFAULT_FILTER_TYPE_2,
        help=(
            "Exact filter_type_2 to keep. Use empty string to disable this filter. "
            f"Default: {DEFAULT_FILTER_TYPE_2}"
        ),
    )
    parser.add_argument(
        "--tools",
        default="",
        help="Optional comma-separated raw tool list to include (e.g. 'dapars2,scmapa').",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="Optional debug cap on number of source files (0 means all).",
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


def _sample_from_filename(path: Path) -> str:
    name = path.name
    if name.endswith(DE_APA_DAPARS_PARALLEL_SUFFIX):
        return name.removesuffix(DE_APA_DAPARS_PARALLEL_SUFFIX)
    if name.endswith(DE_APA_MAIN_SUFFIX):
        return name.removesuffix(DE_APA_MAIN_SUFFIX)
    raise ValueError(f"Unexpected performance filename: {path}")


def _source_kind(path: Path) -> str:
    if path.name.endswith(DE_APA_DAPARS_PARALLEL_SUFFIX):
        return "dapars_parallel_sidecar"
    if path.name.endswith(DE_APA_MAIN_SUFFIX):
        return "main"
    return "unknown"


def main() -> None:
    args = parse_args()
    data_project_root, performance_dir = resolve_data_paths(args.data_root)
    out_dir = topic_intermediate_dir(TOPIC)

    selected_tools = {
        name.strip() for name in str(args.tools).split(",") if name and name.strip()
    }

    main_files = sorted(performance_dir.glob(f"*/*{DE_APA_MAIN_SUFFIX}"))
    sidecar_files = sorted(performance_dir.glob(f"*/*{DE_APA_DAPARS_PARALLEL_SUFFIX}"))
    all_files = sorted(main_files + sidecar_files)
    if selected_tools:
        all_files = [p for p in all_files if p.parent.name in selected_tools]

    if not all_files:
        raise FileNotFoundError(
            "No DE-APA performance files found under "
            f"{performance_dir} for selected tools: {sorted(selected_tools)}"
        )

    candidates: dict[tuple[str, str], dict[str, Path]] = {}
    for path in all_files:
        tool = path.parent.name
        sample = _sample_from_filename(path)
        key = (tool, sample)
        if key not in candidates:
            candidates[key] = {}
        candidates[key][_source_kind(path)] = path

    selected_file_map: dict[tuple[str, str], Path] = {}
    missing_required_sources: list[str] = []
    for (tool, sample), kinds in sorted(candidates.items()):
        if tool in MAIN_ONLY_TOOLS:
            chosen = kinds.get("main")
            if chosen is None:
                missing_required_sources.append(f"{tool}/{sample}: missing main *_de_apa_performance.tsv")
                continue
        else:
            chosen = kinds.get("dapars_parallel_sidecar")
            if chosen is None:
                missing_required_sources.append(
                    f"{tool}/{sample}: missing *_de_apa_dapars_parallel_performance.tsv"
                )
                continue
        selected_file_map[(tool, sample)] = chosen

    if missing_required_sources:
        preview = "\n".join(missing_required_sources[:10])
        raise ValueError(
            "Missing required DE-APA source files for selected policy "
            "(sidecar for non-dapars2/scmapa; main for dapars2/scmapa):\n"
            f"{preview}"
        )

    files = [selected_file_map[key] for key in sorted(selected_file_map)]
    if args.max_files and args.max_files > 0:
        files = files[: args.max_files]

    frames: list[pd.DataFrame] = []
    source_files: list[str] = []
    source_kind_by_tool_sample: dict[str, str] = {}

    for idx, table_path in enumerate(files, start=1):
        tool_name = table_path.parent.name
        sample_name = _sample_from_filename(table_path)
        df = read_table(table_path)
        df = _fill_identity_columns(df, tool=tool_name, sample=sample_name)
        df["tool_raw"] = tool_name
        df["source_kind"] = _source_kind(table_path)
        df["source_file"] = _safe_relative(table_path, data_project_root)
        frames.append(df)
        source_files.append(_safe_relative(table_path, data_project_root))
        source_kind_by_tool_sample[f"{tool_name}/{sample_name}"] = _source_kind(table_path)
        if idx % 2000 == 0:
            print(f"Loaded {idx}/{len(files)} DE-APA performance files ...")

    merged = pd.concat(frames, ignore_index=True, sort=False)

    required_cols = {"filter_type_1", "filter_type_2", "precision", "recall", "f1"}
    missing = sorted(required_cols.difference(merged.columns))
    if missing:
        raise KeyError(f"Missing required columns in DE-APA performance data: {missing}")

    merged["filter_type_1"] = merged["filter_type_1"].astype("string")
    merged["filter_type_2"] = merged["filter_type_2"].astype("string")

    out_frames: list[pd.DataFrame] = []
    selected_filter_by_tool: dict[str, str] = {}
    missing_primary_filter_tools: list[str] = []
    for tool_raw, group in merged.groupby("tool_raw", sort=False):
        wanted_f1, wanted_f2 = str(args.filter_type_1), str(args.filter_type_2)
        selected = group[
            (group["filter_type_1"] == wanted_f1) & (group["filter_type_2"] == wanted_f2)
        ].copy()
        if selected.empty:
            missing_primary_filter_tools.append(str(tool_raw))
            continue

        selected_filter_by_tool[str(tool_raw)] = f"{wanted_f1}|{wanted_f2}"
        out_frames.append(selected)

    out = pd.concat(out_frames, ignore_index=True, sort=False) if out_frames else pd.DataFrame()

    if out.empty:
        raise ValueError(
            "Single-filter table is empty after strict primary filter selection. "
            f"primary={args.filter_type_1!r}|{args.filter_type_2!r}."
        )
    if missing_primary_filter_tools:
        missing_unique = sorted(set(missing_primary_filter_tools))
        raise ValueError(
            "Primary filter pair is missing in selected source files for tools: "
            + ", ".join(missing_unique)
        )

    if "tool" in out.columns:
        out["tool"] = out["tool"].map(TOOL_MAP).fillna(out["tool"])
    if "protocol" in out.columns:
        out["protocol"] = out["protocol"].map(PROTOCOL_MAP).fillna(out["protocol"])

    for metric in ["precision", "recall", "f1", "gt_count", "pd_count"]:
        if metric in out.columns:
            out[metric] = pd.to_numeric(out[metric], errors="coerce")

    sort_cols = [name for name in ["tool", "sample", "match_type"] if name in out.columns]
    if sort_cols:
        out = out.sort_values(sort_cols).reset_index(drop=True)

    out_path = out_dir / "sim_data_differential_apa_single_filter_performance.parquet"
    meta_path = out_dir / "sim_data_differential_apa_single_filter_metadata.json"

    write_parquet(out, out_path)
    write_json(
        {
            "topic": TOPIC,
            "data_root_argument": str(args.data_root),
            "resolved_data_project_root": str(data_project_root),
            "resolved_performance_dir": str(performance_dir),
            "filter_type_1": str(args.filter_type_1),
            "filter_type_2": str(args.filter_type_2),
            "uniform_filter_applied": True,
            "source_selection_policy": (
                "non-dapars2/scmapa use *_de_apa_dapars_parallel_performance.tsv; "
                "dapars2/scmapa use *_de_apa_performance.tsv"
            ),
            "tools_argument": str(args.tools),
            "max_files": int(args.max_files),
            "source_file_count": int(len(files)),
            "rows_merged": int(len(merged)),
            "rows_output": int(len(out)),
            "tool_counts": {k: int(v) for k, v in out["tool"].value_counts().to_dict().items()},
            "selected_filter_by_tool_raw": selected_filter_by_tool,
            "source_kind_by_tool_sample": source_kind_by_tool_sample,
            "columns": list(out.columns),
            "source_files": source_files,
        },
        meta_path,
    )
    print(f"Wrote single-filter APA performance table: {out_path}")


if __name__ == "__main__":
    main()
