#!/usr/bin/env python3
"""Prepare intermediate tables for computation resource consumption plots."""

from __future__ import annotations

import os
from pathlib import Path
import re
import sys

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import TOOL_MAP
from _shared.io import read_table, write_json, write_tsv
from _shared.paths import resolve_project_path, topic_intermediate_dir

TOPIC = "computation_resource_consumption"
INPUT_ENV = "CP_RESOURCE_INPUT_DIR"
INPUT_DEFAULT = "data/result/cp_resource/pf_bam"
OUTPUT_RAW = "cp_resource_consumption_merged_raw.tsv"
OUTPUT_PREPARED = "cp_resource_consumption_prepared.tsv"

FILE_PATTERN = re.compile(
    r"^ps(?P<peak_size>\d+)_gn(?P<gene_number>\d+)_rl(?P<read_length>\d+)"
    r"_bn(?P<barcode_number>\d+)_pas(?P<pas_id>\d+)_rep(?P<replicate>\d+)\.txt$"
)


def resolve_source_dir() -> Path:
    source = os.environ.get(INPUT_ENV, INPUT_DEFAULT)
    source_path = Path(source).expanduser()
    if not source_path.is_absolute():
        source_path = resolve_project_path(str(source_path))
    return source_path.resolve()


def parse_file(tool_raw: str, file_path: Path) -> dict[str, object]:
    matched = FILE_PATTERN.match(file_path.name)
    if not matched:
        raise ValueError(f"Unexpected filename format: {file_path}")

    table = read_table(file_path)
    if table.empty:
        raise ValueError(f"Empty profile file: {file_path}")

    row = table.iloc[0].to_dict()
    row.update({key: int(value) for key, value in matched.groupdict().items()})
    row["tool_raw"] = tool_raw
    row["source_file"] = str(file_path)
    return row


def main() -> None:
    source_dir = resolve_source_dir()
    if not source_dir.exists():
        raise FileNotFoundError(
            f"Missing input directory: {source_dir}. "
            f"Set {INPUT_ENV} to override the source path."
        )

    files: list[tuple[str, Path]] = []
    for tool_dir in sorted(path for path in source_dir.iterdir() if path.is_dir()):
        for file_path in sorted(tool_dir.glob("*.txt")):
            files.append((tool_dir.name, file_path))

    if not files:
        raise FileNotFoundError(f"No profile text files found in {source_dir}")

    merged = pd.DataFrame.from_records(parse_file(tool_raw, path) for tool_raw, path in files)
    merged = merged.rename(columns={"s": "runtime", "max_rss": "max_mem"})
    merged["runtime_min"] = pd.to_numeric(merged["runtime"], errors="coerce") / 60.0
    merged["tool"] = merged["tool_raw"].map(TOOL_MAP).fillna(merged["tool_raw"])

    numeric_columns = [
        "peak_size",
        "gene_number",
        "read_length",
        "barcode_number",
        "pas_id",
        "replicate",
        "runtime",
        "runtime_min",
        "max_mem",
        "mean_load",
        "io_in",
        "io_out",
        "cpu_time",
    ]
    for column in numeric_columns:
        if column in merged.columns:
            merged[column] = pd.to_numeric(merged[column], errors="coerce")

    sort_columns = [
        "tool",
        "peak_size",
        "gene_number",
        "read_length",
        "barcode_number",
        "pas_id",
        "replicate",
    ]
    merged = merged.sort_values(sort_columns).reset_index(drop=True)

    prepared_columns = [
        "tool",
        "tool_raw",
        "peak_size",
        "gene_number",
        "read_length",
        "barcode_number",
        "pas_id",
        "replicate",
        "runtime",
        "runtime_min",
        "max_mem",
        "mean_load",
        "io_in",
        "io_out",
        "cpu_time",
    ]
    prepared = merged[prepared_columns].copy()

    out_dir = topic_intermediate_dir(TOPIC)
    raw_path = out_dir / OUTPUT_RAW
    prepared_path = out_dir / OUTPUT_PREPARED
    write_tsv(merged, raw_path)
    write_tsv(prepared, prepared_path)
    write_json(
        {
            "topic": TOPIC,
            "source_dir": str(source_dir),
            "input_env": INPUT_ENV,
            "input_files": int(len(files)),
            "rows": int(len(prepared)),
            "tools_raw": sorted(prepared["tool_raw"].dropna().unique().tolist()),
            "tools_mapped": sorted(prepared["tool"].dropna().unique().tolist()),
            "columns": prepared.columns.tolist(),
            "output_raw": str(raw_path),
            "output_prepared": str(prepared_path),
        },
        out_dir / "metadata.json",
    )
    print(f"Wrote merged raw table: {raw_path}")
    print(f"Wrote prepared table: {prepared_path}")


if __name__ == "__main__":
    main()
