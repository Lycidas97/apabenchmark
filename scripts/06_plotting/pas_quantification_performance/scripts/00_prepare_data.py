#!/usr/bin/env python3
"""Prepare intermediate table for PAS quantification performance plots."""

from __future__ import annotations

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[2]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import PROTOCOL_MAP, TOOL_MAP
from _shared.io import read_table, write_json, write_tsv
from _shared.paths import resolve_project_path, topic_intermediate_dir

TOPIC = "pas_quantification_performance"
SOURCE_RELATIVE = "data/result/merged_tool_performance/pas_quantification_performance.tsv"
OUTPUT_NAME = "pas_quantification_performance_prepared.tsv"


def main() -> None:
    source_path = resolve_project_path(SOURCE_RELATIVE)
    if not source_path.exists():
        raise FileNotFoundError(f"Missing input table: {source_path}")

    df = read_table(source_path)
    if "tool" in df.columns:
        df["tool"] = df["tool"].map(TOOL_MAP).fillna(df["tool"])
    if "protocol" in df.columns:
        df["protocol"] = df["protocol"].map(PROTOCOL_MAP).fillna(df["protocol"])

    # Match reference notebook outlier clipping.
    if "mape_pas" in df.columns:
        df = df[df["mape_pas"] <= df["mape_pas"].quantile(0.999)]
    if "mape_pas_ct" in df.columns:
        df = df[df["mape_pas_ct"] <= df["mape_pas_ct"].quantile(0.999)]

    out_dir = topic_intermediate_dir(TOPIC)
    out_path = out_dir / OUTPUT_NAME
    write_tsv(df, out_path)
    write_json(
        {
            "topic": TOPIC,
            "source": str(source_path),
            "rows": int(len(df)),
            "columns": list(df.columns),
        },
        out_dir / "metadata.json",
    )
    print(f"Wrote intermediate table: {out_path}")


if __name__ == "__main__":
    main()
