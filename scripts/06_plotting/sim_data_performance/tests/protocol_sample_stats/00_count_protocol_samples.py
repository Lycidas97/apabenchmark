#!/usr/bin/env python3
"""Count unique sim-data samples by protocol for plotting QA."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[2]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import PROTOCOL_MAP
from _shared.paths import topic_intermediate_dir

TOPIC = "sim_data_performance"
DEFAULT_INPUT_NAME = "sim_data_performance_prepared.parquet"
DEFAULT_METADATA_NAME = "metadata.json"

SINGLE_CELL_PROTOCOLS = {
    "10X Chromium",
    "Drop-seq",
    "Microwell-seq",
}

SPATIAL_PROTOCOLS = {
    "10X Visium",
    "10X Visium HD",
    "Stereo-seq",
    "Slide-seq V2",
    "ST",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize protocol/sample usage from sim_data_performance prepared parquet."
    )
    parser.add_argument(
        "--input-parquet",
        type=Path,
        default=topic_intermediate_dir(TOPIC) / DEFAULT_INPUT_NAME,
        help="Input prepared parquet path.",
    )
    parser.add_argument(
        "--metadata-json",
        type=Path,
        default=topic_intermediate_dir(TOPIC) / DEFAULT_METADATA_NAME,
        help="Metadata json produced by 00_prepare_data.py. Used as fast source of sample list.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SCRIPT_DIR / "output",
        help="Output directory for summary tables.",
    )
    return parser.parse_args()


def normalize_protocol(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["sample"] = out["sample"].astype(str).str.strip()
    out["protocol"] = out["protocol"].astype(str).str.strip()
    out["protocol"] = out["protocol"].replace({"": pd.NA, "nan": pd.NA, "None": pd.NA})

    inferred = out["sample"].str.split("_").str[0]
    out["protocol"] = out["protocol"].fillna(inferred)
    out["protocol"] = out["protocol"].map(PROTOCOL_MAP).fillna(out["protocol"])
    return out


def build_sample_protocol(df: pd.DataFrame) -> pd.DataFrame:
    required = {"sample", "protocol"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns in input: {sorted(missing)}")

    pairs = normalize_protocol(df[list(required)]).dropna(subset=["sample", "protocol"]).copy()
    pairs = pairs[pairs["sample"] != ""]

    # Each sample should map to one protocol; keep first if repeated.
    pairs = pairs.sort_values(["sample", "protocol"])
    sample_protocol = pairs.drop_duplicates(subset=["sample"], keep="first").reset_index(drop=True)
    return sample_protocol


def sample_protocol_from_metadata(metadata_json: Path) -> pd.DataFrame:
    if not metadata_json.exists():
        raise FileNotFoundError(f"Metadata json not found: {metadata_json}")

    payload = json.loads(metadata_json.read_text(encoding="utf-8"))
    source_files_by_domain = payload.get("source_files_by_domain", {})
    if not isinstance(source_files_by_domain, dict):
        raise KeyError("Invalid metadata: source_files_by_domain must be a dict.")

    domain_suffix = {
        "match_performance": "_match_performance.tsv",
        "pas_quantify_performance": "_pas_quantify_performance.tsv",
        "de_apa_performance": "_de_apa_performance.tsv",
    }

    rows: list[dict[str, str]] = []
    for metric_domain, paths in source_files_by_domain.items():
        suffix = domain_suffix.get(metric_domain)
        if suffix is None or not isinstance(paths, list):
            continue
        for rel_path in paths:
            name = Path(str(rel_path)).name
            if not name.endswith(suffix):
                continue
            sample = name[: -len(suffix)]
            protocol = sample.split("_")[0] if sample else ""
            rows.append({"sample": sample, "protocol": protocol})

    if not rows:
        raise ValueError("No sample entries extracted from metadata source_files_by_domain.")

    return build_sample_protocol(pd.DataFrame(rows))


def classify_protocol(protocol: str) -> str:
    if protocol in SINGLE_CELL_PROTOCOLS:
        return "single_cell"
    if protocol in SPATIAL_PROTOCOLS:
        return "spatial"
    return "unknown"


def write_outputs(sample_protocol: pd.DataFrame, output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)

    protocol_counts = (
        sample_protocol.groupby("protocol", as_index=False)
        .size()
        .rename(columns={"size": "sample_count"})
        .sort_values(["sample_count", "protocol"], ascending=[False, True])
        .reset_index(drop=True)
    )
    protocol_counts["category"] = protocol_counts["protocol"].map(classify_protocol)

    total_samples = int(sample_protocol["sample"].nunique())
    single_cell_samples = int(
        sample_protocol[sample_protocol["protocol"].isin(SINGLE_CELL_PROTOCOLS)]["sample"].nunique()
    )
    spatial_samples = int(
        sample_protocol[sample_protocol["protocol"].isin(SPATIAL_PROTOCOLS)]["sample"].nunique()
    )
    unknown_samples = total_samples - single_cell_samples - spatial_samples

    protocol_counts.to_parquet(output_dir / "protocol_sample_counts.parquet", index=False)
    sample_protocol.to_parquet(output_dir / "sample_protocol_map.parquet", index=False)

    summary = {
        "total_unique_samples": total_samples,
        "single_cell_unique_samples": single_cell_samples,
        "spatial_unique_samples": spatial_samples,
        "unknown_unique_samples": unknown_samples,
        "single_cell_protocols": sorted(SINGLE_CELL_PROTOCOLS),
        "spatial_protocols": sorted(SPATIAL_PROTOCOLS),
        "protocols_seen": protocol_counts["protocol"].tolist(),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True),
        encoding="utf-8",
    )

    lines = [
        "# sim_data protocol/sample summary",
        "",
        f"- total_unique_samples: {total_samples}",
        f"- single_cell_unique_samples: {single_cell_samples}",
        f"- spatial_unique_samples: {spatial_samples}",
        f"- unknown_unique_samples: {unknown_samples}",
        "",
        "## protocol sample counts",
        "",
        "| protocol | category | sample_count |",
        "|---|---|---:|",
    ]
    for _, row in protocol_counts.iterrows():
        lines.append(f"| {row['protocol']} | {row['category']} | {int(row['sample_count'])} |")
    (output_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    return summary


def main() -> None:
    args = parse_args()
    sample_protocol: pd.DataFrame
    source_used: str

    if args.metadata_json.exists():
        sample_protocol = sample_protocol_from_metadata(args.metadata_json)
        source_used = "metadata_json"
    else:
        if not args.input_parquet.exists():
            raise FileNotFoundError(
                f"Neither metadata nor input parquet exists:\n"
                f"  metadata_json={args.metadata_json}\n"
                f"  input_parquet={args.input_parquet}"
            )
        df = pd.read_parquet(args.input_parquet, columns=["sample", "protocol"])
        sample_protocol = build_sample_protocol(df)
        source_used = "input_parquet"

    summary = write_outputs(sample_protocol, args.output_dir)

    print(f"source_used={source_used}")
    print(f"metadata_json={args.metadata_json}")
    print(f"input={args.input_parquet}")
    print(f"output_dir={args.output_dir}")
    print(f"total_unique_samples={summary['total_unique_samples']}")
    print(f"single_cell_unique_samples={summary['single_cell_unique_samples']}")
    print(f"spatial_unique_samples={summary['spatial_unique_samples']}")
    print(f"unknown_unique_samples={summary['unknown_unique_samples']}")


if __name__ == "__main__":
    main()
