#!/usr/bin/env python3
"""Prepare intermediate table for peak model parameter plots."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import PROTOCOL_MAP
from _shared.io import write_json, write_tsv
from _shared.paths import resolve_project_path, topic_intermediate_dir

TOPIC = "peak_model_params"
OUTPUT_NAME = "peak_model_params_prepared.tsv"
SOURCE_CANDIDATES = [
    "data/int_data/model_result",
    "data/model_result",
]
MODEL_SUFFIX = "_model_results"
SHAPE_SUMMARY_SUFFIX = "_peak_shape_summary"
SHAPE_BYPAS_SUFFIX = "_peak_shape_by_pas"
KEY_COLS = ["sample_full_id", "sample_id", "protocol", "species", "tissue"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare intermediate table for peak model params plotting")
    parser.add_argument(
        "--data-root",
        type=Path,
        default=None,
        help="Project root containing data/ or the data directory itself.",
    )
    parser.add_argument(
        "--source-dir",
        type=Path,
        default=None,
        help="Directory containing *_model_results.json and optional *_peak_shape_*.feather files",
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


def pick_source_dir(source_dir: Path | None, data_root: Path | None = None) -> Path:
    if source_dir is not None:
        path = source_dir.resolve()
        if not path.exists():
            raise FileNotFoundError(f"Missing input directory: {path}")
        return path

    tried: list[Path] = []
    root = resolve_data_root(data_root)
    if root is not None:
        for relative in SOURCE_CANDIDATES:
            path = root / relative
            tried.append(path)
            if path.exists():
                return path

    for relative in SOURCE_CANDIDATES:
        path = resolve_project_path(relative)
        tried.append(path)
        if path.exists():
            return path

    tried_msg = "\n".join(f"- {p}" for p in tried)
    raise FileNotFoundError(
        "Could not find default input directory for peak model params.\n"
        f"Tried:\n{tried_msg}\n"
        "Use --source-dir to specify the model_result directory."
    )


def strip_suffix(stem: str, suffix: str) -> str:
    if stem.endswith(suffix):
        return stem[: -len(suffix)]
    return stem


def sample_ids_from_full(full_sample_id: str) -> tuple[str, str, str, str]:
    parts = full_sample_id.split("_")
    sample_id = "_".join(parts[:4]) if len(parts) >= 4 else full_sample_id
    protocol = parts[0] if len(parts) > 0 else np.nan
    species = parts[1] if len(parts) > 1 else np.nan
    tissue = parts[2] if len(parts) > 2 else np.nan
    return sample_id, protocol, species, tissue


def empty_with_keys(extra_cols: list[str] | None = None) -> pd.DataFrame:
    cols = KEY_COLS.copy()
    if extra_cols:
        cols.extend(extra_cols)
    return pd.DataFrame(columns=cols)


def build_model_df(source_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for path in sorted(source_dir.glob("*_model_results.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        full_sample_id = strip_suffix(path.stem, MODEL_SUFFIX)
        sample_id, protocol, species, tissue = sample_ids_from_full(full_sample_id)

        custom = payload.get("custom_model", {})
        normal = payload.get("normal_dist", {})
        custom_params = custom.get("params", {})
        normal_params = normal.get("params", {})

        rows.append(
            {
                "sample_full_id": full_sample_id,
                "sample_id": sample_id,
                "protocol": protocol,
                "species": species,
                "tissue": tissue,
                "model_result_path": str(path),
                "custom_mu1": custom_params.get("mu1", np.nan),
                "custom_sigma1": custom_params.get("sigma1", np.nan),
                "custom_sigma2": custom_params.get("sigma2", np.nan),
                "custom_a": custom_params.get("a", np.nan),
                "custom_b": custom_params.get("b", np.nan),
                "custom_log_likelihood": custom.get("log_likelihood", np.nan),
                "custom_AIC": custom.get("AIC", np.nan),
                "custom_BIC": custom.get("BIC", np.nan),
                "normal_mu": normal_params.get("mu", np.nan),
                "normal_sigma": normal_params.get("sigma", np.nan),
                "normal_log_likelihood": normal.get("log_likelihood", np.nan),
                "normal_AIC": normal.get("AIC", np.nan),
                "normal_BIC": normal.get("BIC", np.nan),
                "data_size": payload.get("data_size", np.nan),
                "apex_mean": payload.get("apex_mean", np.nan),
                "apex_std": payload.get("apex_std", np.nan),
                "weighted_mean": payload.get("weighted_mean", np.nan),
                "weighted_std": payload.get("weighted_std", np.nan),
            }
        )

    if not rows:
        return empty_with_keys(
            [
                "model_result_path",
                "custom_mu1",
                "custom_sigma1",
                "custom_sigma2",
                "custom_a",
                "custom_b",
                "custom_log_likelihood",
                "custom_AIC",
                "custom_BIC",
                "normal_mu",
                "normal_sigma",
                "normal_log_likelihood",
                "normal_AIC",
                "normal_BIC",
                "data_size",
                "apex_mean",
                "apex_std",
                "weighted_mean",
                "weighted_std",
            ]
        )

    df = pd.DataFrame(rows)
    df["delta_AIC_custom"] = df["custom_AIC"] - df[["custom_AIC", "normal_AIC"]].min(axis=1)
    df["delta_AIC_normal"] = df["normal_AIC"] - df[["custom_AIC", "normal_AIC"]].min(axis=1)
    df["delta_AIC"] = df["custom_AIC"] - df["normal_AIC"]
    df["delta_BIC"] = df["custom_BIC"] - df["normal_BIC"]
    return df


def summarize_by_pas(df: pd.DataFrame) -> dict[str, Any]:
    row: dict[str, Any] = {"shape_by_pas_n_rows": int(len(df))}
    numeric = df.select_dtypes(include=[np.number])
    for col in numeric.columns:
        series = numeric[col].dropna()
        if series.empty:
            continue
        row[f"shape_by_pas_{col}_mean"] = float(series.mean())
        row[f"shape_by_pas_{col}_median"] = float(series.median())
        row[f"shape_by_pas_{col}_std"] = float(series.std(ddof=1)) if len(series) > 1 else np.nan
    return row


def build_shape_summary_df(source_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for path in sorted(source_dir.glob("*_peak_shape_summary.feather")):
        frame = pd.read_feather(path)
        if frame.empty:
            continue
        full_sample_id = strip_suffix(path.stem, SHAPE_SUMMARY_SUFFIX)
        sample_id, protocol, species, tissue = sample_ids_from_full(full_sample_id)
        row = {f"shape_{k}": v for k, v in frame.iloc[0].to_dict().items()}
        row.update(
            {
                "sample_full_id": full_sample_id,
                "sample_id": sample_id,
                "protocol": protocol,
                "species": species,
                "tissue": tissue,
                "shape_summary_path": str(path),
            }
        )
        rows.append(row)

    if not rows:
        return empty_with_keys(["shape_summary_path"])
    return pd.DataFrame(rows)


def build_shape_by_pas_summary_df(source_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for path in sorted(source_dir.glob("*_peak_shape_by_pas.feather")):
        frame = pd.read_feather(path)
        if frame.empty:
            continue
        full_sample_id = strip_suffix(path.stem, SHAPE_BYPAS_SUFFIX)
        sample_id, protocol, species, tissue = sample_ids_from_full(full_sample_id)
        row = summarize_by_pas(frame)
        row.update(
            {
                "sample_full_id": full_sample_id,
                "sample_id": sample_id,
                "protocol": protocol,
                "species": species,
                "tissue": tissue,
                "shape_by_pas_path": str(path),
            }
        )
        rows.append(row)

    if not rows:
        return empty_with_keys(["shape_by_pas_path", "shape_by_pas_n_rows"])
    return pd.DataFrame(rows)


def merge_tables(model_df: pd.DataFrame, shape_summary_df: pd.DataFrame, shape_by_pas_df: pd.DataFrame) -> pd.DataFrame:
    merged = pd.merge(model_df, shape_summary_df, on=KEY_COLS, how="outer")
    merged = pd.merge(merged, shape_by_pas_df, on=KEY_COLS, how="outer")
    merged["protocol"] = merged["protocol"].map(PROTOCOL_MAP).fillna(merged["protocol"])
    return merged.sort_values(["sample_id", "sample_full_id"], kind="stable").reset_index(drop=True)


def main() -> None:
    args = parse_args()
    source_dir = pick_source_dir(args.source_dir, args.data_root)

    model_df = build_model_df(source_dir)
    shape_summary_df = build_shape_summary_df(source_dir)
    shape_by_pas_df = build_shape_by_pas_summary_df(source_dir)
    merged = merge_tables(model_df, shape_summary_df, shape_by_pas_df)

    out_dir = topic_intermediate_dir(TOPIC)
    out_path = out_dir / OUTPUT_NAME
    write_tsv(merged, out_path)
    write_json(
        {
            "topic": TOPIC,
            "source_dir": str(source_dir),
            "rows": int(len(merged)),
            "columns": list(merged.columns),
            "n_model_rows": int(len(model_df)),
            "n_shape_summary_rows": int(len(shape_summary_df)),
            "n_shape_by_pas_rows": int(len(shape_by_pas_df)),
            "model_paths": int(len(list(source_dir.glob("*_model_results.json")))),
            "shape_summary_paths": int(len(list(source_dir.glob("*_peak_shape_summary.feather")))),
            "shape_by_pas_paths": int(len(list(source_dir.glob("*_peak_shape_by_pas.feather")))),
        },
        out_dir / "metadata.json",
    )
    print(f"Wrote intermediate table: {out_path}")


if __name__ == "__main__":
    main()
