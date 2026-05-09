"""I/O helpers for plotting scripts."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix == ".feather":
        return pd.read_feather(path)
    if suffix == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"Unsupported table suffix: {suffix}")


def write_tsv(df: pd.DataFrame, path: Path) -> Path:
    ensure_dir(path.parent)
    df.to_csv(path, sep="\t", index=False)
    return path


def write_parquet(df: pd.DataFrame, path: Path) -> Path:
    ensure_dir(path.parent)
    df.to_parquet(path, index=False)
    return path


def write_json(data: dict, path: Path) -> Path:
    ensure_dir(path.parent)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
    return path
