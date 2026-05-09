"""Shared helpers for sim_data_performance panel scripts."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib as mpl
from matplotlib.patches import Patch
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import PALETTE_HEX, PROTOCOL_MAP, PROTOCOL_ORDER, TOOL_ORDER
from _shared.io import read_table
from _shared.paths import topic_intermediate_dir

TOPIC = "sim_data_performance"
MAIN_INPUT = "sim_data_performance_prepared.parquet"
APA_RESIDUAL_INPUT = "sim_data_apa_residuals.parquet"
APA_RESIDUAL_BY_FILTER_INPUT = "sim_data_apa_residuals_by_filter.parquet"
SINGLE_FILTER_APA_INPUT = "sim_data_differential_apa_single_filter_performance.parquet"
HATCH_PATTERN = "///////"
BOXPLOT_TRIM_LOWER_Q = 0.01
BOXPLOT_TRIM_UPPER_Q = 0.99

MATCH_TYPE_ORDER = [
    "pas_25",
    "pas_50",
    "pas_75",
    "pas_100",
    "gt_25",
    "gt_50",
    "gt_75",
    "gt_100",
    "pd_25",
    "pd_50",
    "pd_75",
    "pd_100",
    "te",
]
PD_MATCH_TYPES = {"pd_25", "pd_50", "pd_75", "pd_100"}

REFERENCE_APA_FILTER_SET = [
    "wilcox_DWUI_0.05|dexseq_log2fc_0.5",
    "wilcox_DWUI_0.05|dexseq_log2fc_1",
    "wilcox_DWUI_0.05|dexseq_log2fc_1.25",
]

TOOL_COLOR_OVERRIDES = {
    "DaPars2": "#3b9ab2",
    "scMAPA": "#c76d9e",
}


def _fallback_sorted(values: list[str]) -> list[str]:
    return sorted(values)


def resolve_match_type_order(values: list[str]) -> list[str]:
    observed = [v for v in values if isinstance(v, str)]
    in_order = [v for v in MATCH_TYPE_ORDER if v in observed]
    leftovers = [v for v in observed if v not in in_order]
    seen = set()
    unique_leftovers = []
    for item in _fallback_sorted(leftovers):
        if item not in seen:
            seen.add(item)
            unique_leftovers.append(item)
    return in_order + unique_leftovers


def match_type_palette(hue_order: list[str]) -> dict[str, str]:
    palette_map = {
        "pas_25": PALETTE_HEX[0],
        "pas_50": PALETTE_HEX[1],
        "pas_75": PALETTE_HEX[2],
        "pas_100": PALETTE_HEX[3],
        "gt_25": PALETTE_HEX[0],
        "gt_50": PALETTE_HEX[1],
        "gt_75": PALETTE_HEX[2],
        "gt_100": PALETTE_HEX[3],
        "pd_25": PALETTE_HEX[0],
        "pd_50": PALETTE_HEX[1],
        "pd_75": PALETTE_HEX[2],
        "pd_100": PALETTE_HEX[3],
        "te": PALETTE_HEX[4],
    }
    return {name: palette_map.get(name, PALETTE_HEX[-1]) for name in hue_order}


def build_hatched_legend_handles(hue_order: list[str]) -> list[Patch]:
    color_map = match_type_palette(hue_order)
    handles: list[Patch] = []
    for name in hue_order:
        patch = Patch(facecolor=color_map[name], edgecolor="black", label=name)
        if name in PD_MATCH_TYPES:
            patch.set_hatch(HATCH_PATTERN)
        handles.append(patch)
    return handles


def apply_pd_hatching(ax, hue_order: list[str]) -> None:
    if not hue_order:
        return
    n_hues = len(hue_order)
    patch_with_x: list[tuple[float, mpl.patches.PathPatch]] = []
    for child in ax.get_children():
        if isinstance(child, mpl.patches.PathPatch):
            verts = child.get_path().vertices
            if verts is None or len(verts) < 5:
                continue
            patch_with_x.append((float(verts[:, 0].min()), child))
    patch_with_x.sort(key=lambda item: item[0])
    for idx, (_, patch) in enumerate(patch_with_x):
        hue_name = hue_order[idx % n_hues]
        if hue_name in PD_MATCH_TYPES:
            patch.set_hatch(HATCH_PATTERN)
            patch.set_edgecolor("black")
            patch.set_linewidth(0.6)


def black_boxplot_props(linewidth: float = 0.5) -> dict[str, object]:
    """Return seaborn/matplotlib boxplot line-style kwargs with pure black outlines."""
    lw = float(linewidth)
    return {
        "linecolor": "black",
        "boxprops": {"edgecolor": "black", "linewidth": lw},
        "whiskerprops": {"color": "black", "linewidth": lw},
        "capprops": {"color": "black", "linewidth": lw},
        "medianprops": {"color": "black", "linewidth": lw},
        "flierprops": {
            "markeredgecolor": "black",
            "markerfacecolor": "white",
            "markersize": 2.0,
            "markeredgewidth": max(0.4, lw),
        },
    }


def load_main_prepared() -> pd.DataFrame:
    path = topic_intermediate_dir(TOPIC) / MAIN_INPUT
    if not path.exists():
        raise FileNotFoundError(
            f"Missing prepared table: {path}. Run 00_prepare_data.py first."
        )
    return read_table(path)


def load_metric_domain(metric_domain: str) -> pd.DataFrame:
    df = load_main_prepared()
    if "metric_domain" not in df.columns:
        raise KeyError("Column 'metric_domain' is required in prepared table.")
    out = df[df["metric_domain"] == metric_domain].copy()
    if out.empty:
        raise ValueError(f"No rows found for metric_domain={metric_domain!r}")
    return out


def load_apa_residual_table() -> pd.DataFrame:
    path = topic_intermediate_dir(TOPIC) / APA_RESIDUAL_INPUT
    if not path.exists():
        raise FileNotFoundError(
            f"Missing APA residual table: {path}. Run 02_prepare_apa_residuals.py first."
        )
    return read_table(path)


def load_apa_residual_by_filter_table() -> pd.DataFrame:
    path = topic_intermediate_dir(TOPIC) / APA_RESIDUAL_BY_FILTER_INPUT
    if not path.exists():
        raise FileNotFoundError(
            f"Missing APA residual-by-filter table: {path}. Run 02_prepare_apa_residuals.py first."
        )
    return read_table(path)


def load_single_filter_apa_table() -> pd.DataFrame:
    path = topic_intermediate_dir(TOPIC) / SINGLE_FILTER_APA_INPUT
    if not path.exists():
        raise FileNotFoundError(
            f"Missing single-filter APA table: {path}. "
            "Run 03_prepare_differential_apa_single_filter_performance.py first."
        )
    return read_table(path)


def _prepare_filter_comb(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "filter_type_comb" not in out.columns and {"filter_type_1", "filter_type_2"}.issubset(out.columns):
        out["filter_type_comb"] = (
            out["filter_type_1"].astype(str) + "|" + out["filter_type_2"].astype(str)
        )
    return out


def _pick_apa_filter_set(df: pd.DataFrame) -> list[str]:
    required = {"filter_type_comb", "f1_residuals", "gt_recall"}
    if not required.issubset(df.columns):
        return []
    clean = df[list(required)].dropna().copy()
    if clean.empty:
        return []

    by_filter = (
        clean.groupby("filter_type_comb", dropna=False)[["f1_residuals", "gt_recall"]]
        .mean(numeric_only=True)
        .reset_index()
    )
    by_filter["filter_type_comb"] = by_filter["filter_type_comb"].astype(str)

    y_min = float(by_filter["gt_recall"].min())
    y_max = float(by_filter["gt_recall"].max())
    third_1 = y_min + (y_max - y_min) / 3.0
    third_2 = y_min + 2.0 * (y_max - y_min) / 3.0

    sections = [
        by_filter[by_filter["gt_recall"] <= third_1],
        by_filter[(by_filter["gt_recall"] > third_1) & (by_filter["gt_recall"] <= third_2)],
        by_filter[by_filter["gt_recall"] > third_2],
    ]

    selected: list[str] = []
    for section in sections:
        if section.empty:
            continue
        # Keep the same selection rule as residual scatter markers (true top1).
        combo = str(section.nlargest(1, "f1_residuals").iloc[0]["filter_type_comb"])
        if combo not in selected:
            selected.append(combo)

    if selected:
        return selected[:3]

    available = set(clean["filter_type_comb"].astype(str).unique())
    preferred = [name for name in REFERENCE_APA_FILTER_SET if name in available]
    if preferred:
        return preferred
    counts = clean["filter_type_comb"].astype(str).value_counts().sort_values(ascending=False)
    return list(counts.head(3).index)


def build_apa_plot_table() -> pd.DataFrame:
    df = build_apa_plot_table_by_filter()
    metric_cols = [name for name in ["precision", "recall", "f1"] if name in df.columns]
    grouped = (
        df.groupby(["tool", "sample", "match_type"], as_index=False)[metric_cols]
        .mean(numeric_only=True)
    )
    grouped["protocol"] = grouped["sample"].astype(str).str.split("_").str[0]
    grouped["protocol"] = grouped["protocol"].map(PROTOCOL_MAP).fillna(grouped["protocol"])
    return grouped


def build_apa_plot_table_by_filter() -> pd.DataFrame:
    """Return APA residual table constrained to the selected filter combinations."""
    df = _prepare_filter_comb(load_apa_residual_table())
    filter_set = _pick_apa_filter_set(df)
    if filter_set:
        df = df[df["filter_type_comb"].isin(filter_set)].copy()
    if df.empty:
        raise ValueError("APA residual table is empty after filter selection.")

    required = {"tool", "sample", "match_type", "filter_type_comb", "precision", "recall", "f1"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise KeyError(f"Missing APA columns for filter-split plotting: {missing}")

    out = df[list(required)].copy()
    out["protocol"] = out["sample"].astype(str).str.split("_").str[0]
    out["protocol"] = out["protocol"].map(PROTOCOL_MAP).fillna(out["protocol"])
    return out


def available_tool_order(df: pd.DataFrame) -> list[str]:
    present = set(df["tool"].dropna().astype(str).unique()) if "tool" in df.columns else set()
    ordered = [name for name in TOOL_ORDER if name in present]
    leftovers = sorted(name for name in present if name not in ordered)
    return ordered + leftovers


def available_protocol_order(df: pd.DataFrame) -> list[str]:
    present = set(df["protocol"].dropna().astype(str).unique()) if "protocol" in df.columns else set()
    ordered = [name for name in PROTOCOL_ORDER if name in present]
    leftovers = sorted(name for name in present if name not in ordered)
    return ordered + leftovers


def tool_palette(tool_order: list[str]) -> dict[str, str]:
    palette: dict[str, str] = {}
    for idx, name in enumerate(tool_order):
        palette[name] = TOOL_COLOR_OVERRIDES.get(name, PALETTE_HEX[idx % len(PALETTE_HEX)])
    return palette


def format_filter_type_comb(text: str) -> str:
    pieces = text.split("|")
    out_parts = []
    for part in pieces:
        if "_" in part:
            name, value = part.rsplit("_", 1)
            out_parts.append(f"{name} ({value})")
        else:
            out_parts.append(part)
    return " | ".join(out_parts)


def trim_metric_to_quantile(
    df: pd.DataFrame,
    metric: str,
    *,
    lower_q: float = BOXPLOT_TRIM_LOWER_Q,
    upper_q: float = BOXPLOT_TRIM_UPPER_Q,
) -> pd.DataFrame:
    """Keep rows where metric is within [q(lower_q), q(upper_q)]."""
    if metric not in df.columns:
        return df.iloc[0:0].copy()

    out = df.copy()
    out[metric] = pd.to_numeric(out[metric], errors="coerce")
    out = out[out[metric].notna()].copy()
    if out.empty:
        return out

    low = float(out[metric].quantile(lower_q))
    high = float(out[metric].quantile(upper_q))
    if pd.isna(low) or pd.isna(high):
        return out.iloc[0:0].copy()
    if low > high:
        low, high = high, low
    return out[out[metric].between(low, high, inclusive="both")].copy()
