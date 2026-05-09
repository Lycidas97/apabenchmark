#!/usr/bin/env python3
"""Render dapars2-pipeline raw recall vs sim F1 scatter with x/y error bars."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

SCRIPT_DIR = Path(__file__).resolve().parent
PLOT_ROOT = SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.constants import TOOL_MAP, TOOL_ORDER_CP
from _shared.io import read_table, write_parquet, write_tsv
from _shared.paths import resolve_project_path, topic_figure_dir, topic_intermediate_dir, topic_root
from _shared.style import apply_reference_style, figsize_from_mm, figsize_from_preset, save_figure, tool_palette

TOPIC = "raw_data_performance"
RAW_INPUT = topic_intermediate_dir(TOPIC) / "raw_dapars2_pipeline_recall_boxplot_points.tsv"
SIM_PREP_INPUT = topic_intermediate_dir("sim_data_performance") / "sim_data_performance_prepared.parquet"
FILTER_CONFIG_PATH = topic_root(TOPIC) / "config" / "recall_gt_filter.json"
PLOT_CONFIG_PATH = topic_root(TOPIC) / "config" / "plot_panels.json"
OUTPUT_NAME = "raw_dapars2_pipeline_recall_vs_sim_f1_scatter_errorbar.pdf"
SCATTER_DATA_NAME = "raw_dapars2_pipeline_recall_vs_sim_f1_scatter_points.tsv"
DEFAULT_DATA_ROOT = resolve_project_path("")
MAIN_ONLY_TOOLS = {"dapars2", "scmapa"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scatter dapars2-pipeline raw recall vs sim F1 with x/y error bars."
    )
    parser.add_argument(
        "--data-root",
        default=str(DEFAULT_DATA_ROOT),
        help="Project root containing data/result/performance (used when sim intermediate is missing).",
    )
    parser.add_argument(
        "--plot-config",
        default=str(PLOT_CONFIG_PATH),
        help="Path to plot_panels.json.",
    )
    parser.add_argument(
        "--filter-config",
        default=str(FILTER_CONFIG_PATH),
        help="Path to recall_gt_filter.json (optional extra pair source for sim F1).",
    )
    return parser.parse_args()


def _load_filter_pairs_from_raw(raw: pd.DataFrame) -> set[tuple[str, str]]:
    if not {"filter_type_1", "filter_type_2"}.issubset(raw.columns):
        return set()
    return {
        (str(f1), str(f2))
        for f1, f2 in raw[["filter_type_1", "filter_type_2"]].drop_duplicates().itertuples(index=False, name=None)
    }


def _load_filter_pairs_from_config(filter_config_path: Path) -> set[tuple[str, str]]:
    if not filter_config_path.exists():
        return set()
    cfg = json.loads(filter_config_path.read_text(encoding="utf-8"))
    out: set[tuple[str, str]] = set()
    for item in cfg.get("filter_pairs", []):
        f1 = str(item.get("filter_type_1", "")).strip()
        f2 = str(item.get("filter_type_2", "")).strip()
        if f1 and f2:
            out.add((f1, f2))
    return out


def load_raw_recall(raw_input: Path) -> pd.DataFrame:
    if not raw_input.exists():
        raise FileNotFoundError(f"Missing raw filtered recall table: {raw_input}")
    raw = read_table(raw_input)
    required = {"tool", "recall"}
    missing = sorted(required - set(raw.columns))
    if missing:
        raise KeyError(f"Missing required raw columns: {missing}")
    raw = raw.copy()
    raw["recall"] = pd.to_numeric(raw["recall"], errors="coerce")
    raw = raw[raw["recall"].notna()].copy()
    return raw


def _load_sim_from_prepared(
    filter_pairs: set[tuple[str, str]],
    sim_prepared_input: Path,
    sim_metric_domain: str,
) -> pd.DataFrame | None:
    if not sim_prepared_input.exists():
        return None
    sim = read_table(sim_prepared_input)
    required = {"metric_domain", "tool", "f1", "filter_type_1", "filter_type_2"}
    missing = sorted(required - set(sim.columns))
    if missing:
        return None
    sim = sim[sim["metric_domain"] == sim_metric_domain].copy()
    if filter_pairs:
        pair_keys = {f"{f1}|{f2}" for (f1, f2) in filter_pairs}
        sim_key = sim["filter_type_1"].astype(str) + "|" + sim["filter_type_2"].astype(str)
        sim = sim[sim_key.isin(pair_keys)].copy()
    sim["f1"] = pd.to_numeric(sim["f1"], errors="coerce")
    sim = sim[sim["f1"].notna()].copy()
    return sim


def _resolve_project_root(data_root_arg: str) -> Path:
    root = Path(data_root_arg).expanduser().resolve()
    if (root / "data" / "result" / "performance").exists():
        return root
    if (root / "result" / "performance").exists():
        return root.parent
    raise FileNotFoundError(
        f"Cannot resolve project root from --data-root: {root}. "
        "Expected <root>/data/result/performance or <data_dir>/result/performance"
    )


def _load_sim_from_performance_files(
    filter_pairs: set[tuple[str, str]],
    data_root_arg: str,
    sim_fallback_subdir: str,
    sim_fallback_glob: str,
) -> pd.DataFrame:
    project_root = _resolve_project_root(data_root_arg)
    perf_dir = project_root / sim_fallback_subdir
    files = sorted(perf_dir.glob(sim_fallback_glob))
    if not files:
        raise FileNotFoundError(f"No sim de_apa tables found under: {perf_dir}")

    keep_cols = ["tool", "f1", "filter_type_1", "filter_type_2", "sample"]
    parts: list[pd.DataFrame] = []
    pair_keys = {f"{f1}|{f2}" for (f1, f2) in filter_pairs} if filter_pairs else set()

    for path in files:
        try:
            cur = pd.read_csv(path, sep="\t", usecols=lambda c: c in keep_cols)
        except Exception:
            continue
        need = {"tool", "f1", "filter_type_1", "filter_type_2"}
        if not need.issubset(cur.columns):
            continue
        if pair_keys:
            key = cur["filter_type_1"].astype(str) + "|" + cur["filter_type_2"].astype(str)
            cur = cur[key.isin(pair_keys)]
        if cur.empty:
            continue
        cur["f1"] = pd.to_numeric(cur["f1"], errors="coerce")
        cur = cur[cur["f1"].notna()]
        if not cur.empty:
            parts.append(cur)

    if not parts:
        raise ValueError("No sim rows matched configured filter pairs.")
    return pd.concat(parts, ignore_index=True)


def _load_sim_from_dapars_parallel_policy(
    filter_pairs: set[tuple[str, str]],
    data_root_arg: str,
    sim_fallback_subdir: str,
) -> pd.DataFrame:
    """Match main-branch sim single-filter policy:
    - dapars2/scmapa use *_de_apa_performance.tsv
    - other tools use *_de_apa_dapars_parallel_performance.tsv
    """
    if not filter_pairs:
        raise ValueError("No filter pairs provided for dapars-parallel sim loading.")

    project_root = _resolve_project_root(data_root_arg)
    perf_dir = project_root / sim_fallback_subdir
    tool_dirs = sorted([p for p in perf_dir.iterdir() if p.is_dir()])
    if not tool_dirs:
        raise FileNotFoundError(f"No tool directories found under: {perf_dir}")

    keep_cols = ["tool", "f1", "filter_type_1", "filter_type_2", "sample"]
    pair_keys = {f"{f1}|{f2}" for (f1, f2) in filter_pairs}
    parts: list[pd.DataFrame] = []

    for tool_dir in tool_dirs:
        tool_raw = tool_dir.name
        if tool_raw in MAIN_ONLY_TOOLS:
            pattern = "*_de_apa_performance.tsv"
        else:
            pattern = "*_de_apa_dapars_parallel_performance.tsv"

        files = sorted(tool_dir.glob(pattern))
        for path in files:
            try:
                cur = pd.read_csv(path, sep="\t", usecols=lambda c: c in keep_cols)
            except Exception:
                continue
            need = {"tool", "f1", "filter_type_1", "filter_type_2"}
            if not need.issubset(cur.columns):
                continue
            key = cur["filter_type_1"].astype(str) + "|" + cur["filter_type_2"].astype(str)
            cur = cur[key.isin(pair_keys)]
            if cur.empty:
                continue
            cur["f1"] = pd.to_numeric(cur["f1"], errors="coerce")
            cur = cur[cur["f1"].notna()]
            if not cur.empty:
                parts.append(cur)

    if not parts:
        raise ValueError("No sim rows matched dapars-parallel policy and configured filter pairs.")
    return pd.concat(parts, ignore_index=True)


def _aggregate_metric_cfg(
    df: pd.DataFrame,
    value_col: str,
    prefix: str,
    groupby_cols: list[str],
    errorbar_stat: str,
) -> pd.DataFrame:
    if errorbar_stat not in {"std", "sem"}:
        raise ValueError(f"Unsupported errorbar_stat={errorbar_stat!r}. Use 'std' or 'sem'.")

    grouped = (
        df.groupby(groupby_cols, as_index=False)[value_col]
        .agg(["mean", "std", "count"])
        .rename(
            columns={
                "mean": f"{prefix}_mean",
                "std": f"{prefix}_sd",
                "count": f"{prefix}_n",
            }
        )
    )
    grouped[f"{prefix}_sd"] = grouped[f"{prefix}_sd"].fillna(0.0)
    if errorbar_stat == "sem":
        n = grouped[f"{prefix}_n"].clip(lower=1)
        grouped[f"{prefix}_sd"] = grouped[f"{prefix}_sd"] / n.pow(0.5)
    return grouped


def _load_or_build_sim_cache(
    *,
    cache_path: Path,
    use_cache: bool,
    refresh_cache: bool,
    filter_pairs: set[tuple[str, str]],
    sim_prepared_input: Path,
    sim_metric_domain: str,
    data_root: str,
    sim_fallback_subdir: str,
    sim_fallback_glob: str,
    sim_source_mode: str,
) -> pd.DataFrame:
    if use_cache and cache_path.exists() and not refresh_cache:
        sim = read_table(cache_path)
        if not sim.empty and "sample" in sim.columns:
            return sim
        if not sim.empty and "sample" not in sim.columns:
            print(f"Sim cache missing sample column, rebuilding: {cache_path}")

    if sim_source_mode == "dapars_parallel_policy":
        sim = _load_sim_from_dapars_parallel_policy(
            filter_pairs,
            data_root,
            sim_fallback_subdir,
        )
    else:
        sim = _load_sim_from_prepared(filter_pairs, sim_prepared_input, sim_metric_domain)
        if sim is None or sim.empty:
            sim = _load_sim_from_performance_files(
                filter_pairs,
                data_root,
                sim_fallback_subdir,
                sim_fallback_glob,
            )

    keep_cols = [c for c in ["tool", "f1", "filter_type_1", "filter_type_2", "sample"] if c in sim.columns]
    sim = sim[keep_cols].copy()
    if use_cache:
        write_parquet(sim, cache_path)
        print(f"Wrote sim cache: {cache_path}")
    return sim


def _resolve_figsize(scatter_cfg: dict, figsize_preset: str) -> tuple[float, float]:
    mm = scatter_cfg.get("figsize_mm")
    if isinstance(mm, list) and len(mm) == 2:
        return figsize_from_mm(float(mm[0]), float(mm[1]))
    return figsize_from_preset(figsize_preset)


def _with_suffix(file_name: str, suffix: str) -> str:
    p = Path(file_name)
    return f"{p.stem}_{suffix}{p.suffix}"


def _species_mask(df: pd.DataFrame, species: str) -> pd.Series | None:
    if "sample" not in df.columns:
        return None
    sample = df["sample"].astype(str).str.lower()
    return sample.str.contains(str(species).lower(), na=False)


def main() -> None:
    args = parse_args()
    apply_reference_style()

    plot_cfg = json.loads(Path(args.plot_config).read_text(encoding="utf-8"))
    scatter_cfg = plot_cfg.get("dapars2_pipeline_scatter", {})

    raw_input = topic_intermediate_dir(TOPIC) / str(scatter_cfg.get("raw_input_name", RAW_INPUT.name))
    f1_prepared_topic = str(scatter_cfg.get("f1_prepared_topic", "sim_data_performance"))
    sim_prepared_input = topic_intermediate_dir(f1_prepared_topic) / str(
        scatter_cfg.get("sim_prepared_input", SIM_PREP_INPUT.name)
    )
    sim_metric_domain = str(scatter_cfg.get("sim_metric_domain", "de_apa_performance"))
    sim_fallback_glob = str(scatter_cfg.get("sim_fallback_glob", "*/*_de_apa_performance.tsv"))
    sim_fallback_subdir = str(scatter_cfg.get("sim_fallback_subdir", "data/result/performance"))
    sim_source_mode = str(scatter_cfg.get("sim_source_mode", "prepared_or_fallback")).strip().lower()
    use_sim_cache = bool(scatter_cfg.get("use_sim_cache", True))
    refresh_sim_cache = bool(scatter_cfg.get("refresh_sim_cache", False))
    sim_cache_name = str(scatter_cfg.get("sim_cache_name", "raw_dapars2_scatter_sim_f1_cache.parquet"))
    sim_cache_path = topic_intermediate_dir(TOPIC) / sim_cache_name

    output_figure = str(scatter_cfg.get("output_figure", OUTPUT_NAME))
    output_points = str(scatter_cfg.get("output_points", SCATTER_DATA_NAME))
    figsize_preset = str(scatter_cfg.get("figsize_preset", "single_column"))
    x_label = str(scatter_cfg.get("x_label", "Raw Recall (DaPars2 pipeline, mean ± SD)"))
    y_label = str(scatter_cfg.get("y_label", "Sim F1 (mean ± SD)"))
    xmin = float(scatter_cfg.get("xmin", 0.0))
    xmax = float(scatter_cfg.get("xmax", 1.02))
    ymin = float(scatter_cfg.get("ymin", 0.0))
    ymax = float(scatter_cfg.get("ymax", 1.02))
    marker_size = float(scatter_cfg.get("marker_size", 4))
    capsize = float(scatter_cfg.get("capsize", 2))
    annotate_tools = bool(scatter_cfg.get("annotate_tools", True))
    show_legend = bool(scatter_cfg.get("show_legend", True))
    legend_loc = str(scatter_cfg.get("legend_loc", "center left"))
    legend_bbox_anchor = scatter_cfg.get("legend_bbox_anchor", [1.02, 0.5])
    legend_fontsize = float(scatter_cfg.get("legend_fontsize", 6))
    legend_title_fontsize = float(scatter_cfg.get("legend_title_fontsize", legend_fontsize))
    legend_title = str(scatter_cfg.get("legend_title", "")).strip()
    legend_ncol = int(scatter_cfg.get("legend_ncol", 1))
    show_correlation = bool(scatter_cfg.get("show_correlation", True))
    correlation_method = str(scatter_cfg.get("correlation_method", "pearson")).strip().lower()
    correlation_xy = scatter_cfg.get("correlation_xy", [0.03, 0.97])
    correlation_fontsize = float(scatter_cfg.get("correlation_fontsize", 7))
    axis_label_fontsize = scatter_cfg.get("axis_label_fontsize")
    tick_label_fontsize = scatter_cfg.get("tick_label_fontsize")
    annotation_fontsize = float(scatter_cfg.get("annotation_fontsize", 6))
    figure_bottom = float(scatter_cfg.get("figure_bottom", 0.18))
    errorbar_stat = str(scatter_cfg.get("errorbar_stat", "std"))
    groupby_cols = [str(v) for v in scatter_cfg.get("groupby", ["tool"])]
    if not groupby_cols:
        groupby_cols = ["tool"]

    raw_all = load_raw_recall(raw_input)
    sim_filter_mode = str(scatter_cfg.get("sim_filter_pairs_mode", "union_raw_and_global")).strip().lower()
    raw_pairs = _load_filter_pairs_from_raw(raw_all)
    global_pairs = _load_filter_pairs_from_config(Path(args.filter_config))
    if sim_filter_mode == "raw_only":
        filter_pairs = raw_pairs
    elif sim_filter_mode == "global_only":
        filter_pairs = global_pairs
    else:
        filter_pairs = raw_pairs | global_pairs

    sim_all = _load_or_build_sim_cache(
        cache_path=sim_cache_path,
        use_cache=use_sim_cache,
        refresh_cache=refresh_sim_cache,
        filter_pairs=filter_pairs,
        sim_prepared_input=sim_prepared_input,
        sim_metric_domain=sim_metric_domain,
        data_root=args.data_root,
        sim_fallback_subdir=sim_fallback_subdir,
        sim_fallback_glob=sim_fallback_glob,
        sim_source_mode=sim_source_mode,
    )
    if filter_pairs:
        sim_pair_keys = {f"{f1}|{f2}" for (f1, f2) in filter_pairs}
        sim_key = sim_all["filter_type_1"].astype(str) + "|" + sim_all["filter_type_2"].astype(str)
        sim_all = sim_all[sim_key.isin(sim_pair_keys)].copy()
    if sim_all.empty:
        raise ValueError("Sim cache/source has no rows for dapars2 pipeline filter pairs.")

    variants = [None, "human", "mouse"]
    for species in variants:
        raw = raw_all
        sim = sim_all
        if species is not None:
            raw_mask = _species_mask(raw_all, species)
            sim_mask = _species_mask(sim_all, species)
            if raw_mask is None or sim_mask is None:
                print(f"Skip species split ({species}): missing 'sample' column in raw/sim.")
                continue
            raw = raw_all[raw_mask].copy()
            sim = sim_all[sim_mask].copy()
            if raw.empty or sim.empty:
                print(f"Skip species split ({species}): no matched rows.")
                continue

        raw_stats = _aggregate_metric_cfg(raw, "recall", "raw_recall", groupby_cols, errorbar_stat)
        sim_stats = _aggregate_metric_cfg(sim, "f1", "sim_f1", groupby_cols, errorbar_stat)

        raw_stats["tool"] = raw_stats["tool"].map(TOOL_MAP).fillna(raw_stats["tool"])
        sim_stats["tool"] = sim_stats["tool"].map(TOOL_MAP).fillna(sim_stats["tool"])

        merge_cols = [col for col in groupby_cols if col in raw_stats.columns and col in sim_stats.columns]
        if not merge_cols:
            raise ValueError("No shared groupby columns found between raw and sim aggregates.")
        points = pd.merge(raw_stats, sim_stats, on=merge_cols, how="inner")
        if points.empty:
            print(f"Skip variant ({species or 'all'}): no overlapping groups after aggregation.")
            continue

        points_name = output_points if species is None else _with_suffix(output_points, species)
        out_points = topic_intermediate_dir(TOPIC) / points_name
        write_tsv(points.sort_values("tool"), out_points)

        present_tools = set(points["tool"].astype(str))
        legend_order = [t for t in TOOL_ORDER_CP if t in present_tools]
        legend_order.extend(sorted(t for t in present_tools if t not in legend_order))
        points["tool"] = pd.Categorical(points["tool"], categories=legend_order, ordered=True)
        points = points.sort_values("tool").reset_index(drop=True)

        fig, ax = plt.subplots(figsize=_resolve_figsize(scatter_cfg, figsize_preset))
        fig.subplots_adjust(bottom=figure_bottom)
        colors = tool_palette(legend_order)

        for _, row in points.iterrows():
            tool = str(row["tool"])
            x = float(row["raw_recall_mean"])
            y = float(row["sim_f1_mean"])
            xerr = float(row["raw_recall_sd"])
            yerr = float(row["sim_f1_sd"])
            c = colors.get(tool, "#929292")
            ax.errorbar(
                x,
                y,
                xerr=xerr,
                yerr=yerr,
                fmt="o",
                color=c,
                ecolor=c,
                elinewidth=0.8,
                capsize=capsize,
                markersize=marker_size,
                alpha=0.95,
                label=tool,
            )
            if annotate_tools:
                ax.annotate(
                    tool,
                    (x, y),
                    xytext=(4, 3),
                    textcoords="offset points",
                    fontsize=annotation_fontsize,
                )

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_xlabel(x_label, fontsize=axis_label_fontsize)
        ax.set_ylabel(y_label, fontsize=axis_label_fontsize)
        if tick_label_fontsize is not None:
            ax.tick_params(axis="both", labelsize=float(tick_label_fontsize))
        ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.7)

        if show_legend:
            if not isinstance(legend_bbox_anchor, list) or len(legend_bbox_anchor) != 2:
                legend_bbox = [1.02, 0.5]
            else:
                legend_bbox = legend_bbox_anchor
            handles = [
                Line2D([0], [0], marker="o", linestyle="", color=colors.get(tool, "#929292"), label=tool, markersize=4)
                for tool in legend_order
            ]
            fig.legend(
                handles=handles,
                labels=legend_order,
                loc=legend_loc,
                bbox_to_anchor=(float(legend_bbox[0]), float(legend_bbox[1])),
                frameon=False,
                fontsize=legend_fontsize,
                title=(legend_title if legend_title else None),
                title_fontsize=legend_title_fontsize,
                ncol=max(1, legend_ncol),
            )

        if show_correlation and len(points) >= 2:
            corr_method = correlation_method
            if corr_method not in {"pearson", "spearman"}:
                corr_method = "pearson"
            x_ser = pd.to_numeric(points["raw_recall_mean"], errors="coerce")
            y_ser = pd.to_numeric(points["sim_f1_mean"], errors="coerce")
            valid = x_ser.notna() & y_ser.notna()
            if int(valid.sum()) >= 2:
                r = x_ser[valid].corr(y_ser[valid], method=corr_method)
                if not np.isnan(r):
                    text = f"{corr_method.title()} r = {r:.3f}"
                    if not isinstance(correlation_xy, list) or len(correlation_xy) != 2:
                        corr_xy = [0.03, 0.97]
                    else:
                        corr_xy = correlation_xy
                    ax.text(
                        float(corr_xy[0]),
                        float(corr_xy[1]),
                        text,
                        transform=ax.transAxes,
                        fontsize=correlation_fontsize,
                        va="top",
                        ha="left",
                        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.7, "pad": 1.5},
                    )

        fig_name = output_figure if species is None else _with_suffix(output_figure, species)
        out_fig = topic_figure_dir(TOPIC) / fig_name
        save_figure(fig, out_fig)
        plt.close(fig)

        print(f"Wrote scatter data: {out_points}")
        print(f"Wrote figure: {out_fig}")


if __name__ == "__main__":
    main()
