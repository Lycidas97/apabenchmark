#!/usr/bin/env python3
"""Run sim data performance plotting workflow with intermediate-file reuse."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
from typing import Iterable

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_PLOT_CONFIG = SCRIPT_DIR.parent / "config" / "plot_params.json"

PREPARE_SCRIPTS = [
    "00_prepare_data.py",
    "01_prepare_dapars2_pipeline_data.py",
    "02_prepare_apa_residuals.py",
    "03_prepare_differential_apa_single_filter_performance.py",
    "04_prepare_pairwise_winrate.py",
]

PANEL_SCRIPTS = [
    "10_panel_pas_detect_overall.py",
    "11_panel_pas_detect_by_protocol.py",
    "20_panel_pas_quantification_overall.py",
    "21_panel_pas_quantification_by_protocol.py",
    "30_panel_apa_detect_overall.py",
    "31_panel_apa_detect_by_protocol.py",
    "32_panel_apa_detect_single_filter_overall.py",
    "33_panel_apa_detect_single_filter_by_protocol.py",
    "40_panel_filter_residuals_heatmap.py",
    "41_panel_max_residuals_vs_gt_recall.py",
    "42_panel_f1_vs_stringency_regression.py",
    "50_panel_pairwise_winrate_heatmaps.py",
]

INTERMEDIATE_DIR = SCRIPT_DIR.parent / "data" / "intermediate"

PREPARE_OUTPUTS = {
    "00_prepare_data.py": [
        "sim_data_performance_prepared.parquet",
    ],
    "01_prepare_dapars2_pipeline_data.py": [
        "sim_data_dapars2_pipeline_manifest.parquet",
        "sim_data_dapars2_pipeline_summary.parquet",
    ],
    "02_prepare_apa_residuals.py": [
        "sim_data_apa_residuals.parquet",
        "sim_data_apa_residuals_by_filter.parquet",
    ],
    "03_prepare_differential_apa_single_filter_performance.py": [
        "sim_data_differential_apa_single_filter_performance.parquet",
    ],
    "04_prepare_pairwise_winrate.py": [
        "sim_data_pairwise_winrate.parquet",
    ],
}

DEFAULT_REQUIRED_PREPARE = [
    "00_prepare_data.py",
    "02_prepare_apa_residuals.py",
    "03_prepare_differential_apa_single_filter_performance.py",
    "04_prepare_pairwise_winrate.py",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run sim data performance workflow scripts.")
    parser.add_argument(
        "--data-root",
        default=None,
        help="Forwarded to prepare scripts that accept --data-root (optional).",
    )
    parser.add_argument(
        "--filter-type-1-suffix",
        default="0.05",
        help="Forwarded to 02_prepare_apa_residuals.py (default: 0.05).",
    )
    parser.add_argument(
        "--gt-count-total",
        type=float,
        default=5000.0,
        help="Forwarded to 02_prepare_apa_residuals.py (default: 5000).",
    )
    parser.add_argument(
        "--f1-lower-quantile",
        type=float,
        default=0.01,
        help="Forwarded to 02_prepare_apa_residuals.py (default: 0.01).",
    )
    parser.add_argument(
        "--f1-upper-quantile",
        type=float,
        default=0.99,
        help="Forwarded to 02_prepare_apa_residuals.py (default: 0.99).",
    )
    parser.add_argument(
        "--include-stage-files",
        action="store_true",
        help="Forwarded to 00_prepare_data.py (default: exclude stage files).",
    )
    parser.add_argument(
        "--max-files-per-domain",
        type=int,
        default=0,
        help="Forwarded to 00_prepare_data.py (debug cap, default: 0 for all).",
    )
    parser.add_argument(
        "--plot-config",
        default=str(DEFAULT_PLOT_CONFIG),
        help=(
            "JSON config for panel size/layout/style overrides. "
            "Forwarded via SIM_DATA_PERF_PLOT_CONFIG env var."
        ),
    )
    parser.add_argument(
        "--force-prepare",
        action="store_true",
        help=(
            "Force rerun prepare scripts even when intermediate outputs already exist. "
            "Default behavior reuses existing intermediate files."
        ),
    )
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help=(
            "Skip all prepare scripts and render panels only. "
            "Requires existing intermediate files."
        ),
    )
    parser.add_argument(
        "--with-pipeline-manifest",
        action="store_true",
        help=(
            "Also run 01_prepare_dapars2_pipeline_data.py (not required for plotting panels)."
        ),
    )
    return parser.parse_args()


def _output_paths(script_name: str) -> list[Path]:
    names = PREPARE_OUTPUTS.get(script_name, [])
    return [INTERMEDIATE_DIR / name for name in names]


def _all_outputs_exist(script_name: str) -> bool:
    paths = _output_paths(script_name)
    return bool(paths) and all(path.exists() for path in paths)


def _build_cmd(name: str, args: argparse.Namespace) -> list[str]:
    script_path = SCRIPT_DIR / name
    cmd = [sys.executable, str(script_path)]
    if name in set(PREPARE_SCRIPTS) and args.data_root:
        cmd.extend(["--data-root", args.data_root])
    if name == "00_prepare_data.py" and args.include_stage_files:
        cmd.append("--include-stage-files")
    if name == "00_prepare_data.py" and args.max_files_per_domain > 0:
        cmd.extend(["--max-files-per-domain", str(args.max_files_per_domain)])
    if name == "02_prepare_apa_residuals.py":
        cmd.extend(
            [
                "--filter-type-1-suffix",
                str(args.filter_type_1_suffix),
                "--gt-count-total",
                str(args.gt_count_total),
                "--f1-lower-quantile",
                str(args.f1_lower_quantile),
                "--f1-upper-quantile",
                str(args.f1_upper_quantile),
            ]
        )
    return cmd


def _run_scripts(names: Iterable[str], *, args: argparse.Namespace, env: dict[str, str]) -> None:
    for name in names:
        cmd = _build_cmd(name, args)
        print(f"[run_all] Running: {name}")
        subprocess.check_call(cmd, env=env)


def _missing_required_outputs(scripts: Iterable[str]) -> list[str]:
    missing: list[str] = []
    for name in scripts:
        for path in _output_paths(name):
            if not path.exists():
                missing.append(str(path))
    return missing


def main() -> None:
    args = parse_args()
    child_env = os.environ.copy()
    child_env["SIM_DATA_PERF_PLOT_CONFIG"] = str(Path(args.plot_config).expanduser().resolve())
    required_prepare = list(DEFAULT_REQUIRED_PREPARE)
    if args.with_pipeline_manifest:
        required_prepare.append("01_prepare_dapars2_pipeline_data.py")

    if args.plots_only:
        missing = _missing_required_outputs(required_prepare)
        if missing:
            preview = "\n".join(f"- {item}" for item in missing[:8])
            raise FileNotFoundError(
                "plots-only mode requested, but required intermediate files are missing:\n"
                f"{preview}"
            )
        print("[run_all] plots-only mode: skipping all prepare scripts.")
    else:
        to_run_prepare: list[str] = []
        reused_prepare: list[str] = []
        for name in required_prepare:
            if args.force_prepare or (not _all_outputs_exist(name)):
                to_run_prepare.append(name)
            else:
                reused_prepare.append(name)
        if reused_prepare:
            print(f"[run_all] Reusing intermediates, skip prepare: {', '.join(reused_prepare)}")
        if to_run_prepare:
            _run_scripts(to_run_prepare, args=args, env=child_env)
        else:
            print("[run_all] All required intermediates already exist; no prepare scripts run.")

    _run_scripts(PANEL_SCRIPTS, args=args, env=child_env)


if __name__ == "__main__":
    main()
