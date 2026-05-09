#!/usr/bin/env python3
"""Run raw-data performance plotting workflow."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run raw-data performance workflow scripts.")
    parser.add_argument(
        "--data-root",
        default=None,
        help="Forwarded to prepare scripts that accept --data-root (optional).",
    )
    parser.add_argument(
        "--raw-run-id",
        default="default",
        help="Forwarded to 00_prepare_data.py (default: default).",
    )
    parser.add_argument(
        "--run-root",
        default="",
        help="Forwarded to 00_prepare_data.py; explicit raw run root override.",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Forwarded to 01_prepare_recall_gt_data.py (optional).",
    )
    parser.add_argument(
        "--plot-config",
        default=None,
        help="Forwarded to panel scripts (optional).",
    )
    parser.add_argument(
        "--dapars2-config",
        default=None,
        help="Forwarded to 02_prepare_dapars2_pipeline_data.py (optional).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    prep_cmd = [sys.executable, str(SCRIPT_DIR / "00_prepare_data.py")]
    if args.data_root:
        prep_cmd.extend(["--data-root", args.data_root])
    if args.raw_run_id:
        prep_cmd.extend(["--raw-run-id", str(args.raw_run_id)])
    if args.run_root:
        prep_cmd.extend(["--run-root", args.run_root])
    subprocess.check_call(prep_cmd)

    subset_cmd = [sys.executable, str(SCRIPT_DIR / "01_prepare_recall_gt_data.py")]
    if args.config:
        subset_cmd.extend(["--config", args.config])
    subprocess.check_call(subset_cmd)

    dapars2_cmd = [sys.executable, str(SCRIPT_DIR / "02_prepare_dapars2_pipeline_data.py")]
    if args.data_root:
        dapars2_cmd.extend(["--data-root", args.data_root])
    if args.raw_run_id:
        dapars2_cmd.extend(["--raw-run-id", str(args.raw_run_id)])
    if args.run_root:
        dapars2_cmd.extend(["--run-root", args.run_root])
    if args.dapars2_config:
        dapars2_cmd.extend(["--config", args.dapars2_config])
    subprocess.check_call(dapars2_cmd)

    box_cmd = [sys.executable, str(SCRIPT_DIR / "10_panel_raw_recall_boxplot.py")]
    if args.config:
        box_cmd.extend(["--filter-config", args.config])
    if args.plot_config:
        box_cmd.extend(["--plot-config", args.plot_config])
    subprocess.check_call(box_cmd)

    scatter_cmd = [sys.executable, str(SCRIPT_DIR / "20_panel_raw_recall_vs_sim_f1_scatter.py")]
    if args.data_root:
        scatter_cmd.extend(["--data-root", args.data_root])
    if args.config:
        scatter_cmd.extend(["--filter-config", args.config])
    if args.plot_config:
        scatter_cmd.extend(["--plot-config", args.plot_config])
    subprocess.check_call(scatter_cmd)

    recall_pd_scatter_cmd = [sys.executable, str(SCRIPT_DIR / "21_panel_raw_recall_vs_pd_scatter.py")]
    if args.config:
        recall_pd_scatter_cmd.extend(["--filter-config", args.config])
    if args.plot_config:
        recall_pd_scatter_cmd.extend(["--plot-config", args.plot_config])
    subprocess.check_call(recall_pd_scatter_cmd)

    dapars2_panel_cmd = [sys.executable, str(SCRIPT_DIR / "30_panel_dapars2_pipeline_recall_boxplot.py")]
    if args.config:
        dapars2_panel_cmd.extend(["--filter-config", args.config])
    if args.dapars2_config:
        dapars2_panel_cmd.extend(["--dapars2-config", args.dapars2_config])
    if args.plot_config:
        dapars2_panel_cmd.extend(["--plot-config", args.plot_config])
    subprocess.check_call(dapars2_panel_cmd)

    dapars2_scatter_cmd = [sys.executable, str(SCRIPT_DIR / "31_panel_dapars2_pipeline_recall_vs_sim_f1_scatter.py")]
    if args.data_root:
        dapars2_scatter_cmd.extend(["--data-root", args.data_root])
    if args.plot_config:
        dapars2_scatter_cmd.extend(["--plot-config", args.plot_config])
    subprocess.check_call(dapars2_scatter_cmd)

    heatmap_main_cmd = [sys.executable, str(SCRIPT_DIR / "40_panel_recall_pairwise_winrate_heatmap.py")]
    if args.plot_config:
        heatmap_main_cmd.extend(["--plot-config", args.plot_config])
    heatmap_main_cmd.extend(["--panel-key", "heatmap_main_recall"])
    subprocess.check_call(heatmap_main_cmd)

    heatmap_dapars2_cmd = [sys.executable, str(SCRIPT_DIR / "40_panel_recall_pairwise_winrate_heatmap.py")]
    if args.plot_config:
        heatmap_dapars2_cmd.extend(["--plot-config", args.plot_config])
    heatmap_dapars2_cmd.extend(["--panel-key", "heatmap_dapars2_recall"])
    subprocess.check_call(heatmap_dapars2_cmd)

    dapars2_combo_cmd = [sys.executable, str(SCRIPT_DIR / "51_panel_dapars2_recall_scatter_heatmap_combined.py")]
    if args.plot_config:
        dapars2_combo_cmd.extend(["--plot-config", args.plot_config])
    subprocess.check_call(dapars2_combo_cmd)

    main_combo_cmd = [sys.executable, str(SCRIPT_DIR / "50_panel_main_recall_scatter_heatmap_combined.py")]
    if args.plot_config:
        main_combo_cmd.extend(["--plot-config", args.plot_config])
    subprocess.check_call(main_combo_cmd)


if __name__ == "__main__":
    main()
