#!/usr/bin/env python3
"""Run full peak model parameter plotting workflow."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

SCRIPT_DIR = Path(__file__).resolve().parent


SCRIPTS = [
    "00_prepare_data.py",
    "05_prepare_protocol_normalized_peak_profiles.py",
    "10_panel_aic_bic.py",
    "20_panel_model_params.py",
    "25_panel_shape_core_params.py",
    "40_panel_protocol_normalized_peak_profiles.py",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run full peak model parameter plotting workflow.")
    parser.add_argument(
        "--data-root",
        default=None,
        help="Project root containing data/ or the data directory itself.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for name in SCRIPTS:
        script_path = SCRIPT_DIR / name
        cmd = [sys.executable, str(script_path)]
        if args.data_root and name in {
            "00_prepare_data.py",
            "05_prepare_protocol_normalized_peak_profiles.py",
        }:
            cmd.extend(["--data-root", args.data_root])
        subprocess.check_call(cmd)


if __name__ == "__main__":
    main()
