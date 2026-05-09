#!/usr/bin/env python3
"""Run full computation resource consumption plotting workflow."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys

SCRIPT_DIR = Path(__file__).resolve().parent

SCRIPTS = [
    "00_prepare_data.py",
    "10_panel_runtime.py",
    "20_panel_max_mem.py",
    "30_panel_cpu_load.py",
]


def main() -> None:
    for name in SCRIPTS:
        script_path = SCRIPT_DIR / name
        subprocess.check_call([sys.executable, str(script_path)])


if __name__ == "__main__":
    main()
