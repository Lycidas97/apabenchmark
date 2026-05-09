"""Runtime plot config loader for sim_data_performance panels."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
TOPIC_ROOT = SCRIPT_DIR.parent
DEFAULT_CONFIG_PATH = TOPIC_ROOT / "config" / "plot_params.json"
PLOT_CONFIG_ENV = "SIM_DATA_PERF_PLOT_CONFIG"

MM = 1 / 25.4

_DEFAULT_SIZE_PRESETS_MM = {
    "small_square": (30.0, 30.0),
    "single_column": (75.0, 100.0),
    "narrow_tall": (60.0, 170.0),
    "mid_wide": (150.0, 30.0),
    "full_wide": (185.0, 100.0),
    "residual_scatter": (45.0, 70.0),
    "residual_scatter_main": (45.0, 70.0),
    "f1_stringency_regression": (45.0, 70.0),
}

_CACHE: dict[str, Any] | None = None


def _nested_get(data: dict[str, Any], key_path: str, default: Any) -> Any:
    current: Any = data
    for part in key_path.split("."):
        if not isinstance(current, dict) or part not in current:
            return default
        current = current[part]
    return current


def config_path() -> Path:
    raw = os.environ.get(PLOT_CONFIG_ENV, "").strip()
    return Path(raw).expanduser().resolve() if raw else DEFAULT_CONFIG_PATH


def load_config() -> dict[str, Any]:
    global _CACHE
    if _CACHE is not None:
        return _CACHE

    path = config_path()
    if not path.exists():
        _CACHE = {}
        return _CACHE

    loaded = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(loaded, dict):
        raise ValueError(f"Plot config must be a JSON object: {path}")
    _CACHE = loaded
    return _CACHE


def cfg(key_path: str, default: Any) -> Any:
    return _nested_get(load_config(), key_path, default)


def cfg_float(key_path: str, default: float) -> float:
    value = cfg(key_path, default)
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def cfg_int(key_path: str, default: int) -> int:
    value = cfg(key_path, default)
    try:
        return int(value)
    except (TypeError, ValueError):
        return int(default)


def cfg_bool(key_path: str, default: bool) -> bool:
    value = cfg(key_path, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"1", "true", "yes", "y", "on"}:
            return True
        if lowered in {"0", "false", "no", "n", "off"}:
            return False
    return bool(default)


def figsize_from_preset(name: str) -> tuple[float, float]:
    size_map = dict(_DEFAULT_SIZE_PRESETS_MM)
    override_map = cfg("size_presets_mm", {})
    if isinstance(override_map, dict):
        for key, value in override_map.items():
            if isinstance(value, (list, tuple)) and len(value) == 2:
                try:
                    size_map[str(key)] = (float(value[0]), float(value[1]))
                except (TypeError, ValueError):
                    continue

    width_mm, height_mm = size_map[name]
    return (width_mm * MM, height_mm * MM)


def figsize_for_figure(output_name: str, fallback_preset: str) -> tuple[float, float]:
    """Resolve per-figure size in mm with preset fallback.

    Priority:
    1) figure_sizes_mm.<output_name> = [width_mm, height_mm]
    2) size_presets_mm.<fallback_preset>
    """
    figure_map = cfg("figure_sizes_mm", {})
    if isinstance(figure_map, dict):
        value = figure_map.get(output_name)
        if isinstance(value, (list, tuple)) and len(value) == 2:
            try:
                return figsize_from_mm(float(value[0]), float(value[1]))
            except (TypeError, ValueError):
                pass
    return figsize_from_preset(fallback_preset)


def figsize_from_mm(width_mm: float, height_mm: float) -> tuple[float, float]:
    return (float(width_mm) * MM, float(height_mm) * MM)


def apply_runtime_rcparams() -> None:
    plt.rcParams["savefig.dpi"] = cfg_int("render.savefig_dpi", 300)
    rc_overrides = cfg("rcParams", {})
    if isinstance(rc_overrides, dict):
        for key, value in rc_overrides.items():
            try:
                plt.rcParams[str(key)] = value
            except Exception:
                continue
