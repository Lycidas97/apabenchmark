"""Shared plotting style, matching reference notebooks."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns

from .constants import PALETTE_HEX, TOOL_COLOR_HEX

MM = 1 / 25.4

SIZE_PRESETS_MM = {
    "small_square": (30, 30),
    "single_column": (75, 100),
    "narrow_tall": (60, 170),
    "mid_wide": (150, 30),
    "full_wide": (185, 100),
    "residual_scatter": (45, 70),
}


def apply_reference_style() -> str:
    """Apply global style settings used by reference notebooks."""
    style_used = "default"
    try:
        import scienceplots  # noqa: F401

        plt.style.use(["science", "nature"])
        style_used = "science+nature"
    except Exception:
        plt.style.use("default")

    # Keep text as editable fonts in vector outputs (PDF/SVG).
    plt.rcParams["text.usetex"] = False
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Nimbus Sans", "Nimbus Sans L"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["xtick.labelsize"] = 7
    plt.rcParams["ytick.labelsize"] = 7
    plt.rcParams["axes.labelsize"] = 8
    plt.rcParams["xtick.top"] = False
    plt.rcParams["ytick.right"] = False
    # Match reference notebooks: regular panels use inward major ticks.
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["xtick.bottom"] = True
    plt.rcParams["ytick.left"] = True
    plt.rcParams["xtick.minor.bottom"] = False
    plt.rcParams["lines.linewidth"] = 0.5
    plt.rcParams["legend.fontsize"] = 8
    plt.rcParams["hatch.linewidth"] = 0.5
    plt.rcParams["savefig.dpi"] = 300
    # Keep text as editable font objects in vector outputs.
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"
    return style_used


def palette(n_colors: int | None = None):
    n = n_colors if n_colors is not None else len(PALETTE_HEX)
    return sns.color_palette(PALETTE_HEX, n)


def tool_palette(
    tool_order: list[str],
    overrides: dict[str, str] | None = None,
    fallback: str = "#929292",
) -> dict[str, str]:
    colors = dict(TOOL_COLOR_HEX)
    if overrides:
        colors.update(overrides)
    return {tool: colors.get(tool, fallback) for tool in tool_order}


def figsize_from_mm(width_mm: float, height_mm: float) -> tuple[float, float]:
    return (width_mm * MM, height_mm * MM)


def figsize_from_preset(name: str) -> tuple[float, float]:
    width_mm, height_mm = SIZE_PRESETS_MM[name]
    return figsize_from_mm(width_mm, height_mm)


def save_figure(fig: plt.Figure, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", dpi=plt.rcParams.get("savefig.dpi", 300))


def draw_placeholder(ax, title: str, message: str) -> None:
    ax.axis("off")
    ax.text(0.02, 0.92, title, fontsize=9, ha="left", va="top", transform=ax.transAxes)
    ax.text(0.02, 0.75, message, fontsize=8, ha="left", va="top", transform=ax.transAxes)
