"""Path helpers for plotting scripts."""

from __future__ import annotations

from pathlib import Path


def project_root() -> Path:
    return Path(__file__).resolve().parents[3]


def plot_root() -> Path:
    return project_root() / "scripts" / "06_plotting"


def topic_root(topic: str) -> Path:
    return plot_root() / topic


def topic_data_dir(topic: str) -> Path:
    return topic_root(topic) / "data"


def topic_intermediate_dir(topic: str) -> Path:
    return topic_data_dir(topic) / "intermediate"


def topic_figure_dir(topic: str) -> Path:
    return topic_root(topic) / "figures"


def resolve_project_path(relative_path: str) -> Path:
    return project_root() / relative_path
