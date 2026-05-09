#!/usr/bin/env python3
"""Render runtime panel for computation resource consumption."""

from __future__ import annotations

from panel_common import load_prepared_table, render_metric_panel


def main() -> None:
    df = load_prepared_table()
    render_metric_panel(
        df=df,
        metric="runtime_min",
        y_label="Run Time (min)",
        output_name="runtime.pdf",
    )


if __name__ == "__main__":
    main()
