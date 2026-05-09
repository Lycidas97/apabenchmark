#!/usr/bin/env python3
"""Render max memory panel for computation resource consumption."""

from __future__ import annotations

from panel_common import load_prepared_table, render_metric_panel


def main() -> None:
    df = load_prepared_table()
    render_metric_panel(
        df=df,
        metric="max_mem",
        y_label="Max memory (MB)",
        output_name="max_mem.pdf",
    )


if __name__ == "__main__":
    main()
