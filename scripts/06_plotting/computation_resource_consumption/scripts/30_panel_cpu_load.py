#!/usr/bin/env python3
"""Render mean CPU load panel for computation resource consumption."""

from __future__ import annotations

from panel_common import load_prepared_table, render_metric_panel


def main() -> None:
    df = load_prepared_table()
    render_metric_panel(
        df=df,
        metric="mean_load",
        y_label="Mean cpu load (%)",
        output_name="cpu_load.pdf",
    )


if __name__ == "__main__":
    main()
