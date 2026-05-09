# raw_data_performance

This topic tracks plotting data prep for benchmark raw-data runs from
`05_benchmark_analysis` (`pipeline_mode=raw` or `raw_parallel`).

## Objective

1. Record upstream raw-data organization and table schemas.
2. Keep a stable topic folder layout for future panel scripts.
3. Aggregate raw run summaries into topic-level intermediate files for plotting.
4. Provide recall/gt-focused filtered data via a config-driven rule set.

## Upstream Data Organization (from `05_benchmark_analysis`)

Raw run root pattern:

- `data/result/raw/{run_id}/...`
- default `run_id` in workflow is `default`.

Core outputs produced by Snakemake raw rules:

- `raw_pair_manifest/sample_top3_celltypes.tsv`
- `raw_pair_manifest/sample_pair_manifest.tsv`
- `raw_pair_manifest/prepare_metadata.tsv`
- `raw_qc/pas_filter_overview.tsv`
- `raw_performance/summary/raw_performance_match_metrics.tsv`
- `raw_performance/summary/raw_performance_pas_quantify_metrics.tsv`
- `raw_performance/summary/raw_performance_de_apa_metrics.tsv`
- `raw_performance/summary/raw_performance_missing_failed.tsv`
- `raw_dapars2_parallel/summary/raw_dapars2_parallel_metrics.tsv`
- `raw_dapars2_parallel/summary/raw_dapars2_parallel_missing_failed.tsv`
- `raw_dapars2_parallel/summary/raw_dapars2_parallel_{gt,nondapars,dapars}_failures.tsv`
- Alias (copied from parallel summary when available):
  - `raw_dapars2_format/summary/raw_dapars2_metrics.tsv`
  - `raw_dapars2_format/summary/raw_dapars2_missing_failed.tsv`

### Important schema notes

- Pair manifest (`sample_pair_manifest.tsv`):
  - `sample`, `pair_rank`, `pair_label`, `group1`, `group2`,
    `group1_clean`, `group2_clean`, `group1_count`, `group2_count`
- QC (`pas_filter_overview.tsv`):
  - `sample`, `pair_label`, `tool`, `selected_cells`, `pas_total`,
    `pas_kept_nonzero`, `pas_dropped_zero`, `bed_rows_kept`,
    `barcode_transform`
- Raw performance summaries:
  - match: `sample`, `tool`, `protocol`, `match_type`, `tp`, `fp`, `fn`,
    `precision`, `recall`, `f1`
  - quantify: `sample`, `tool`, `protocol`, `match_type`, `cor_pas`,
    `rmse_pas`, `mae_pas`, `mape_pas`, `rmse_pas_ct`, `mae_pas_ct`,
    `mape_pas_ct`
  - de_apa / dapars2 metrics: `sample`, `tool`, `protocol`, `match_type`,
    `filter_type_1`, `filter_type_2`, `precision`, `recall`, `f1`,
    `gt_count`, `pd_count`
- Failure summaries:
  - `sample`, `pair_label`, `tool`, `stage`, `message`

Note: in summary metric tables, `sample` is typically encoded as
`{sample}__{pair_label}`.

## Topic Directory Layout

```text
scripts/06_plotting/raw_data_performance/
  README.md
  config/
    dapars2_pipeline.json
    recall_gt_filter.json
    plot_panels.json
  figures/
    .gitkeep
  reference/
    .gitkeep
  scripts/
    00_prepare_data.py
    01_prepare_recall_gt_data.py
    10_panel_raw_recall_boxplot.py
    20_panel_raw_recall_vs_sim_f1_scatter.py
    30_panel_dapars2_pipeline_recall_boxplot.py
    run_all.py
```

Runtime outputs:

```text
scripts/06_plotting/raw_data_performance/data/intermediate/
  raw_data_performance_prepared.parquet
  raw_data_pair_manifest.parquet
  raw_data_qc_overview.parquet
  raw_data_failures.parquet
  metadata.json
  raw_data_recall_gt_filtered.parquet
  raw_data_recall_gt_distribution.tsv
  raw_data_recall_gt_metadata.json
  raw_recall_vs_sim_f1_scatter_points.tsv
  raw_dapars2_pipeline_recall_boxplot_points.tsv
  raw_data_dapars2_pipeline_manifest.parquet
  raw_data_dapars2_pipeline_summary.parquet
  raw_data_dapars2_pipeline_metadata.json
```

## Scripts

`00_prepare_data.py`:

- reads one raw run root and aggregates all summary metrics.
- outputs topic-level base tables and `metadata.json`.

`01_prepare_recall_gt_data.py`:

- reads `raw_data_performance_prepared.parquet`.
- keeps only configured `metric_domains` and configured
  `(filter_type_1, filter_type_2)` pairs.
- converts `recall`/`gt_count` to numeric and drops invalid rows.
- applies gt filters from config:
  - `drop_gt_zero`
  - `min_gt_count`
- writes:
  - `raw_data_recall_gt_filtered.parquet`
  - `raw_data_recall_gt_distribution.tsv`
  - `raw_data_recall_gt_metadata.json`

`02_prepare_dapars2_pipeline_data.py`:

- prepares a dedicated raw DaPars2/scMAPA pipeline manifest, aligned with
  sim-data separate-pipeline handling.
- default config focuses on `dapars2` and `scmapa` so their pipeline outputs can
  be tracked separately from other tools.
- outputs:
  - `raw_data_dapars2_pipeline_manifest.parquet`
  - `raw_data_dapars2_pipeline_summary.parquet`
  - `raw_data_dapars2_pipeline_metadata.json`

`10_panel_raw_recall_boxplot.py`:

- renders raw recall boxplot by tool.
- hue uses configured `filter_pairs` from `config/recall_gt_filter.json`.
- output: `figures/raw_recall_boxplot.pdf`.

`20_panel_raw_recall_vs_sim_f1_scatter.py`:

- renders raw recall vs sim F1 scatter by tool.
- x-axis: raw recall mean, y-axis: sim F1 mean.
- x/y error bars are group-level standard deviations.
- sim source priority:
  - `scripts/06_plotting/sim_data_performance/data/intermediate/sim_data_performance_prepared.parquet`
  - fallback: `data/result/performance/*/*_de_apa_performance.tsv`
- outputs:
  - `figures/raw_recall_vs_sim_f1_scatter_errorbar.pdf`
  - `data/intermediate/raw_recall_vs_sim_f1_scatter_points.tsv`

`30_panel_dapars2_pipeline_recall_boxplot.py`:

- renders a parallel recall boxplot for dedicated DaPars2/scMAPA-aligned
  pipeline data (`metric_domain=raw_dapars2_performance`).
- applies the same gt/recall filtering rules from `recall_gt_filter.json`
  for comparability with main raw plots.
- outputs:
  - `figures/raw_dapars2_pipeline_recall_boxplot.pdf`
  - `data/intermediate/raw_dapars2_pipeline_recall_boxplot_points.tsv`

## Config

`config/recall_gt_filter.json` controls recall/gt filtering and is intended to
be consumed by plotting scripts.

Current defaults:

- `metric_domains`:
  - `raw_de_apa_performance`
  - `raw_dapars2_performance`
- `filter_pairs`:
  - `dexseq_0.05 | dexseq_log2fc_1`
  - `dexseq_0.05 | dexseq_log2fc_1.25`
  - `fisher_0.05 | MPRO_0.5`
- `gt_rules`:
  - `drop_gt_zero: true`
  - `min_gt_count: 1`
- `recall_rules`:
  - `drop_all_tool_zero_recall_sample_pairs: true`

`config/plot_panels.json` controls panel plotting parameters.

Commonly changed keys:

- `boxplot`:
  - `input_name`, `output_name`, `figsize_preset`
  - `ymin`, `ymax`, axis labels, tick rotation
- `scatter`:
  - `output_figure`, `output_points`, `figsize_preset`
  - `xmin/xmax/ymin/ymax`, axis labels
  - `errorbar_stat` (`std` or `sem`)
  - `groupby` (default `['tool']`)
  - `sim_prepared_input` and fallback source config

`config/dapars2_pipeline.json` controls dedicated pipeline aggregation:

- `pipeline_root` (default: `raw_dapars2_parallel`)
- `tools_mode` (`selected` or `all`)
- `selected_tools` (default: `dapars2`, `scmapa`)
- domain toggles (`converted_inputs`, `differential_apa`, `performance`)

## CLI Example

```bash
python scripts/06_plotting/raw_data_performance/scripts/00_prepare_data.py \
  --data-root /path/to/apabenchmark_final \
  --raw-run-id default
```

```bash
python scripts/06_plotting/raw_data_performance/scripts/01_prepare_recall_gt_data.py
  --config scripts/06_plotting/raw_data_performance/config/recall_gt_filter.json
```

```bash
python scripts/06_plotting/raw_data_performance/scripts/02_prepare_dapars2_pipeline_data.py \
  --run-root /path/to/apabenchmark_final/data/result/raw/top5_split5_pair5 \
  --config scripts/06_plotting/raw_data_performance/config/dapars2_pipeline.json
```

Or:

```bash
python scripts/06_plotting/raw_data_performance/scripts/run_all.py \
  --data-root /path/to/apabenchmark_final \
  --raw-run-id default \
  --config scripts/06_plotting/raw_data_performance/config/recall_gt_filter.json \
  --plot-config scripts/06_plotting/raw_data_performance/config/plot_panels.json \
  --dapars2-config scripts/06_plotting/raw_data_performance/config/dapars2_pipeline.json
```
