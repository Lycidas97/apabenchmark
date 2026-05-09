# 06_plotting

This directory is organized by plotting topics.

For release reproduction, prefer the top-level `scripts/run_full_repro.sh`
entrypoint, which runs the connected simulation, peak-model, and computation
resource plotting topics. Use topic-level commands for debugging or regenerating
individual figure panels.

## Design Rules

- One topic per folder.
- Legacy notebooks are kept only under `reference/`.
- Topic scripts must read project-level raw inputs and generate topic-level intermediate tables in `data/intermediate/`.
- Each panel has its own script and explicit output path under `figures/`.
- Plot style (font sizes, linewidth, palette, size presets, DPI) must be applied from `_shared/style.py`.

## Shared Layer

- `_shared/constants.py`: tool/protocol maps, ordering, palette.
- `_shared/paths.py`: project and topic path resolution.
- `_shared/io.py`: table and metadata read/write helpers.
- `_shared/style.py`: reference-notebook visual style and figure helpers.

## Topic Workflow

For each topic folder:

1. Run `scripts/00_prepare_data.py`.
2. Run panel scripts in numeric order (`10_*`, `20_*`, ...).
3. Outputs are written to `figures/`.

Example:

```bash
python scripts/06_plotting/pas_detect_performance/scripts/00_prepare_data.py
python scripts/06_plotting/pas_detect_performance/scripts/10_panel_overall.py
python scripts/06_plotting/pas_detect_performance/scripts/20_panel_metric_by_protocol.py
```

## Topic Glossary

- `sim_data_performance`: simulation benchmark performance across tools and
  protocols.
- `raw_data_performance`: optional raw-data benchmark and paired-validation
  plots.
- `pas_detect_performance`: PAS detection-oriented panels.
- `pas_quantification_performance`: PAS quantification-oriented panels.
- `apa_detect_performance`: differential APA detection-oriented panels.
- `peak_model_params`: peak model and coverage-shape parameter panels.
- `computation_resource_consumption`: runtime, memory, and CPU resource panels.

These topic directory names are stable workflow paths for the v0.1 pre-release;
runners, configs, and tests refer to them directly.
