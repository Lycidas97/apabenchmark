Generating PAS Annotations
==========================

This stage integrates mouse and human PAS annotation sources and writes the
filtered annotation files consumed by downstream stages.

Prerequisites
-------------

- Create and activate the `apasim` environment from `scripts/envs/apasim.yml`.
- Place the required raw annotation inputs under `data/raw_data/annotations/`
  as described in `docs/full_repro_data_flow.md`.
- Run from this directory, or use the top-level `scripts/run_full_repro.sh`.

Entrypoint
----------

```bash
cd scripts/01_generate_annotation
bash run_annotation.sh
```

Run order:

1. `mouse_pas_integration.py`
2. `human_pas_integration.py`
3. `filter_pas.py`
4. `summarize_annotation_stage_counts.py`

Outputs
-------

Generated files are written to `data/int_data/annotations/`, including
`mouse_integrated_pas.bed`, `human_integrated_pas.bed`, `*_pas_500.bed`,
`*_pas_1000.bed`, `*_pas_single.bed`, and
`supplementary_annotation_stage_counts.tsv`.

The source notebooks are retained as references; the scripts are the public
reproduction entrypoint.
