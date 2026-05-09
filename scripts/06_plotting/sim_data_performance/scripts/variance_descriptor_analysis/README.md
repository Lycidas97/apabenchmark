# Variance descriptor analysis

Prepare DE-APA performance tables joined to sample-level peak-shape descriptors
from `peak_model_params`, then analyze workflow-specific associations with
protocol and peak morphology after removing benchmark-instance baseline
difficulty.

This is an optional supplementary analysis utility. It is not part of the
connected `run_all.py` plotting workflow.

## Commands

From the repository root:

```bash
python3 scripts/06_plotting/sim_data_performance/scripts/variance_descriptor_analysis/00_prepare_descriptor_performance_table.py
python3 scripts/06_plotting/sim_data_performance/scripts/variance_descriptor_analysis/10_fit_peak_morphology_models.py
python3 scripts/06_plotting/sim_data_performance/scripts/variance_descriptor_analysis/20_make_supplementary_tables.py
```

If the sim-data intermediate table lives in another worktree, pass it explicitly:

```bash
python3 scripts/06_plotting/sim_data_performance/scripts/variance_descriptor_analysis/00_prepare_descriptor_performance_table.py \
  --apa-residuals /path/to/sim_data_apa_residuals.parquet \
  --single-filter-performance /path/to/sim_data_differential_apa_single_filter_performance.parquet
python3 scripts/06_plotting/sim_data_performance/scripts/variance_descriptor_analysis/10_fit_peak_morphology_models.py
python3 scripts/06_plotting/sim_data_performance/scripts/variance_descriptor_analysis/20_make_supplementary_tables.py
```

## Inputs

- `sim_data_performance/data/intermediate/sim_data_apa_residuals.parquet`
- `sim_data_performance/data/intermediate/sim_data_differential_apa_single_filter_performance.parquet`
- `peak_model_params/data/intermediate/peak_model_params_prepared.tsv`

## Output

Generated files are written under `output/`, which is intentionally ignored:

- `de_apa_descriptor_performance.parquet`
- `de_apa_descriptor_performance_metadata.json`
- `descriptor_summary.tsv`
- `unmatched_sim_samples.tsv`
- `de_apa_single_filter_descriptor_performance.parquet`
- `de_apa_single_filter_descriptor_performance_metadata.json`
- `single_filter_descriptor_summary.tsv`
- `single_filter_unmatched_sim_samples.tsv`
- `protocol_response_by_workflow.tsv`
- `protocol_response_estimates.tsv`
- `descriptor_response_by_workflow.tsv`
- `descriptor_explains_protocol_response.tsv`
- `sample_collapsed_protocol_response_by_workflow.tsv`
- `sample_collapsed_descriptor_response_by_workflow.tsv`
- `complete_case_protocol_response_by_workflow.tsv`
- `complete_case_descriptor_response_by_workflow.tsv`
- `study_loo_protocol_response_by_workflow.tsv`
- `study_loo_descriptor_response_by_workflow.tsv`
- `study_loo_protocol_sensitivity_summary.tsv`
- `study_loo_descriptor_sensitivity_summary.tsv`
- `morphology_coverage_diagnostics.tsv`
- `transform_consistency_summary.tsv`
- `pooled_protocol_descriptor_model_comparison.tsv`
- `morphology_descriptor_diagnostics.tsv`
- `morphology_instance_summary.tsv`
- `morphology_model_metadata.json`
- `supplementary_table_S8A_protocol_associated_variation.tsv`
- `supplementary_table_S8B_peak_morphology_associations.tsv`
- `supplementary_table_S8C_sensitivity_diagnostics.tsv`
- `Supplementary_Table_S8_protocol_peak_morphology_performance_variation.xlsx`

Both performance tables keep `precision`, `recall`, and `f1` as separate columns.
The three main-text DE-APA criteria are selected with the same helper used by
the plotting scripts, then averaged within each `tool x sample x match_type`.

The single-filter table follows the existing DaPars2/scMAPA-compatible plotting
pipeline: DaPars2/scMAPA use main `*_de_apa_performance.tsv` rows, other tools
use `*_de_apa_dapars_parallel_performance.tsv` sidecar rows, and all tools are
restricted to `scmapa_adPval_0.05 | scmapa_or_0.1`.
For model instance matching, `Sample_1_vs_Sample_2` is harmonized to `T1_vs_T2`
in `match_type_harmonized`; the original `match_type` is retained.

The peak descriptors are the three parameters plotted in
`peak_model_params/scripts/25_panel_shape_core_params.py`.

## Models

The model script fits each regime, metric, and response transform separately.
Clipped logit is the primary response transform; identity-scale rows are
sensitivity checks.

Regimes:

- `main_text_three_criteria`: selected main-text DE-APA criteria, averaged
  within `tool x sample x match_type`.
- `dapars2_scmapa_compatible_single_filter`: harmonized single-filter endpoint
  including DaPars2 and scMAPA.

Instance fixed effects:

- main criteria: `sample + match_type`
- single-filter: `sample + match_type_harmonized + filter_type_1 + filter_type_2`

Core response:

```text
R_i,m = g(Y_i,m) - mean_{m' != m} g(Y_i,m')
```

`R_i,m` is the transformed performance of workflow `m` relative to the other
workflows in the same benchmark instance. All primary workflow-specific models
use this response and do not use `scAPAtrap` as a reference.

Primary workflow-specific outputs:

- `protocol_response_by_workflow.tsv`: tests `R_i,m ~ protocol` separately for
  each workflow.
- `protocol_response_estimates.tsv`: protocol-level mean relative response for
  each workflow.
- `descriptor_response_by_workflow.tsv`: tests `R_i,m ~ Z_k` separately for
  each workflow and each descriptor.
- `descriptor_explains_protocol_response.tsv`: compares protocol-only,
  descriptor-only, descriptor-plus-protocol, and protocol-plus-within-protocol
  descriptor residual models for each workflow.

Peak descriptors:

- `upstream_asymmetry`
- `effective_width`
- `mean_distance_to_pas`

Secondary pooled output:

`pooled_protocol_descriptor_model_comparison.tsv` keeps the previous pooled
workflow interaction model comparisons as secondary global heterogeneity tests.
It should not be used for workflow-specific interpretation.

FDR:

- protocol-associated variation: BH within `regime x metric x response_transform`
  across workflows.
- peak-morphology associations: BH within `regime x metric x response_transform`
  across all workflow-by-descriptor tests.
- descriptor explanation: separate BH corrections for protocol,
  all-descriptor, protocol-after-descriptor, and within-protocol descriptor
  tests.

## Sensitivity and diagnostics

The model script also writes lightweight guardrail outputs for inference:

- `morphology_coverage_diagnostics.tsv` reports empirical sample counts,
  samples per protocol, valid workflow counts per instance (`K_i`),
  complete-case instance counts, workflow-missing instance counts, and duplicate
  `instance_id x workflow` checks.
- `sample_collapsed_*_by_workflow.tsv` averages `R_i,m` within each
  `empirical_sample_id x workflow` before refitting protocol and descriptor
  association models. This reduces sample-level pseudoreplication from multiple
  benchmark instances derived from the same empirical sample.
- `complete_case_*_by_workflow.tsv` keeps only benchmark instances containing
  every expected workflow in the regime, then recomputes the leave-one-out
  baseline and refits the workflow-specific models. This checks sensitivity to
  changes in the comparator set used in `R_i,m`.
- `study_loo_*` tables use study groups derived from the available source grouping field and refit the
  sample-collapsed models after leaving out one study at a time. This checks
  whether the main associations are driven by one study group.
- `transform_consistency_summary.tsv` compares `logit_clip` and identity-scale
  results by descriptor slope direction, slope rank, absolute slope rank,
  descriptor `R^2` rank, and protocol-associated `R^2` rank.

The sensitivity tables are not intended to replace the primary instance-level
`logit_clip` outputs. They are used to check that the main direction and
effect-size patterns are not artifacts of row-level replication, incomplete
workflow coverage, or response transformation.

## Compact supplementary tables

`20_make_supplementary_tables.py` builds paper-facing compact tables from the
full TSV outputs. These tables use primary `logit_clip` rows only. The script
also writes one Excel workbook for the paper supplement:

`Supplementary_Table_S8_protocol_peak_morphology_performance_variation.xlsx`

Recommended caption:

```text
Supplementary Table S8. Protocol- and peak-morphology-associated variation in
workflow performance.
```

The workbook contains a notes sheet plus three data sheets. Sheet C uses an
abbreviated name because Excel sheet names are limited to 31 characters.

- `supplementary_table_S8A_protocol_associated_variation.tsv` combines
  `protocol_response_by_workflow.tsv` and `protocol_response_estimates.tsv`.
  It reports protocol-associated variation in workflow performance plus the
  protocols with the highest and lowest mean relative performance for each
  workflow. The primary sample-level P value and FDR columns come from the
  sample-collapsed protocol model; instance-level statistics are retained as
  descriptive sensitivity results.
- `supplementary_table_S8B_peak_morphology_associations.tsv` is a compact
  version of `descriptor_response_by_workflow.tsv`, reporting each
  workflow-specific standalone peak-morphology association, descriptor `R^2`,
  and FDR. The primary sample-level slope, P value, and FDR columns come from
  the sample-collapsed descriptor model; instance-level slope, P value, FDR,
  and descriptor `R^2` are retained as explicitly named columns.
- `supplementary_table_S8C_sensitivity_diagnostics.tsv` combines
  descriptor/protocol variance decomposition with sample-collapsed,
  complete-case, study leave-one-out, transform-consistency, and coverage
  diagnostics. This table is intended as a compact sensitivity summary; the
  full component TSV files remain the source of record for detailed checks.

Important interpretation notes:

- Primary P value and FDR columns in sheets A and B come from
  sample-collapsed models fitted after averaging relative performance within
  `empirical_sample_id x workflow`.
- Slopes were estimated on the clipped-logit instance-adjusted
  relative-performance scale. Positive slopes indicate that a workflow
  performed better relative to other workflows as the corresponding descriptor
  increased.
- Study groups are derived from the available source grouping field for leave-one-study-out sensitivity.
- `R^2` values were calculated on the clipped-logit instance-adjusted
  relative-performance scale and should not be interpreted as variance explained
  in raw F1, precision, or recall.
- Descriptor rows are individual or standalone descriptor associations. They
  should not be described as independent contributions unless explicitly
  referring to conditional analyses.

Workbook sheets:

- `S8-notes`
- `S8A Protocol variation`
- `S8B Morphology assoc`
- `S8C Robustness diag`
