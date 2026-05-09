# Filter criteria cross-validation sensitivity analysis

This is an optional developer sensitivity analysis for DE-APA filter-criteria
cross-validation. It lives under the sim-data-performance plotting scripts
directory, but is not part of the connected `run_all.py` plotting workflow.

## Objective

Run a post-hoc sensitivity analysis using existing DE-APA performance matrices.
No tools are rerun. The analysis separates criterion selection from tool
evaluation to test whether the selected criteria and held-out tool rankings are
stable.

## Existing data

Primary input:

- `scripts/06_plotting/sim_data_performance/data/intermediate/sim_data_apa_residuals.parquet`

Required columns are already present:

- `sample`
- `tool`
- `protocol`
- `match_type`
- `filter_type_1`
- `filter_type_2`
- `filter_type_comb`
- `precision`
- `recall`
- `f1`
- `gt_recall`

Useful summary input:

- `scripts/06_plotting/sim_data_performance/data/intermediate/sim_data_apa_residuals_by_filter.parquet`

The full residual table has about 22M rows, so the first script writes a
row-level cache with only the required columns plus derived grouping metadata.
It also precomputes the split-invariant residuals once within each
`sample x tool x match_type` group. It does not average across simulated
samples or criteria before fitting residuals.

## Grouping unit

Preferred split unit:

- empirical sample background

Current table does not have an explicit `empirical_background` column. The
script derives it from `sample` by removing the simulation suffix:

```text
<background>_<genome>_pas<id>_gn<count>_rep<id>
```

Example:

```text
Chromium_human_embryolimb_5386STDY7557335_hg38_pas1_gn5000_rep1
-> Chromium_human_embryolimb_5386STDY7557335
```

Fallback:

- if the parser cannot match the suffix, it uses the full `sample` as the split
  unit, which is still valid and avoids row-level leakage.

## Analysis design

### Repeated grouped half-split validation

For each split:

1. Split grouped units 50:50 into train/test.
2. On train only:
   - recompute criterion stringency as mean `gt_recall`;
   - reuse the precomputed `F1 ~ stringency` residuals from each
     `sample x tool x match_type` group;
   - compute mean F1 residual per criterion on the train subset;
   - divide criteria into strict/moderate/lenient thirds by stringency;
   - choose the highest-residual criterion in each tier.
3. On held-out test:
   - apply the train-selected criterion;
   - compute tool precision/recall/F1 and F1 rank;
   - compute criterion-level train-held-out residual concordance within each
     training-defined tier;
   - also write same-tier raw-F1 and residual oracle gaps as auxiliary table
     columns;
   - compute rank concordance against full-data ranking.

Default:

- 100 random 50:50 grouped splits
- seed 1

### Leave-protocol-out validation

For each protocol:

1. Train on the other protocols.
2. Select criteria on train.
3. Evaluate selected criteria on the held-out protocol.
4. Compute the same ranking and held-out gap outputs.

## Outputs

Generated under:

- `scripts/06_plotting/sim_data_performance/scripts/filter_criteria_cv_sensitivity/output/`

Tables:

- `cv_selected_criteria.tsv`
- `cv_tool_metrics.tsv`
- `cv_regret.tsv`
- `cv_criterion_concordance.tsv`
- `cv_rank_concordance.tsv`
- `cv_top_frequency.tsv`
- `cv_metadata.json`

Figures:

- `figures/filter_criteria_cv_sensitivity.pdf`

## Commands

From repo root:

```bash
python3 scripts/06_plotting/sim_data_performance/scripts/filter_criteria_cv_sensitivity/00_prepare_cv_tables.py
python3 scripts/06_plotting/sim_data_performance/scripts/filter_criteria_cv_sensitivity/10_plot_cv_summary.py
```

For a faster dry run:

```bash
python3 scripts/06_plotting/sim_data_performance/scripts/filter_criteria_cv_sensitivity/00_prepare_cv_tables.py --n-random-splits 10
```

## Supplement figure summary

The plotting script writes a sensitivity summary with these panels:

Panel A:

- workflow schematic for grouped split and leave-protocol-out validation.

Panel B:

- criteria family selection frequency by tier.

Panel C:

- criterion residual concordance:
  Spearman correlation between train-set and held-out-set mean residuals across
  all criteria in the same train-defined tier.

Panel D:

- held-out tool rank frequency heatmap.

Panel E:

- rank concordance with full-data ranking, using Spearman rho and Kendall tau.

## Analysis conventions

1. `sample` parsing is used as the empirical-background grouping rule.
2. Criterion tiers follow the existing figure logic: thirds of the `gt_recall`
   range, where low `gt_recall` is strict and high `gt_recall` is lenient.
3. Criterion concordance uses train-set tier boundaries for candidate criteria.
4. Full-data reference ranking uses the current main-text criteria selected by
   the existing plotting helper.
5. Ranking stability uses F1 as the primary metric; precision and recall remain
   available in the output table for secondary checks.
