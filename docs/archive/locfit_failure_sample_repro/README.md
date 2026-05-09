# locfit numeric failure repro (single sample)

This archived test folder reproduces a historical stage20 failure on one
sample. It is retained for debugging provenance and is not part of the active
public reproduction entrypoints.

- tool: `sierra`
- sample: `Microwell_human_embryo_SRR9887755_mm10_pas1_gn5000_rep2`
- known error: `newsplit: out of vertex space`

## Usage

From repo root:

```bash
bash docs/archive/locfit_failure_sample_repro/run_stage20_repro.sh
```

To run the full wrapper (stage10 -> stage20 -> stage30) for the same sample:

```bash
bash docs/archive/locfit_failure_sample_repro/run_wrapper_repro.sh
```

The script runs `20_stage_apa_and_te.R` against the existing stage10 bundle and writes logs to:

`docs/archive/locfit_failure_sample_repro/logs/`

## Notes

- This test reads existing files under `data/result/performance/` and does not copy raw data.
- Override the bundle path if needed:

```bash
STAGE10_BUNDLE=/path/to/custom_stage10_bundle.rds \
bash docs/archive/locfit_failure_sample_repro/run_stage20_repro.sh
```
