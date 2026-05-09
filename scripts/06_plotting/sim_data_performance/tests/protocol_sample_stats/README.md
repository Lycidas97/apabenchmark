# protocol_sample_stats

Small test utility to summarize which protocols appear in sim-data plotting
samples and how many unique samples belong to each protocol.

This is a developer utility for auditing plotting inputs, not a release
validation requirement for the connected workflow.

## What it reports

- protocol list used by sim-data plotting
- unique sample count per protocol
- total unique sample count for single-cell protocols
- total unique sample count for spatial transcriptomics protocols

## Run

```bash
python scripts/06_plotting/sim_data_performance/tests/protocol_sample_stats/00_count_protocol_samples.py
```

Optional input override:

```bash
python scripts/06_plotting/sim_data_performance/tests/protocol_sample_stats/00_count_protocol_samples.py \
  --input-parquet scripts/06_plotting/sim_data_performance/data/intermediate/sim_data_performance_prepared.parquet
```

## Outputs

- `output/protocol_sample_counts.parquet`
- `output/sample_protocol_map.parquet`
- `output/summary.json`
- `output/summary.md`
