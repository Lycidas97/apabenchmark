# 02 Peak Analysis

This stage extracts non-overlapping peaks from raw BAM inputs and summarizes
peak model and coverage-shape statistics consumed by simulation and plotting.

## Recommended Entrypoint

For release reproduction, use the top-level runner:

```bash
./scripts/run_full_repro.sh \
  --data-root /path/to/apabenchmark_final \
  --mm10-fasta /path/to/GRCm38.p6.genome.fa \
  --cores 8
```

Use this directory directly only for stage-level debugging:

```bash
cd scripts/02_peak_analysis
conda run -n apasim snakemake --cores 8 --use-conda
```

## Inputs

- `data/int_data/annotations/mouse_pas_500.bed` from Stage 01.
- `data/raw_data/bam_for_peak_extraction/{sample}.bam`
- `data/raw_data/bam_for_peak_extraction/{sample}.bam.bai`

The full input contract is maintained in
`docs/full_repro_data_flow.md` and `docs/external_inputs_manifest.md`.
Dataset-level public sources for the prepared empirical BAM inputs are listed in
`docs/raw_dataset_inventory.md`.

## Outputs

Outputs are written under `data/int_data`, including:

- Peak tables and processed BAM-derived intermediates.
- `*_model_results.json`
- `*_peak_shape_by_pas.feather`
- `*_peak_shape_summary.feather`

## Notes

Snakemake creates rule-specific conda environments when `--use-conda` is used.
Do not commit generated logs, peak tables, or runtime outputs.
