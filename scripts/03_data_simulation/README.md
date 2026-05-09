# 03 Data Simulation

This stage consumes integrated PAS annotations and peak summaries, then
generates simulated PAS sets, BAMs, expression matrices, and CPRSB resource
benchmark inputs.

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
cd scripts/03_data_simulation
APABENCHMARK_MM10_GENOME_FASTA=/path/to/GRCm38.p6.genome.fa \
conda run -n apasim snakemake --cores 8 --use-conda
```

## Inputs

- `data/int_data/annotations/mouse_integrated_pas.bed` from Stage 01.
- Peak analysis outputs from Stage 02 under `data/int_data`.
- Mouse genome FASTA supplied by `--mm10-fasta` or
  `APABENCHMARK_MM10_GENOME_FASTA`.

The full input contract is maintained in
`docs/full_repro_data_flow.md` and `docs/external_inputs_manifest.md`.

## Outputs

Outputs are written under `data/sim_data`, including simulated PAS files,
simulated BAMs, expression matrices, and CPRSB resource benchmark inputs.

## Notes

Stages 01-04 use `../../data` paths relative to their stage directories. For a
full upstream run, keep `data/` beside this repository checkout.
