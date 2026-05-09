# APABenchmark

APABenchmark is a reproducible workflow for preparing APA benchmark inputs,
running benchmarked tools, summarizing performance, and generating figure-ready
outputs.

## Quick Start

From the repository root, bootstrap the checkout:

```bash
./scripts/bootstrap.sh
```

Add the prepared empirical BAM/BAI files:

```text
data/raw_data/bam_for_peak_extraction/{sample}.bam
data/raw_data/bam_for_peak_extraction/{sample}.bam.bai
```

Run a preflight check:

```bash
./scripts/run_full_repro.sh \
  --data-root "$PWD" \
  --mm10-fasta /path/to/GRCm38.p6.genome.fa \
  --cores 8 \
  --dry-run
```

Run the workflow:

```bash
./scripts/run_full_repro.sh \
  --data-root "$PWD" \
  --mm10-fasta /path/to/GRCm38.p6.genome.fa \
  --cores 8
```

`APABENCHMARK_MM10_GENOME_FASTA` can be used instead of `--mm10-fasta`.

## What Bootstrap Does

`scripts/bootstrap.sh` prepares the standard local layout:

- creates the workflow conda environments
- creates the expected `data/` directories
- downloads public annotation inputs
- downloads benchmark Singularity/Apptainer images
- builds tool reference files from public annotations

Useful options:

```bash
./scripts/bootstrap.sh --skip-sif
./scripts/bootstrap.sh --skip-tool-references
./scripts/bootstrap.sh --skip-scape
./scripts/bootstrap.sh --force
```

If you only need the conda environments, run:

```bash
./scripts/setup_envs.sh
```

## Inputs

Bootstrap prepares the public/reference side of the workflow. The empirical
BAM/BAI inputs are large prepared files and should be placed locally under
`data/raw_data/bam_for_peak_extraction/`.

Source URLs, versions, and checksums are listed in
[docs/external_inputs_manifest.md](docs/external_inputs_manifest.md). Public
dataset-level accessions for empirical inputs are summarized in
[docs/raw_dataset_inventory.md](docs/raw_dataset_inventory.md).

## Workflow

```text
public annotations + empirical BAMs + tool references
  -> 01_generate_annotation
  -> 02_peak_analysis
  -> 03_data_simulation
  -> 04_benchmark
  -> 05_benchmark_analysis
  -> 06_plotting
```

Stage-level commands, inputs, and outputs are documented in
[docs/full_repro_data_flow.md](docs/full_repro_data_flow.md) and in the
per-stage README files under `scripts/`.

## Outputs

Generated files are written under `data/`:

- `data/int_data/`: integrated annotations, processed BAM-derived intermediates,
  peak tables, and peak-shape summaries
- `data/sim_data/`: simulated PAS sets, simulated BAMs, expression matrices, and
  CPRSB inputs
- `data/result/benchmark/`: benchmarked tool outputs and resource profiles
- `data/result/performance/`: performance summaries
- `data/result/differential_apa/`: differential APA outputs
- `data/result/`: plotting-source tables and downstream analysis products

Generated figures, logs, pid files, and large intermediate tables are ignored by
git.

## Layout

```text
apabenchmark/
  docs/
  scripts/
    bootstrap.sh
    setup_envs.sh
    run_full_repro.sh
    01_generate_annotation/
    02_peak_analysis/
    03_data_simulation/
    04_benchmark/
    05_benchmark_analysis/
    06_plotting/
  data/                  # external or generated; ignored by git
```

See [docs/directory_structure.md](docs/directory_structure.md) for the full
directory contract.

## Citation

The manuscript citation is not yet available. Until the manuscript is public,
please cite the repository metadata in [CITATION.cff](CITATION.cff).

## License

This repository is distributed under the MIT License. See [LICENSE](LICENSE).
