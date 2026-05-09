# Full Reproducibility Data Flow

This document defines the connected APA benchmark data flow from prepared raw
inputs to generated analysis tables and figures.

## Objective

```text
raw annotations + raw BAMs + references
  -> 01_generate_annotation
  -> 02_peak_analysis
  -> 03_data_simulation
  -> 04_benchmark
  -> 05_benchmark_analysis
  -> 06_plotting
```

`data/` is the project data plane. It is ignored by git and must be supplied or
generated beside this repository.

Directory naming conventions are summarized in
`docs/directory_structure.md`. The v0.1 pre-release keeps the legacy-compatible
data directory names below as the workflow contract:

- `data/raw_data/`: external inputs and prepared non-git payloads.
- `data/int_data/`: intermediate outputs from upstream stages.
- `data/sim_data/`: simulated PAS/BAM/expression and CPRSB inputs.
- `data/result/`: downstream benchmark, analysis, plotting-source, and resource
  profile outputs.

Stages 01-04 use `../../data` paths from their stage directories, so a full
upstream run requires `data/` at the repository root. Stages 05-06 also accept
`APABENCHMARK_FINAL_ROOT` or `--data-root`.

## External Inputs

For release tracking, source URLs, versions, and checksums are listed in
`docs/external_inputs_manifest.md`.

Required annotation inputs:

```text
data/raw_data/annotations/
  gencode.vM25.polyAs.gtf
  gencode.vM25.annotation.bed
  atlas.clusters.3.0.GRCm38.GENCODE_M25.bed.gz
  MousePas/mm10.PAS.main.tsv
  gencode.v40.polyAs.gtf
  gencode.v40.annotation.bed
  atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz
  HumanPas/hg38.PAS.main.tsv
```

Additional full-run inputs:

```text
data/raw_data/bam_for_peak_extraction/{sample}.bam
data/raw_data/bam_for_peak_extraction/{sample}.bam.bai
data/raw_data/tool_reference/*
scripts/04_benchmark/sif/*.sif
APABENCHMARK_MM10_GENOME_FASTA=/path/to/GRCm38.p6.genome.fa
```

Raw-data benchmark panels require separately prepared raw benchmark outputs and
raw ground-truth directories under `data/result/benchmark/raw_bam`,
`data/raw_data/raw_bam`, and `data/result/benchmark/analysis/raw_pair_manifest`.

## Environment Contract

```bash
conda env create -f scripts/envs/apasim.yml
conda env create -f scripts/envs/dexseq.yml
conda env create -f scripts/envs/scmapa_r.yml
conda env create -f scripts/envs/benchmark_plotting.yml
```

The workflow also requires `snakemake`, `samtools`, `bedtools`, `umi_tools`, and
Singularity or Apptainer. Stage 01 uses the `apasim` environment through the full
runner.

## Stage Entrypoints

### 01 Generate Annotation

```bash
cd scripts/01_generate_annotation
bash run_annotation.sh
```

Primary outputs:

```text
data/int_data/annotations/mouse_integrated_pas.bed
data/int_data/annotations/human_integrated_pas.bed
data/int_data/annotations/mouse_pas_500.bed
data/int_data/annotations/mouse_pas_1000.bed
data/int_data/annotations/mouse_pas_single.bed
data/int_data/annotations/human_pas_500.bed
data/int_data/annotations/human_pas_1000.bed
data/int_data/annotations/human_pas_single.bed
data/int_data/annotations/supplementary_annotation_stage_counts.tsv
```

### 02 Peak Analysis

```bash
cd scripts/02_peak_analysis
snakemake --cores 8 --use-conda
```

Consumes `mouse_pas_500.bed` and raw BAMs, then produces processed BAMs, peak
tables, and peak-model summaries under `data/int_data`.

### 03 Data Simulation

```bash
cd scripts/03_data_simulation
APABENCHMARK_MM10_GENOME_FASTA=/path/to/GRCm38.p6.genome.fa \
snakemake --cores 8 --use-conda
```

Consumes integrated mouse PAS annotations and peak tables, then produces
simulated PAS sets, simulated BAMs, expression matrices, and CPRSB resource
benchmark inputs under `data/sim_data`.

### 04 Benchmark Tool Runs

```bash
cd scripts/04_benchmark
bash download_sif.sh
bash generate_sample_config.sh /path/to/apabenchmark_final ./snakemake_profile/sample.yaml
snakemake --cores 8 --profile snakemake_profile \
  --singularity-args "--bind /path/to/apabenchmark_final:/path/to/apabenchmark_final --cpus 8"
```

Consumes simulated BAMs, tool references, and SIF images; produces tool PAS
outputs and resource profiles under `data/result/benchmark` and
`data/result/cp_resource`.

### 05 Benchmark Analysis

```bash
APABENCHMARK_FINAL_ROOT=/path/to/apabenchmark_final \
APABENCHMARK_PIPELINE_MODE=both \
conda run -n apasim \
  snakemake -s scripts/05_benchmark_analysis/Snakefile --cores 8
```

Produces performance tables, DaPars2-format inputs, differential APA outputs,
and optional raw-mode analysis outputs under `data/result`.

### 06 Plot Scripts

```bash
conda run -n benchmark_plotting \
  python scripts/06_plotting/sim_data_performance/scripts/run_all.py \
    --data-root /path/to/apabenchmark_final
```

Additional connected plotting topics:

```bash
conda run -n benchmark_plotting \
  python scripts/06_plotting/raw_data_performance/scripts/run_all.py \
    --data-root /path/to/apabenchmark_final --raw-run-id default

conda run -n benchmark_plotting \
  python scripts/06_plotting/peak_model_params/scripts/run_all.py

APABENCHMARK_CP_RESOURCE_ROOT=/path/to/apabenchmark_final/data/result/cp_resource/pf_bam \
conda run -n benchmark_plotting \
  python scripts/06_plotting/computation_resource_consumption/scripts/run_all.py
```

## Full-Flow Runner

Use the top-level coordinator:

```bash
./scripts/run_full_repro.sh \
  --data-root /path/to/apabenchmark_final \
  --mm10-fasta /path/to/GRCm38.p6.genome.fa \
  --cores 8
```

Useful options:

```text
--skip-01 ... --skip-05
--skip-plots
--with-raw-plots
--pipeline-mode sim|raw|both|raw_parallel
--dry-run
```

## Known Continuity Notes

- Stage 04 needs valid SIF images and a Singularity/Apptainer runtime.
- Raw-data analysis depends on externally prepared raw benchmark outputs and raw
  ground-truth manifests.
- Legacy plotting topics under `pas_detect_performance`,
  `pas_quantification_performance`, and `apa_detect_performance` are preserved as
  reference script topics; the connected simulation plotting path is
  `sim_data_performance`.
