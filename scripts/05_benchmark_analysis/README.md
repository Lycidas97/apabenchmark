# 05_benchmark_analysis

Postprocess benchmark workflow for the `apabenchmark_final` data layout.

For release reproduction, prefer the top-level `scripts/run_full_repro.sh`
entrypoint. Use this README for Stage 05 debugging, raw-mode reruns, or
developer-level DAG checks.

## Inputs (final layout)
- Tool results: `data/result/benchmark/sim_bam/{tool}/{sample}/...`
- Ground truth matrix: `data/sim_data/bam/{sample}.bam.expr.tsv`
- Ground truth PAS: `data/sim_data/sim_pas/{genome}_sim_pas_{gene}_rep{rep}.bed`
- Raw pair manifest: `data/result/benchmark/analysis/raw_pair_manifest/{sample_top3_celltypes.tsv,sample_pair_manifest.tsv}`
- Raw source data: `data/result/benchmark/raw_bam/{tool}/{sample}/...` and `data/raw_data/raw_bam/{sample}/...`

## Outputs
- Performance: `data/result/performance/{tool}/{sample}_*.tsv`
- Quality report: `data/result/performance/{tool}/{sample}_quality_report.tsv`
- DaPars2 format: `data/result/dapars2_format/...`
- Differential APA: `data/result/differential_apa/...`
- Ground truth differential APA: `data/result/differential_apa_ground_truth/...`
- Raw benchmark run root (when `pipeline_mode=raw|both`):
  `data/result/raw/{run_id}/raw_{pair_manifest,pair_inputs,qc,performance,dapars2_format}/...`

## Run
```bash
cd /path/to/apabenchmark_final
conda run -n apasim \
  snakemake -s scripts/05_benchmark_analysis/Snakefile --cores 8
```

## Run raw benchmark pipeline
```bash
cd /path/to/apabenchmark_final
APABENCHMARK_PIPELINE_MODE=raw \
conda run -n apasim \
  snakemake -s scripts/05_benchmark_analysis/Snakefile --cores 8
```

Raw mode and `raw_parallel` mode share the same pair-input preparation step
(`run_raw_top_samples_split5.py --prepare-only`), then execute per `sample x pair x tool`
combination rules in Snakemake (no dependency on `scripts/04_benchmark`).

- `pipeline_mode=raw`: pair-level sierra benchmark + pair/tool-level DaPars2-parallel branch, then aggregate summaries.
- `pipeline_mode=raw_parallel`: pair/tool-level DaPars2-parallel branch only, then aggregate summaries.

`pipeline_mode=both` will run both sim and raw targets.

This workflow uses fixed named conda environments inside rules (`dexseq` and `scmapa_r`)
via `conda run -n ...`. Do not pass `--use-conda`.

## Full DAG dry-run
```bash
conda run -n apasim \
  snakemake -s scripts/05_benchmark_analysis/Snakefile -n --cores 1
```

## Optional env vars
- `APABENCHMARK_FINAL_ROOT` (project root containing `data/`)
- `APABENCHMARK_BENCHMARK5_DIR` (override workflow directory)
- `APABENCHMARK_GENE_TAG` (default auto-detect)
- `APABENCHMARK_DAPARS2_EXCLUDE_FILE` (optional exclude list for single-group DaPars2)
- `APABENCHMARK_CONDA_EXE` (default: `conda` on `PATH`)
- `APABENCHMARK_DEXSEQ_ENV` (default: `dexseq`)
- `APABENCHMARK_DEXSEQ_SPLIT_NUM` (default: `6` for `sim`; default follows `raw_split_num` for `raw/raw_parallel`)
- `APABENCHMARK_SCMAPA_ENV` (default: `scmapa_r`)
- `APABENCHMARK_BEDTOOLS_BIN` (default: `bedtools`; can also be a full binary path or bin directory)
- `APABENCHMARK_PIPELINE_MODE` (`sim`/`raw`/`both`/`raw_parallel`, default: `sim`)
- `APABENCHMARK_RAW_RUN_ID` (default: `default`, sanitized to `[A-Za-z0-9_.-]`)
- `APABENCHMARK_RAW_OUTPUT_BASE` (default: `data/result/raw`)
- `APABENCHMARK_RAW_TOP_N` (default: inferred from `raw_pair_manifest/sample_top3_celltypes.tsv`, i.e. full available samples)
- `APABENCHMARK_RAW_SPLIT_NUM` (default: `5`)
- `APABENCHMARK_RAW_PAIR_RANK_MAX` (default: inferred max `pair_rank` from `raw_pair_manifest/sample_pair_manifest.tsv`, i.e. full available pairs)
