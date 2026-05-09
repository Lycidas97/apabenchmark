# scripts directory

This directory contains executable scripts used by the workflow.

Current scripts:
- `calculate_benchmark_performance.R`
- `10_stage_pas_match_quantify.R`
- `20_stage_apa_and_te.R`
- `30_stage_export_results.R`
- `common_performance_functions.R`
- `calculate_te_gap.py`
- `calculate_dapars2_de_apa_performance.R`
- `generate_dapars2_single_group_exclude.py`
- `monitor_retry_process.sh`
- `raw_parallel_aggregate.py`
- `raw_parallel_common.py`
- `raw_parallel_eval_one.py`
- `raw_parallel_prepare_dapars_one.py`
- `raw_parallel_prepare_dapars.py`
- `raw_parallel_prepare_gt_one.py`
- `raw_parallel_prepare_gt.py`
- `raw_parallel_prepare_nondapars_one.py`
- `raw_parallel_prepare_nondapars.py`
- `raw_parallel_run_eval.py`
- `raw_performance_aggregate.py`
- `raw_performance_eval_one.py`
- `run_raw_top_samples_split5.py`

Path conventions:
- Use `APABENCHMARK_FINAL_ROOT` (or `APABENCHMARK_ROOT`) to resolve data paths.
- `Snakefile` resolves script paths relative to this workflow directory (or `APABENCHMARK_BENCHMARK5_DIR` if provided).
- `APABENCHMARK_BEDTOOLS_BIN` can override the `bedtools` binary path (default uses `bedtools` from `PATH`).
