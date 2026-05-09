#!/usr/bin/env python3
import argparse
import os
from raw_parallel_common import read_tsv, write_tsv, run_cmd, touch_done, ensure_parent


FAIL_HEADERS = ["sample", "pair_label", "tool", "stage", "message"]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run-root", required=True)
    p.add_argument("--script-root", required=True)
    p.add_argument("--conda-exe", required=True)
    p.add_argument("--dexseq-env", required=True)
    p.add_argument("--split-num", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--pair-label", required=True)
    p.add_argument("--tool", required=True)
    p.add_argument("--state-file", required=True)
    p.add_argument("--failure-file", required=True)
    args = p.parse_args()

    pair_manifest = os.path.join(args.run_root, "raw_pair_manifest", "sample_pair_manifest.tsv")
    pair_rows = read_tsv(pair_manifest)
    pair = next((r for r in pair_rows if r["sample"] == args.sample and r["pair_label"] == args.pair_label), None)
    if pair is None:
        raise RuntimeError(f"pair_not_found:{args.sample}:{args.pair_label}")

    sample = pair["sample"]
    pair_label = pair["pair_label"]
    sample_pair_id = f"{sample}__{pair_label}"

    raw_pair_inputs = os.path.join(args.run_root, "raw_pair_inputs")
    raw_perf_dir = os.path.join(args.run_root, "raw_performance")
    out_dir = os.path.join(raw_perf_dir, "by_sample_pair_tool")
    log_dir = os.path.join(raw_perf_dir, "logs")

    gt_pair_expr = os.path.join(raw_pair_inputs, "ground_truth", sample, pair_label, "gt_expr_pair.tsv")
    gt_pair_pas = os.path.join(raw_pair_inputs, "ground_truth", sample, pair_label, "gt_pas_pair.bed")
    pd_pair_expr = os.path.join(raw_pair_inputs, args.tool, sample, pair_label, "pd_pas_counts.tsv")
    pd_pair_pas = os.path.join(raw_pair_inputs, args.tool, sample, pair_label, "pd_pas.bed")

    stage10_perf_script = os.path.join(args.script_root, "scripts", "10_stage_pas_match_quantify.R")
    stage20_perf_script = os.path.join(args.script_root, "scripts", "20_stage_apa_and_te.R")
    stage30_perf_script = os.path.join(args.script_root, "scripts", "30_stage_export_results.R")

    failures = []
    try:
        for need in [gt_pair_expr, gt_pair_pas, pd_pair_expr, pd_pair_pas, stage10_perf_script, stage20_perf_script, stage30_perf_script]:
            if not os.path.exists(need):
                raise RuntimeError(f"missing_input:{need}")

        run_cmd(
            [
                args.conda_exe,
                "run",
                "-n",
                args.dexseq_env,
                "Rscript",
                stage10_perf_script,
                "--pd_pas",
                pd_pair_pas,
                "--pd_mtx",
                pd_pair_expr,
                "--gt_pas",
                gt_pair_pas,
                "--gt_mtx",
                gt_pair_expr,
                "--tool",
                args.tool,
                "--sample",
                sample_pair_id,
                "--output_dir",
                out_dir,
                "--dexseq_split_num",
                str(args.split_num),
            ],
            log_path=os.path.join(log_dir, args.tool, sample, pair_label, "run_performance_stage10.log"),
        )

        stage10_bundle = os.path.join(out_dir, args.tool, f"{sample_pair_id}_stage10_bundle.rds")
        if not os.path.exists(stage10_bundle):
            raise RuntimeError(f"missing_stage10_bundle:{stage10_bundle}")

        run_cmd(
            [
                args.conda_exe,
                "run",
                "-n",
                args.dexseq_env,
                "Rscript",
                stage20_perf_script,
                "--stage10_bundle",
                stage10_bundle,
            ],
            log_path=os.path.join(log_dir, args.tool, sample, pair_label, "run_performance_stage20.log"),
        )

        stage20_bundle = os.path.join(out_dir, args.tool, f"{sample_pair_id}_stage20_bundle.rds")
        if not os.path.exists(stage20_bundle):
            raise RuntimeError(f"missing_stage20_bundle:{stage20_bundle}")

        run_cmd(
            [
                args.conda_exe,
                "run",
                "-n",
                args.dexseq_env,
                "Rscript",
                stage30_perf_script,
                "--stage20_bundle",
                stage20_bundle,
            ],
            log_path=os.path.join(log_dir, args.tool, sample, pair_label, "run_performance_stage30.log"),
        )
    except Exception as e:
        failures.append(
            {
                "sample": sample,
                "pair_label": pair_label,
                "tool": args.tool,
                "stage": "calculate_benchmark_performance",
                "message": str(e),
            }
        )

    # Keep deterministic outputs for downstream aggregation.
    for suffix in ["match_performance", "pas_quantify_performance", "de_apa_performance"]:
        out_path = os.path.join(out_dir, args.tool, f"{sample_pair_id}_{suffix}.tsv")
        if not os.path.exists(out_path):
            ensure_parent(out_path)
            with open(out_path, "w", encoding="utf-8") as fh:
                fh.write("")

    write_tsv(args.failure_file, failures, FAIL_HEADERS)
    touch_done(args.state_file)


if __name__ == "__main__":
    main()
