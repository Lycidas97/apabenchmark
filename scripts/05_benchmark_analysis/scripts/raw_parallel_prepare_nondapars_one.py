#!/usr/bin/env python3
import argparse
import os
from raw_parallel_common import write_tsv, run_cmd, touch_done


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run-root", required=True)
    p.add_argument("--script-root", required=True)
    p.add_argument("--conda-exe", required=True)
    p.add_argument("--scmapa-env", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--pair-label", required=True)
    p.add_argument("--tool", required=True)
    p.add_argument("--state-file", required=True)
    p.add_argument("--failure-file", required=True)
    args = p.parse_args()

    sample = args.sample
    pair_label = args.pair_label

    raw_root = os.path.join(args.run_root, "raw_dapars2_parallel")
    conv_dir = os.path.join(raw_root, "converted_inputs", args.tool)
    log_dir = os.path.join(raw_root, "logs", args.tool)
    integrated_converter = os.path.join(args.script_root, "scmapa", "integrated_pas_to_dapars2_converter.R")

    pd_pair_expr = os.path.join(args.run_root, "raw_pair_inputs", args.tool, sample, pair_label, "pd_pas_counts.tsv")
    pd_pair_pas = os.path.join(args.run_root, "raw_pair_inputs", args.tool, sample, pair_label, "pd_pas.bed")
    gt_pair_pas = os.path.join(args.run_root, "raw_pair_inputs", "ground_truth", sample, pair_label, "gt_pas_pair.bed")
    gt_pair_expr = os.path.join(args.run_root, "raw_pair_inputs", "ground_truth", sample, pair_label, "gt_expr_pair.tsv")
    out_path = os.path.join(conv_dir, sample, pair_label, "dapars2_format.txt")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    failures = []
    try:
        run_cmd(
            [
                args.conda_exe,
                "run",
                "-n",
                args.scmapa_env,
                "Rscript",
                integrated_converter,
                "--pas_counts",
                pd_pair_expr,
                "--pas_coords",
                pd_pair_pas,
                "--pas_annotation",
                gt_pair_pas,
                "--celltype_info",
                gt_pair_expr,
                "-o",
                out_path,
                "--min_exp",
                "0",
                "--min_cells",
                "0",
            ],
            log_path=os.path.join(log_dir, sample, pair_label, "convert.log"),
        )
    except Exception as e:
        failures.append(
            {
                "sample": sample,
                "pair_label": pair_label,
                "tool": args.tool,
                "stage": "prepare_nondapars",
                "message": str(e),
            }
        )

    write_tsv(args.failure_file, failures, ["sample", "pair_label", "tool", "stage", "message"])
    touch_done(args.state_file)


if __name__ == "__main__":
    main()
