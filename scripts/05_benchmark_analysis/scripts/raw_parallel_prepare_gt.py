#!/usr/bin/env python3
import argparse
import os
from raw_parallel_common import read_tsv, write_tsv, run_cmd, touch_done, ensure_parent


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run-root", required=True)
    p.add_argument("--script-root", required=True)
    p.add_argument("--conda-exe", required=True)
    p.add_argument("--scmapa-env", required=True)
    p.add_argument("--done-file", required=True)
    args = p.parse_args()

    pair_manifest = os.path.join(args.run_root, "raw_pair_manifest", "sample_pair_manifest.tsv")
    pair_rows = read_tsv(pair_manifest)

    raw_root = os.path.join(args.run_root, "raw_dapars2_parallel")
    gt_dir = os.path.join(raw_root, "ground_truth")
    de_dir = os.path.join(raw_root, "differential_apa", "ground_truth")
    log_dir = os.path.join(raw_root, "logs", "ground_truth")
    os.makedirs(raw_root, exist_ok=True)

    multi_group_converter = os.path.join(args.script_root, "scmapa", "multi_group_PAS_to_DaPars2_converter.R")
    scmapa_script = os.path.join(args.script_root, "scmapa", "parameterized_scmapa_analysis.R")

    failures = []
    for pair in pair_rows:
        sample = pair["sample"]
        pair_label = pair["pair_label"]
        g1_clean = pair["group1_clean"]
        g2_clean = pair["group2_clean"]
        gt_pair_expr = os.path.join(args.run_root, "raw_pair_inputs", "ground_truth", sample, pair_label, "gt_expr_pair.tsv")
        gt_dapars = os.path.join(gt_dir, sample, pair_label, "dapars2_format.txt")
        gt_de_out = os.path.join(de_dir, sample, pair_label)
        gt_apa = os.path.join(gt_de_out, f"{g1_clean}_vs_{g2_clean}", "APA_results.txt")

        try:
            ensure_parent(gt_dapars)
            run_cmd([
                args.conda_exe, "run", "-n", args.scmapa_env, "Rscript", multi_group_converter,
                "-i", gt_pair_expr, "-o", gt_dapars, "--min_exp", "0", "--min_cells", "0",
            ], log_path=os.path.join(log_dir, sample, pair_label, "convert.log"))
            if not os.path.exists(gt_dapars):
                raise RuntimeError(f"missing_gt_dapars:{gt_dapars}")
            os.makedirs(gt_de_out, exist_ok=True)
            run_cmd([
                args.conda_exe, "run", "-n", args.scmapa_env, "Rscript", scmapa_script,
                "-i", gt_dapars, "-o", gt_de_out,
                "--comparison_mode", "single", "--group1", g1_clean, "--group2", g2_clean,
                "--NAcutoff", "0", "--CPMcutoff_L", "0.001", "--CPMcutoff_S", "0.001",
                "--coverageCutoff", "0.001", "--ORcutoff", "0.1", "--adPval", "0.05",
            ], log_path=os.path.join(log_dir, sample, pair_label, "scmapa.log"))
            if not os.path.exists(gt_apa):
                raise RuntimeError(f"missing_gt_apa:{gt_apa}")
        except Exception as e:
            failures.append({"sample": sample, "pair_label": pair_label, "tool": "ground_truth", "stage": "prepare_gt", "message": str(e)})

    write_tsv(
        os.path.join(raw_root, "summary", "raw_dapars2_parallel_gt_failures.tsv"),
        failures,
        ["sample", "pair_label", "tool", "stage", "message"],
    )
    touch_done(args.done_file)


if __name__ == "__main__":
    main()
