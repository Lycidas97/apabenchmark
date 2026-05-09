#!/usr/bin/env python3
import argparse
import os
from raw_parallel_common import read_tsv, write_tsv, run_cmd, touch_done


def concat_tsv(files, out_path, header_fields=None):
    files = [f for f in files if os.path.exists(f) and os.path.getsize(f) > 0]
    if not files:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w", encoding="utf-8", newline="") as out_fh:
            if header_fields:
                out_fh.write("\t".join(header_fields) + "\n")
        return False
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    header = None
    with open(out_path, "w", encoding="utf-8", newline="") as out_fh:
        for f in files:
            with open(f, "r", encoding="utf-8", newline="") as in_fh:
                lines = in_fh.readlines()
            if not lines:
                continue
            if header is None:
                header = lines[0]
                out_fh.write(header)
            for line in lines[1:]:
                out_fh.write(line)
    return True


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run-root", required=True)
    p.add_argument("--script-root", required=True)
    p.add_argument("--conda-exe", required=True)
    p.add_argument("--dexseq-env", required=True)
    p.add_argument("--scmapa-env", required=True)
    p.add_argument("--done-file", required=True)
    args = p.parse_args()

    pair_manifest = os.path.join(args.run_root, "raw_pair_manifest", "sample_pair_manifest.tsv")
    pair_rows = read_tsv(pair_manifest)

    raw_root = os.path.join(args.run_root, "raw_dapars2_parallel")
    conv_dir = os.path.join(raw_root, "converted_inputs")
    de_dir = os.path.join(raw_root, "differential_apa")
    perf_dir = os.path.join(raw_root, "performance")
    log_dir = os.path.join(raw_root, "logs")
    summary_dir = os.path.join(raw_root, "summary")
    os.makedirs(summary_dir, exist_ok=True)

    scmapa_script = os.path.join(args.script_root, "scmapa", "parameterized_scmapa_analysis.R")
    perf_script = os.path.join(args.script_root, "scripts", "calculate_dapars2_de_apa_performance.R")

    failures = []
    for pair in pair_rows:
        sample = pair["sample"]
        pair_label = pair["pair_label"]
        g1_clean = pair["group1_clean"]
        g2_clean = pair["group2_clean"]
        sample_pair_id = f"{sample}__{pair_label}"

        gt_apa = os.path.join(de_dir, "ground_truth", sample, pair_label, f"{g1_clean}_vs_{g2_clean}", "APA_results.txt")
        if not os.path.exists(gt_apa):
            failures.append({"sample": sample, "pair_label": pair_label, "tool": "ground_truth", "stage": "eval", "message": f"missing_gt_apa:{gt_apa}"})
            continue

        for tool in ["sierra", "scmapa", "dapars2"]:
            try:
                tool_dapars = os.path.join(conv_dir, tool, sample, pair_label, "dapars2_format.txt")
                tool_de_out = os.path.join(de_dir, tool, sample, pair_label)
                os.makedirs(tool_de_out, exist_ok=True)
                run_cmd([
                    args.conda_exe, "run", "-n", args.scmapa_env, "Rscript", scmapa_script,
                    "-i", tool_dapars, "-o", tool_de_out,
                    "--comparison_mode", "single", "--group1", g1_clean, "--group2", g2_clean,
                    "--NAcutoff", "0", "--CPMcutoff_L", "0.001", "--CPMcutoff_S", "0.001",
                    "--coverageCutoff", "0.001", "--ORcutoff", "0.1", "--adPval", "0.05",
                ], log_path=os.path.join(log_dir, tool, sample, pair_label, "scmapa.log"))

                pd_apa = os.path.join(tool_de_out, f"{g1_clean}_vs_{g2_clean}", "APA_results.txt")
                perf_out = os.path.join(perf_dir, tool, sample, pair_label, "dapars2_performance.tsv")
                os.makedirs(os.path.dirname(perf_out), exist_ok=True)
                run_cmd([
                    args.conda_exe, "run", "-n", args.dexseq_env, "Rscript", perf_script,
                    "--pd_apa", pd_apa,
                    "--gt_apa", gt_apa,
                    "--tool", tool,
                    "--sample", sample_pair_id,
                    "--match_type", f"{g1_clean}_vs_{g2_clean}",
                    "--filter_type_1", "scmapa_adPval_0.05",
                    "--filter_type_2", "scmapa_or_0.1",
                    "--output", perf_out,
                ], log_path=os.path.join(log_dir, tool, sample, pair_label, "performance.log"))
            except Exception as e:
                failures.append({"sample": sample, "pair_label": pair_label, "tool": tool, "stage": "eval", "message": str(e)})

    perf_files = []
    for pair in pair_rows:
        sample = pair["sample"]
        pair_label = pair["pair_label"]
        for tool in ["sierra", "scmapa", "dapars2"]:
            perf_files.append(os.path.join(perf_dir, tool, sample, pair_label, "dapars2_performance.tsv"))

    concat_tsv(
        perf_files,
        os.path.join(summary_dir, "raw_dapars2_parallel_metrics.tsv"),
        header_fields=["sample", "tool", "protocol", "match_type", "filter_type_1", "filter_type_2", "precision", "recall", "f1", "gt_count", "pd_count"],
    )

    write_tsv(
        os.path.join(summary_dir, "raw_dapars2_parallel_missing_failed.tsv"),
        failures,
        ["sample", "pair_label", "tool", "stage", "message"],
    )

    touch_done(args.done_file)


if __name__ == "__main__":
    main()
