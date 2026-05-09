#!/usr/bin/env python3
import argparse
import os
from raw_parallel_common import read_tsv, write_tsv, run_cmd, touch_done, ensure_parent


METRIC_HEADERS = [
    "sample",
    "tool",
    "protocol",
    "match_type",
    "filter_type_1",
    "filter_type_2",
    "precision",
    "recall",
    "f1",
    "gt_count",
    "pd_count",
]


def write_metric_header(path: str) -> None:
    ensure_parent(path)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(METRIC_HEADERS) + "\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run-root", required=True)
    p.add_argument("--script-root", required=True)
    p.add_argument("--conda-exe", required=True)
    p.add_argument("--dexseq-env", required=True)
    p.add_argument("--scmapa-env", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--pair-label", required=True)
    p.add_argument("--tool", required=True)
    p.add_argument("--perf-out", required=True)
    p.add_argument("--failure-file", required=True)
    p.add_argument("--state-file", required=True)
    args = p.parse_args()

    pair_manifest = os.path.join(args.run_root, "raw_pair_manifest", "sample_pair_manifest.tsv")
    pair_rows = read_tsv(pair_manifest)
    pair = next((r for r in pair_rows if r["sample"] == args.sample and r["pair_label"] == args.pair_label), None)
    if pair is None:
        raise RuntimeError(f"pair_not_found:{args.sample}:{args.pair_label}")

    sample = pair["sample"]
    pair_label = pair["pair_label"]
    g1_clean = pair["group1_clean"]
    g2_clean = pair["group2_clean"]
    sample_pair_id = f"{sample}__{pair_label}"

    raw_root = os.path.join(args.run_root, "raw_dapars2_parallel")
    conv_dir = os.path.join(raw_root, "converted_inputs")
    de_dir = os.path.join(raw_root, "differential_apa")
    log_dir = os.path.join(raw_root, "logs")
    scmapa_script = os.path.join(args.script_root, "scmapa", "parameterized_scmapa_analysis.R")
    perf_script = os.path.join(args.script_root, "scripts", "calculate_dapars2_de_apa_performance.R")

    failures = []
    write_metric_header(args.perf_out)
    gt_apa = os.path.join(de_dir, "ground_truth", sample, pair_label, f"{g1_clean}_vs_{g2_clean}", "APA_results.txt")
    if not os.path.exists(gt_apa):
        failures.append(
            {
                "sample": sample,
                "pair_label": pair_label,
                "tool": "ground_truth",
                "stage": "eval",
                "message": f"missing_gt_apa:{gt_apa}",
            }
        )
        write_tsv(args.failure_file, failures, ["sample", "pair_label", "tool", "stage", "message"])
        touch_done(args.state_file)
        return

    try:
        tool_dapars = os.path.join(conv_dir, args.tool, sample, pair_label, "dapars2_format.txt")
        tool_de_out = os.path.join(de_dir, args.tool, sample, pair_label)
        os.makedirs(tool_de_out, exist_ok=True)
        run_cmd(
            [
                args.conda_exe,
                "run",
                "-n",
                args.scmapa_env,
                "Rscript",
                scmapa_script,
                "-i",
                tool_dapars,
                "-o",
                tool_de_out,
                "--comparison_mode",
                "single",
                "--group1",
                g1_clean,
                "--group2",
                g2_clean,
                "--NAcutoff",
                "0",
                "--CPMcutoff_L",
                "0.001",
                "--CPMcutoff_S",
                "0.001",
                "--coverageCutoff",
                "0.001",
                "--ORcutoff",
                "0.1",
                "--adPval",
                "0.05",
            ],
            log_path=os.path.join(log_dir, args.tool, sample, pair_label, "scmapa.log"),
        )

        pd_apa = os.path.join(tool_de_out, f"{g1_clean}_vs_{g2_clean}", "APA_results.txt")
        run_cmd(
            [
                args.conda_exe,
                "run",
                "-n",
                args.dexseq_env,
                "Rscript",
                perf_script,
                "--pd_apa",
                pd_apa,
                "--gt_apa",
                gt_apa,
                "--tool",
                args.tool,
                "--sample",
                sample_pair_id,
                "--match_type",
                f"{g1_clean}_vs_{g2_clean}",
                "--filter_type_1",
                "scmapa_adPval_0.05",
                "--filter_type_2",
                "scmapa_or_0.1",
                "--output",
                args.perf_out,
            ],
            log_path=os.path.join(log_dir, args.tool, sample, pair_label, "performance.log"),
        )
    except Exception as e:
        failures.append(
            {
                "sample": sample,
                "pair_label": pair_label,
                "tool": args.tool,
                "stage": "eval",
                "message": str(e),
            }
        )

    write_tsv(args.failure_file, failures, ["sample", "pair_label", "tool", "stage", "message"])
    touch_done(args.state_file)


if __name__ == "__main__":
    main()
