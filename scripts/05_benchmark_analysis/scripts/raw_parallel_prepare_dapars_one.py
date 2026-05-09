#!/usr/bin/env python3
import argparse
import os
from raw_parallel_common import read_tsv, write_tsv, touch_done, extract_dapars_pair_input


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--project-root", required=True)
    p.add_argument("--run-root", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--pair-label", required=True)
    p.add_argument("--tool", required=True, choices=["scmapa", "dapars2"])
    p.add_argument("--state-file", required=True)
    p.add_argument("--failure-file", required=True)
    args = p.parse_args()

    pair_manifest = os.path.join(args.run_root, "raw_pair_manifest", "sample_pair_manifest.tsv")
    pair_rows = read_tsv(pair_manifest)
    pair = next((r for r in pair_rows if r["sample"] == args.sample and r["pair_label"] == args.pair_label), None)
    if pair is None:
        raise RuntimeError(f"pair_not_found:{args.sample}:{args.pair_label}")

    raw_bam_root = os.path.join(args.project_root, "data", "result", "benchmark", "raw_bam")
    raw_root = os.path.join(args.run_root, "raw_dapars2_parallel")
    conv_dir = os.path.join(raw_root, "converted_inputs")

    sample = pair["sample"]
    pair_label = pair["pair_label"]
    group1 = pair["group1"]
    group2 = pair["group2"]
    g1_clean = pair["group1_clean"]
    g2_clean = pair["group2_clean"]

    src = os.path.join(raw_bam_root, args.tool, sample, "DaPars2_Result_all_chrs.txt")
    out = os.path.join(conv_dir, args.tool, sample, pair_label, "dapars2_format.txt")

    failures = []
    try:
        extract_dapars_pair_input(src, group1, group2, g1_clean, g2_clean, out)
    except Exception as e:
        failures.append(
            {
                "sample": sample,
                "pair_label": pair_label,
                "tool": args.tool,
                "stage": "prepare_dapars",
                "message": str(e),
            }
        )

    write_tsv(args.failure_file, failures, ["sample", "pair_label", "tool", "stage", "message"])
    touch_done(args.state_file)


if __name__ == "__main__":
    main()
