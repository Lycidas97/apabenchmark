#!/usr/bin/env python3
import argparse
import os
from raw_parallel_common import read_tsv, write_tsv, touch_done, extract_dapars_pair_input


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--project-root", required=True)
    p.add_argument("--run-root", required=True)
    p.add_argument("--done-file", required=True)
    args = p.parse_args()

    raw_bam_root = os.path.join(args.project_root, "data", "result", "benchmark", "raw_bam")
    pair_manifest = os.path.join(args.run_root, "raw_pair_manifest", "sample_pair_manifest.tsv")
    pair_rows = read_tsv(pair_manifest)

    raw_root = os.path.join(args.run_root, "raw_dapars2_parallel")
    conv_dir = os.path.join(raw_root, "converted_inputs")

    failures = []
    for pair in pair_rows:
        sample = pair["sample"]
        pair_label = pair["pair_label"]
        group1 = pair["group1"]
        group2 = pair["group2"]
        g1_clean = pair["group1_clean"]
        g2_clean = pair["group2_clean"]

        for tool in ["scmapa", "dapars2"]:
            src = os.path.join(raw_bam_root, tool, sample, "DaPars2_Result_all_chrs.txt")
            out = os.path.join(conv_dir, tool, sample, pair_label, "dapars2_format.txt")
            try:
                extract_dapars_pair_input(src, group1, group2, g1_clean, g2_clean, out)
            except Exception as e:
                failures.append({"sample": sample, "pair_label": pair_label, "tool": tool, "stage": "prepare_dapars", "message": str(e)})

    write_tsv(
        os.path.join(raw_root, "summary", "raw_dapars2_parallel_dapars_failures.tsv"),
        failures,
        ["sample", "pair_label", "tool", "stage", "message"],
    )
    touch_done(args.done_file)


if __name__ == "__main__":
    main()
