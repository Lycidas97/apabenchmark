#!/usr/bin/env python3
import argparse
import os
from raw_parallel_common import read_tsv, write_tsv


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
FAIL_HEADERS = ["sample", "pair_label", "tool", "stage", "message"]


def concat_metrics(files, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as out_fh:
        out_fh.write("\t".join(METRIC_HEADERS) + "\n")
        for f in files:
            if not os.path.exists(f):
                continue
            with open(f, "r", encoding="utf-8", newline="") as in_fh:
                lines = in_fh.readlines()
            if len(lines) <= 1:
                continue
            for line in lines[1:]:
                out_fh.write(line)


def merge_failures(files):
    rows = []
    for f in files:
        if not os.path.exists(f):
            continue
        rows.extend(read_tsv(f))
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run-root", required=True)
    p.add_argument("--done-file", required=True)
    args = p.parse_args()

    pair_manifest = os.path.join(args.run_root, "raw_pair_manifest", "sample_pair_manifest.tsv")
    pair_rows = read_tsv(pair_manifest)

    raw_root = os.path.join(args.run_root, "raw_dapars2_parallel")
    perf_dir = os.path.join(raw_root, "performance")
    summary_dir = os.path.join(raw_root, "summary")
    fail_dir = os.path.join(raw_root, "failures")
    os.makedirs(summary_dir, exist_ok=True)

    perf_files = []
    eval_fail_files = []
    gt_fail_files = []
    nondapars_fail_files = []
    dapars_fail_files = []
    tools = []
    if os.path.isdir(perf_dir):
        tools = sorted([name for name in os.listdir(perf_dir) if os.path.isdir(os.path.join(perf_dir, name))])

    nondapars_tools = []
    nondapars_root = os.path.join(fail_dir, "prepare_nondapars")
    if os.path.isdir(nondapars_root):
        nondapars_tools = sorted([name for name in os.listdir(nondapars_root) if os.path.isdir(os.path.join(nondapars_root, name))])

    dapars_tools = []
    dapars_root = os.path.join(fail_dir, "prepare_dapars")
    if os.path.isdir(dapars_root):
        dapars_tools = sorted([name for name in os.listdir(dapars_root) if os.path.isdir(os.path.join(dapars_root, name))])

    for pair in pair_rows:
        sample = pair["sample"]
        pair_label = pair["pair_label"]
        for tool in tools:
            perf_files.append(os.path.join(perf_dir, tool, sample, pair_label, "dapars2_performance.tsv"))
            eval_fail_files.append(os.path.join(fail_dir, "eval", tool, sample, pair_label, "failure.tsv"))
        gt_fail_files.append(os.path.join(fail_dir, "prepare_gt", sample, pair_label, "failure.tsv"))
        for tool in nondapars_tools:
            nondapars_fail_files.append(os.path.join(fail_dir, "prepare_nondapars", tool, sample, pair_label, "failure.tsv"))
        for tool in dapars_tools:
            dapars_fail_files.append(os.path.join(fail_dir, "prepare_dapars", tool, sample, pair_label, "failure.tsv"))

    concat_metrics(perf_files, os.path.join(summary_dir, "raw_dapars2_parallel_metrics.tsv"))
    gt_rows = merge_failures(gt_fail_files)
    nondapars_rows = merge_failures(nondapars_fail_files)
    dapars_rows = merge_failures(dapars_fail_files)
    eval_rows = merge_failures(eval_fail_files)
    all_rows = gt_rows + nondapars_rows + dapars_rows + eval_rows

    write_tsv(os.path.join(summary_dir, "raw_dapars2_parallel_gt_failures.tsv"), gt_rows, FAIL_HEADERS)
    write_tsv(os.path.join(summary_dir, "raw_dapars2_parallel_nondapars_failures.tsv"), nondapars_rows, FAIL_HEADERS)
    write_tsv(os.path.join(summary_dir, "raw_dapars2_parallel_dapars_failures.tsv"), dapars_rows, FAIL_HEADERS)
    write_tsv(os.path.join(summary_dir, "raw_dapars2_parallel_missing_failed.tsv"), all_rows, FAIL_HEADERS)

    os.makedirs(os.path.dirname(args.done_file), exist_ok=True)
    with open(args.done_file, "w", encoding="utf-8") as fh:
        fh.write("ok\n")


if __name__ == "__main__":
    main()
