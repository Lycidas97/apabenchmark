#!/usr/bin/env python3
import argparse
import os
from raw_parallel_common import read_tsv, write_tsv


MATCH_HEADERS = ["sample", "tool", "protocol", "match_type", "tp", "fp", "fn", "precision", "recall", "f1"]
QUANT_HEADERS = [
    "sample",
    "tool",
    "protocol",
    "match_type",
    "cor_pas",
    "rmse_pas",
    "mae_pas",
    "mape_pas",
    "rmse_pas_ct",
    "mae_pas_ct",
    "mape_pas_ct",
]
DEAPA_HEADERS = [
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


def concat_tsv(files, out_path, header_fields):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as out_fh:
        out_fh.write("\t".join(header_fields) + "\n")
        for path in files:
            if not os.path.exists(path):
                continue
            with open(path, "r", encoding="utf-8", newline="") as in_fh:
                lines = in_fh.readlines()
            if len(lines) <= 1:
                continue
            for line in lines[1:]:
                out_fh.write(line)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run-root", required=True)
    p.add_argument("--done-file", required=True)
    args = p.parse_args()

    pair_manifest = os.path.join(args.run_root, "raw_pair_manifest", "sample_pair_manifest.tsv")
    pair_rows = read_tsv(pair_manifest)

    perf_dir = os.path.join(args.run_root, "raw_performance")
    by_pair_root = os.path.join(perf_dir, "by_sample_pair_tool")
    fail_root = os.path.join(perf_dir, "failures", "eval")
    summary_dir = os.path.join(perf_dir, "summary")
    os.makedirs(summary_dir, exist_ok=True)

    tools = []
    if os.path.isdir(by_pair_root):
        tools = sorted([name for name in os.listdir(by_pair_root) if os.path.isdir(os.path.join(by_pair_root, name))])

    match_files = []
    quantify_files = []
    deapa_files = []
    fail_files = []
    for tool in tools:
        by_pair_dir = os.path.join(by_pair_root, tool)
        fail_dir = os.path.join(fail_root, tool)
        for pair in pair_rows:
            sample = pair["sample"]
            pair_label = pair["pair_label"]
            sample_pair_id = f"{sample}__{pair_label}"
            match_files.append(os.path.join(by_pair_dir, f"{sample_pair_id}_match_performance.tsv"))
            quantify_files.append(os.path.join(by_pair_dir, f"{sample_pair_id}_pas_quantify_performance.tsv"))
            deapa_files.append(os.path.join(by_pair_dir, f"{sample_pair_id}_de_apa_performance.tsv"))
            fail_files.append(os.path.join(fail_dir, sample, pair_label, "failure.tsv"))

    concat_tsv(match_files, os.path.join(summary_dir, "raw_performance_match_metrics.tsv"), MATCH_HEADERS)
    concat_tsv(quantify_files, os.path.join(summary_dir, "raw_performance_pas_quantify_metrics.tsv"), QUANT_HEADERS)
    concat_tsv(deapa_files, os.path.join(summary_dir, "raw_performance_de_apa_metrics.tsv"), DEAPA_HEADERS)

    failures = []
    for ff in fail_files:
        if os.path.exists(ff):
            failures.extend(read_tsv(ff))
    write_tsv(os.path.join(summary_dir, "raw_performance_missing_failed.tsv"), failures, FAIL_HEADERS)

    os.makedirs(os.path.dirname(args.done_file), exist_ok=True)
    with open(args.done_file, "w", encoding="utf-8") as fh:
        fh.write("ok\n")


if __name__ == "__main__":
    main()
