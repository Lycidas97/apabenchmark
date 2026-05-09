#!/usr/bin/env python3
import argparse
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from itertools import combinations


def clean_cell_type(cell_type: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", cell_type).strip("_")


def ensure_parent(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def read_tsv(path: str):
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: str, rows, fieldnames):
    ensure_parent(path)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_cmd(cmd, log_path=None, env=None):
    if log_path:
        ensure_parent(log_path)
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=env,
        check=False,
    )
    if log_path:
        with open(log_path, "w", encoding="utf-8") as lf:
            lf.write("COMMAND:\n")
            lf.write(" ".join(cmd) + "\n\n")
            lf.write("STDOUT+STDERR:\n")
            lf.write(proc.stdout or "")
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}")
    return proc.stdout


def concat_tsv(files, out_path, header_fields=None):
    files = [f for f in files if os.path.exists(f) and os.path.getsize(f) > 0]
    if not files:
        ensure_parent(out_path)
        with open(out_path, "w", encoding="utf-8", newline="") as out_fh:
            if header_fields:
                out_fh.write("\t".join(header_fields) + "\n")
        return False
    ensure_parent(out_path)
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


def select_top_samples(sample_top3_rows, top_n):
    totals = defaultdict(int)
    for row in sample_top3_rows:
        totals[row["sample"]] += int(row["cell_count"])
    ranked = sorted(totals.items(), key=lambda x: (-x[1], x[0]))
    return [s for s, _ in ranked[:top_n]], totals


def load_sample_celltypes(gt_mtx_path):
    barcode_to_type = {}
    with open(gt_mtx_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if "cell_type" not in header:
            raise RuntimeError(f"cell_type column not found: {gt_mtx_path}")
        ct_idx = header.index("cell_type")
        for row in reader:
            if not row:
                continue
            bc = row[0]
            if not bc:
                continue
            if ct_idx >= len(row):
                continue
            barcode_to_type[bc] = row[ct_idx]
    return barcode_to_type


def infer_barcode_transform(tool_mtx_path, gt_barcode_set):
    candidates = [
        ("identity", lambda x: x),
        ("add_-1", lambda x: x if x.endswith("-1") else f"{x}-1"),
        ("strip_-1", lambda x: x[:-2] if x.endswith("-1") else x),
        ("strip_suffix", lambda x: re.sub(r"-\d+$", "", x)),
        ("strip_suffix_add_-1", lambda x: re.sub(r"-\d+$", "", x) + "-1"),
    ]

    tool_barcodes = []
    with open(tool_mtx_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)  # header row (PAS ids)
        for row in reader:
            if not row:
                continue
            tool_barcodes.append(row[0])

    best_name = None
    best_fn = None
    best_hits = -1
    for name, fn in candidates:
        hits = sum(1 for bc in tool_barcodes if fn(bc) in gt_barcode_set)
        if hits > best_hits:
            best_name = name
            best_fn = fn
            best_hits = hits

    return best_name, best_fn, best_hits, len(tool_barcodes)


def subset_gt_expr_and_filter_pas(gt_mtx_path, selected_cells, out_expr_path):
    with open(gt_mtx_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        ct_idx = header.index("cell_type")
        pas_start = ct_idx + 1
        pas_names = header[pas_start:]
        seen_nonzero = [False] * len(pas_names)
        selected_count = 0
        for row in reader:
            if not row:
                continue
            bc = row[0]
            if bc not in selected_cells:
                continue
            selected_count += 1
            vals = row[pas_start:pas_start + len(pas_names)]
            for i, v in enumerate(vals):
                if not seen_nonzero[i] and v and v != "0":
                    seen_nonzero[i] = True

    keep_idx = [i for i, k in enumerate(seen_nonzero) if k]
    keep_pas = [pas_names[i] for i in keep_idx]

    ensure_parent(out_expr_path)
    with open(gt_mtx_path, "r", encoding="utf-8", newline="") as in_fh, open(
        out_expr_path, "w", encoding="utf-8", newline=""
    ) as out_fh:
        reader = csv.reader(in_fh, delimiter="\t")
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        header = next(reader)
        ct_idx = header.index("cell_type")
        pas_start = ct_idx + 1
        writer.writerow(["barcode", "cell_type"] + [header[pas_start + i] for i in keep_idx])
        for row in reader:
            if not row:
                continue
            bc = row[0]
            if bc not in selected_cells:
                continue
            ct = row[ct_idx]
            vals = row[pas_start:pas_start + len(header[pas_start:])]
            out_vals = [vals[i] if i < len(vals) else "0" for i in keep_idx]
            writer.writerow([bc, ct] + out_vals)

    return selected_count, len(pas_names), len(keep_pas), set(keep_pas)


def filter_gt_pas_bed(gt_pas_path, keep_pas_set, out_pas_path):
    ensure_parent(out_pas_path)
    kept = 0
    with open(gt_pas_path, "r", encoding="utf-8", newline="") as in_fh, open(
        out_pas_path, "w", encoding="utf-8", newline=""
    ) as out_fh:
        reader = csv.reader(in_fh, delimiter="\t")
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        header_written = False
        for row in reader:
            if not row:
                continue
            if not header_written and row[0] == "chr" and len(row) >= 8:
                writer.writerow(row[:8])
                header_written = True
                continue
            if len(row) < 8:
                continue
            pas_id = f"{row[7]}|{row[0]}:{row[1]}:{row[2]}:{row[5]}"
            if pas_id in keep_pas_set or row[3] in keep_pas_set:
                writer.writerow(row[:8])
                kept += 1
    return kept


def subset_tool_expr_and_filter_pas(tool_mtx_path, selected_cells, barcode_transform, out_expr_path):
    with open(tool_mtx_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if len(header) < 2:
            raise RuntimeError(f"Invalid tool matrix header (need barcode + PAS columns): {tool_mtx_path}")
        pas_names = header[1:]
        seen_nonzero = [False] * len(pas_names)
        selected_count = 0
        for row in reader:
            if not row:
                continue
            mapped = barcode_transform(row[0])
            if mapped not in selected_cells:
                continue
            selected_count += 1
            vals = row[1:1 + len(pas_names)]
            for i, v in enumerate(vals):
                if not seen_nonzero[i] and v and v != "0":
                    seen_nonzero[i] = True

    keep_idx = [i for i, k in enumerate(seen_nonzero) if k]
    keep_pas = [pas_names[i] for i in keep_idx]

    ensure_parent(out_expr_path)
    with open(tool_mtx_path, "r", encoding="utf-8", newline="") as in_fh, open(
        out_expr_path, "w", encoding="utf-8", newline=""
    ) as out_fh:
        reader = csv.reader(in_fh, delimiter="\t")
        header = next(reader)
        if len(header) < 2:
            raise RuntimeError(f"Invalid tool matrix header (need barcode + PAS columns): {tool_mtx_path}")
        pas_names = header[1:]
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["barcode"] + [pas_names[i] for i in keep_idx])
        for row in reader:
            if not row:
                continue
            mapped = barcode_transform(row[0])
            if mapped not in selected_cells:
                continue
            vals = row[1:1 + len(pas_names)]
            out_vals = [vals[i] if i < len(vals) else "0" for i in keep_idx]
            writer.writerow([mapped] + out_vals)

    return selected_count, len(pas_names), len(keep_pas), set(keep_pas)


def filter_tool_pas_bed(tool_pas_path, keep_pas_set, out_pas_path):
    ensure_parent(out_pas_path)
    kept = 0
    with open(tool_pas_path, "r", encoding="utf-8", newline="") as in_fh, open(
        out_pas_path, "w", encoding="utf-8", newline=""
    ) as out_fh:
        reader = csv.reader(in_fh, delimiter="\t")
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        for row in reader:
            if not row or len(row) < 4:
                continue
            if row[3] in keep_pas_set:
                writer.writerow(row)
                kept += 1
    return kept


def parse_cluster_label(pdui_col):
    base = os.path.basename(pdui_col)
    if base.startswith("cluster"):
        label = base[len("cluster"):]
        label = re.sub(r"\.(?:trans|intersect)_PDUI$", "", label)
        label = re.sub(r"_PDUI$", "", label)
        return label
    return base


def score_triplet_match(label, pdui_col, target_raw, target_clean):
    label_clean = clean_cell_type(label)
    score = 0
    if label == target_raw:
        score = max(score, 100)
    if label_clean == target_clean:
        score = max(score, 95)
    if target_raw in pdui_col:
        score = max(score, 90)
    if target_clean and target_clean in clean_cell_type(pdui_col):
        score = max(score, 85)

    m = re.match(r"^(\d+)-", target_raw)
    if m:
        num = int(m.group(1))
        if re.search(rf"CLAS_{num:02d}(?:\D|$)", label) or re.search(rf"CLAS_{num:02d}(?:\D|$)", pdui_col):
            score = max(score, 92)
        if re.search(rf"CLAS_{num}(?:\D|$)", label) or re.search(rf"CLAS_{num}(?:\D|$)", pdui_col):
            score = max(score, 88)

    return score


def extract_dapars_pair_input(src_path, group1_raw, group2_raw, group1_clean, group2_clean, out_path):
    ensure_parent(out_path)
    with open(src_path, "r", encoding="utf-8", newline="") as in_fh, open(
        out_path, "w", encoding="utf-8", newline=""
    ) as out_fh:
        reader = csv.reader(in_fh, delimiter="\t")
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")

        header = next(reader)
        if len(header) < 10:
            raise RuntimeError(f"Unexpected DaPars2 header format: {src_path}")

        triplets = []
        for i in range(4, len(header), 3):
            if i + 2 >= len(header):
                break
            pdui_col = header[i + 2]
            triplets.append(
                {
                    "start": i,
                    "long": header[i],
                    "short": header[i + 1],
                    "pdui": pdui_col,
                    "label": parse_cluster_label(pdui_col),
                }
            )

        def pick_triplet(target_raw, target_clean):
            scored = []
            for t in triplets:
                s = score_triplet_match(t["label"], t["pdui"], target_raw, target_clean)
                if s > 0:
                    scored.append((s, t))
            if not scored:
                raise RuntimeError(
                    f"No DaPars2 group match for {target_raw} ({target_clean}) in {src_path}"
                )
            scored.sort(key=lambda x: (-x[0], x[1]["start"]))
            return scored[0][1]

        t1 = pick_triplet(group1_raw, group1_clean)
        t2 = pick_triplet(group2_raw, group2_clean)

        select_idx = [0, 1, 2, 3, t1["start"], t1["start"] + 1, t1["start"] + 2, t2["start"], t2["start"] + 1, t2["start"] + 2]
        out_header = [
            "Gene",
            "fit_value",
            "Predicted_Proximal_APA",
            "Loci",
            f"{group1_clean}_long_exp",
            f"{group1_clean}_short_exp",
            f"{group1_clean}_PDUI",
            f"{group2_clean}_long_exp",
            f"{group2_clean}_short_exp",
            f"{group2_clean}_PDUI",
        ]
        writer.writerow(out_header)

        row_count = 0
        for row in reader:
            if not row:
                continue
            out_row = [row[i] if i < len(row) else "" for i in select_idx]
            writer.writerow(out_row)
            row_count += 1

    return row_count


def find_existing_pair_rows(pair_rows, selected_samples):
    selected = set(selected_samples)
    rows = [r for r in pair_rows if r["sample"] in selected]
    rows.sort(key=lambda r: (r["sample"], int(r["pair_rank"])))
    return rows


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_analysis_root = os.path.abspath(os.path.join(script_dir, "..", "05_benchmark_analysis"))

    parser = argparse.ArgumentParser(description="Run raw benchmark (traditional + dapars) on auto-selected top samples")
    parser.add_argument("--project-root", default=os.environ.get("APABENCHMARK_FINAL_ROOT", os.getcwd()))
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--split-num", type=int, default=5)
    parser.add_argument(
        "--pair-rank-max",
        type=int,
        default=1,
        help="Maximum pair_rank per sample to run (default: 1 for quick multi-sample validation)",
    )
    parser.add_argument("--conda-exe", default=os.environ.get("APABENCHMARK_CONDA_EXE", "conda"))
    parser.add_argument("--dexseq-env", default="dexseq")
    parser.add_argument("--scmapa-env", default="scmapa_r")
    parser.add_argument(
        "--dapars-script-root",
        default=default_analysis_root,
        help="Path containing scripts/10_stage_pas_match_quantify.R, scripts/20_stage_apa_and_te.R, scripts/30_stage_export_results.R, scmapa/*.R, and scripts/calculate_dapars2_de_apa_performance.R",
    )
    parser.add_argument(
        "--output-root",
        default=None,
        help="Output run root (default: <project>/data/result/benchmark/analysis_smoke/top{n}_split{split}_{timestamp})",
    )
    args = parser.parse_args()

    project_root = os.path.abspath(args.project_root)
    data_root = os.path.join(project_root, "data")
    raw_bam_root = os.path.join(data_root, "result", "benchmark", "raw_bam")
    raw_gt_root = os.path.join(data_root, "raw_data", "raw_bam")
    manifest_root = os.path.join(data_root, "result", "benchmark", "analysis", "raw_pair_manifest")

    stage10_perf_script = os.path.join(args.dapars_script_root, "scripts", "10_stage_pas_match_quantify.R")
    stage20_perf_script = os.path.join(args.dapars_script_root, "scripts", "20_stage_apa_and_te.R")
    stage30_perf_script = os.path.join(args.dapars_script_root, "scripts", "30_stage_export_results.R")
    multi_group_converter = os.path.join(args.dapars_script_root, "scmapa", "multi_group_PAS_to_DaPars2_converter.R")
    integrated_converter = os.path.join(args.dapars_script_root, "scmapa", "integrated_pas_to_dapars2_converter.R")
    scmapa_script = os.path.join(args.dapars_script_root, "scmapa", "parameterized_scmapa_analysis.R")
    dapars_perf_script = os.path.join(args.dapars_script_root, "scripts", "calculate_dapars2_de_apa_performance.R")

    for p in [
        stage10_perf_script,
        stage20_perf_script,
        stage30_perf_script,
        multi_group_converter,
        integrated_converter,
        scmapa_script,
        dapars_perf_script,
    ]:
        if not os.path.exists(p):
            raise RuntimeError(f"Required script not found: {p}")

    sample_top3_path = os.path.join(manifest_root, "sample_top3_celltypes.tsv")
    pair_manifest_path = os.path.join(manifest_root, "sample_pair_manifest.tsv")
    sample_top3_rows = read_tsv(sample_top3_path)
    pair_rows_all = read_tsv(pair_manifest_path)

    selected_samples, sample_totals = select_top_samples(sample_top3_rows, args.top_n)
    pair_rows = find_existing_pair_rows(pair_rows_all, selected_samples)
    pair_rows = [r for r in pair_rows if int(r["pair_rank"]) <= args.pair_rank_max]

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = args.output_root or os.path.join(
        data_root,
        "result",
        "benchmark",
        "analysis_smoke",
        f"top{args.top_n}_split{args.split_num}_{ts}",
    )

    raw_pair_manifest_dir = os.path.join(run_root, "raw_pair_manifest")
    raw_pair_inputs_dir = os.path.join(run_root, "raw_pair_inputs")
    raw_qc_dir = os.path.join(run_root, "raw_qc")
    raw_perf_dir = os.path.join(run_root, "raw_performance")
    raw_perf_out_dir = os.path.join(raw_perf_dir, "by_sample_pair_tool")
    raw_perf_log_dir = os.path.join(raw_perf_dir, "logs")
    raw_perf_summary_dir = os.path.join(raw_perf_dir, "summary")

    raw_dapars_dir = os.path.join(run_root, "raw_dapars2_format")
    dapars_gt_dir = os.path.join(raw_dapars_dir, "ground_truth")
    dapars_conv_dir = os.path.join(raw_dapars_dir, "converted_inputs")
    dapars_de_dir = os.path.join(raw_dapars_dir, "differential_apa")
    dapars_perf_dir = os.path.join(raw_dapars_dir, "performance")
    dapars_log_dir = os.path.join(raw_dapars_dir, "logs")
    dapars_summary_dir = os.path.join(raw_dapars_dir, "summary")

    for d in [
        raw_pair_manifest_dir,
        raw_pair_inputs_dir,
        raw_qc_dir,
        raw_perf_out_dir,
        raw_perf_log_dir,
        raw_perf_summary_dir,
        dapars_gt_dir,
        dapars_conv_dir,
        dapars_de_dir,
        dapars_perf_dir,
        dapars_log_dir,
        dapars_summary_dir,
    ]:
        os.makedirs(d, exist_ok=True)

    selected_top3_rows = [r for r in sample_top3_rows if r["sample"] in set(selected_samples)]
    write_tsv(
        os.path.join(raw_pair_manifest_dir, "sample_top3_celltypes.tsv"),
        selected_top3_rows,
        ["sample", "top_rank", "cell_type", "cell_count", "cell_type_clean"],
    )
    write_tsv(
        os.path.join(raw_pair_manifest_dir, "sample_pair_manifest.tsv"),
        pair_rows,
        ["sample", "pair_rank", "pair_label", "group1", "group2", "group1_clean", "group2_clean", "group1_count", "group2_count"],
    )
    prepare_meta_rows = [
        {
            "run_root": run_root,
            "selected_sample_count": str(len(selected_samples)),
            "selected_pair_count": str(len(pair_rows)),
            "top_n": str(args.top_n),
            "split_num": str(args.split_num),
            "pair_rank_max": str(args.pair_rank_max),
            "selected_samples": ",".join(selected_samples),
        }
    ]
    write_tsv(
        os.path.join(raw_pair_manifest_dir, "prepare_metadata.tsv"),
        prepare_meta_rows,
        ["run_root", "selected_sample_count", "selected_pair_count", "top_n", "split_num", "pair_rank_max", "selected_samples"],
    )

    print(f"Run root: {run_root}")
    print(f"Selected samples ({len(selected_samples)}): {', '.join(selected_samples)}")

    qc_rows = []
    perf_failures = []
    dapars_failures = []

    per_sample_celltypes = {}
    per_sample_transform = {}

    for sample in selected_samples:
        gt_mtx = os.path.join(raw_gt_root, sample, "bam.expr.tsv")
        gt_pas = os.path.join(raw_gt_root, sample, "pas.bed")
        sierra_mtx = os.path.join(raw_bam_root, "sierra", sample, "pas_counts.tsv")
        sierra_pas = os.path.join(raw_bam_root, "sierra", sample, "pas.bed")
        scmapa_src = os.path.join(raw_bam_root, "scmapa", sample, "DaPars2_Result_all_chrs.txt")
        dapars2_src = os.path.join(raw_bam_root, "dapars2", sample, "DaPars2_Result_all_chrs.txt")

        required = [gt_mtx, gt_pas, sierra_mtx, sierra_pas, scmapa_src, dapars2_src]
        missing = [p for p in required if not os.path.exists(p)]
        if missing:
            for m in missing:
                msg = f"missing_input:{m}"
                perf_failures.append({"sample": sample, "pair_label": "ALL", "tool": "sierra", "stage": "input_check", "message": msg})
                dapars_failures.append({"sample": sample, "pair_label": "ALL", "tool": "all", "stage": "input_check", "message": msg})
            continue

        if sample not in per_sample_celltypes:
            per_sample_celltypes[sample] = load_sample_celltypes(gt_mtx)

        barcode_to_type = per_sample_celltypes[sample]
        all_gt_barcodes = set(barcode_to_type.keys())

        if sample not in per_sample_transform:
            t_name, t_fn, hits, total = infer_barcode_transform(sierra_mtx, all_gt_barcodes)
            per_sample_transform[sample] = (t_name, t_fn)
            print(f"[{sample}] sierra barcode transform: {t_name} ({hits}/{total} matched)")

        transform_name, transform_fn = per_sample_transform[sample]

        sample_pairs = [r for r in pair_rows if r["sample"] == sample]
        for pair in sample_pairs:
            pair_label = pair["pair_label"]
            group1 = pair["group1"]
            group2 = pair["group2"]
            g1_clean = pair["group1_clean"]
            g2_clean = pair["group2_clean"]
            sample_pair_id = f"{sample}__{pair_label}"

            selected_cells = {bc for bc, ct in barcode_to_type.items() if ct == group1 or ct == group2}
            if not selected_cells:
                msg = "no_cells_for_pair"
                perf_failures.append({"sample": sample, "pair_label": pair_label, "tool": "sierra", "stage": "prepare_pair", "message": msg})
                dapars_failures.append({"sample": sample, "pair_label": pair_label, "tool": "all", "stage": "prepare_pair", "message": msg})
                continue

            gt_pair_expr = os.path.join(raw_pair_inputs_dir, "ground_truth", sample, pair_label, "gt_expr_pair.tsv")
            gt_pair_pas = os.path.join(raw_pair_inputs_dir, "ground_truth", sample, pair_label, "gt_pas_pair.bed")
            pd_pair_expr = os.path.join(raw_pair_inputs_dir, "sierra", sample, pair_label, "pd_pas_counts.tsv")
            pd_pair_pas = os.path.join(raw_pair_inputs_dir, "sierra", sample, pair_label, "pd_pas.bed")

            try:
                gt_selected_n, gt_pas_total, gt_pas_kept, gt_keep_set = subset_gt_expr_and_filter_pas(
                    gt_mtx,
                    selected_cells,
                    gt_pair_expr,
                )
                gt_bed_kept = filter_gt_pas_bed(gt_pas, gt_keep_set, gt_pair_pas)
                qc_rows.append(
                    {
                        "sample": sample,
                        "pair_label": pair_label,
                        "tool": "ground_truth",
                        "selected_cells": str(gt_selected_n),
                        "pas_total": str(gt_pas_total),
                        "pas_kept_nonzero": str(gt_pas_kept),
                        "pas_dropped_zero": str(gt_pas_total - gt_pas_kept),
                        "bed_rows_kept": str(gt_bed_kept),
                        "barcode_transform": "identity",
                    }
                )
            except Exception as e:
                msg = f"gt_prepare_failed:{e}"
                perf_failures.append({"sample": sample, "pair_label": pair_label, "tool": "ground_truth", "stage": "prepare_pair", "message": msg})
                dapars_failures.append({"sample": sample, "pair_label": pair_label, "tool": "ground_truth", "stage": "prepare_pair", "message": msg})
                continue

            try:
                pd_selected_n, pd_pas_total, pd_pas_kept, pd_keep_set = subset_tool_expr_and_filter_pas(
                    sierra_mtx,
                    selected_cells,
                    transform_fn,
                    pd_pair_expr,
                )
                pd_bed_kept = filter_tool_pas_bed(sierra_pas, pd_keep_set, pd_pair_pas)
                qc_rows.append(
                    {
                        "sample": sample,
                        "pair_label": pair_label,
                        "tool": "sierra",
                        "selected_cells": str(pd_selected_n),
                        "pas_total": str(pd_pas_total),
                        "pas_kept_nonzero": str(pd_pas_kept),
                        "pas_dropped_zero": str(pd_pas_total - pd_pas_kept),
                        "bed_rows_kept": str(pd_bed_kept),
                        "barcode_transform": transform_name,
                    }
                )
            except Exception as e:
                msg = f"sierra_prepare_failed:{e}"
                perf_failures.append({"sample": sample, "pair_label": pair_label, "tool": "sierra", "stage": "prepare_pair", "message": msg})
                dapars_failures.append({"sample": sample, "pair_label": pair_label, "tool": "sierra", "stage": "prepare_pair", "message": msg})
                continue

            # Traditional performance (sierra only)
            try:
                perf_log_stage10 = os.path.join(raw_perf_log_dir, "sierra", sample, pair_label, "run_performance_stage10.log")
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
                        "sierra",
                        "--sample",
                        sample_pair_id,
                        "--output_dir",
                        raw_perf_out_dir,
                        "--dexseq_split_num",
                        str(args.split_num),
                    ],
                    log_path=perf_log_stage10,
                )

                stage10_bundle = os.path.join(
                    raw_perf_out_dir,
                    "sierra",
                    f"{sample_pair_id}_stage10_bundle.rds",
                )
                if not os.path.exists(stage10_bundle):
                    raise RuntimeError(f"missing_stage10_bundle:{stage10_bundle}")

                perf_log_stage20 = os.path.join(raw_perf_log_dir, "sierra", sample, pair_label, "run_performance_stage20.log")
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
                    log_path=perf_log_stage20,
                )

                stage20_bundle = os.path.join(
                    raw_perf_out_dir,
                    "sierra",
                    f"{sample_pair_id}_stage20_bundle.rds",
                )
                if not os.path.exists(stage20_bundle):
                    raise RuntimeError(f"missing_stage20_bundle:{stage20_bundle}")

                perf_log_stage30 = os.path.join(raw_perf_log_dir, "sierra", sample, pair_label, "run_performance_stage30.log")
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
                    log_path=perf_log_stage30,
                )
            except Exception as e:
                perf_failures.append({"sample": sample, "pair_label": pair_label, "tool": "sierra", "stage": "calculate_benchmark_performance", "message": str(e)})

            # Dapars line
            try:
                # 1) GT convert and scMAPA
                gt_dapars = os.path.join(dapars_gt_dir, sample, pair_label, "dapars2_format.txt")
                ensure_parent(gt_dapars)
                gt_conv_log = os.path.join(dapars_log_dir, "ground_truth", sample, pair_label, "convert.log")
                run_cmd(
                    [
                        args.conda_exe,
                        "run",
                        "-n",
                        args.scmapa_env,
                        "Rscript",
                        multi_group_converter,
                        "-i",
                        gt_pair_expr,
                        "-o",
                        gt_dapars,
                        "--min_exp",
                        "0",
                        "--min_cells",
                        "0",
                    ],
                    log_path=gt_conv_log,
                )

                gt_de_out = os.path.join(dapars_de_dir, "ground_truth", sample, pair_label)
                os.makedirs(gt_de_out, exist_ok=True)
                gt_scmapa_log = os.path.join(dapars_log_dir, "ground_truth", sample, pair_label, "scmapa.log")
                run_cmd(
                    [
                        args.conda_exe,
                        "run",
                        "-n",
                        args.scmapa_env,
                        "Rscript",
                        scmapa_script,
                        "-i",
                        gt_dapars,
                        "-o",
                        gt_de_out,
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
                    log_path=gt_scmapa_log,
                )
                gt_apa = os.path.join(gt_de_out, f"{g1_clean}_vs_{g2_clean}", "APA_results.txt")

                # 2) Sierra convert and scMAPA
                sierra_dapars = os.path.join(dapars_conv_dir, "sierra", sample, pair_label, "dapars2_format.txt")
                ensure_parent(sierra_dapars)
                sierra_conv_log = os.path.join(dapars_log_dir, "sierra", sample, pair_label, "convert.log")
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
                        sierra_dapars,
                        "--min_exp",
                        "0",
                        "--min_cells",
                        "0",
                    ],
                    log_path=sierra_conv_log,
                )

                sierra_de_out = os.path.join(dapars_de_dir, "sierra", sample, pair_label)
                os.makedirs(sierra_de_out, exist_ok=True)
                sierra_scmapa_log = os.path.join(dapars_log_dir, "sierra", sample, pair_label, "scmapa.log")
                run_cmd(
                    [
                        args.conda_exe,
                        "run",
                        "-n",
                        args.scmapa_env,
                        "Rscript",
                        scmapa_script,
                        "-i",
                        sierra_dapars,
                        "-o",
                        sierra_de_out,
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
                    log_path=sierra_scmapa_log,
                )

                for tool_name, src_file in [("scmapa", scmapa_src), ("dapars2", dapars2_src)]:
                    pair_dapars = os.path.join(dapars_conv_dir, tool_name, sample, pair_label, "dapars2_format.txt")
                    extract_dapars_pair_input(src_file, group1, group2, g1_clean, g2_clean, pair_dapars)

                    tool_de_out = os.path.join(dapars_de_dir, tool_name, sample, pair_label)
                    os.makedirs(tool_de_out, exist_ok=True)
                    tool_scmapa_log = os.path.join(dapars_log_dir, tool_name, sample, pair_label, "scmapa.log")
                    run_cmd(
                        [
                            args.conda_exe,
                            "run",
                            "-n",
                            args.scmapa_env,
                            "Rscript",
                            scmapa_script,
                            "-i",
                            pair_dapars,
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
                        log_path=tool_scmapa_log,
                    )

                # 3) Dapars performance for sierra/scmapa/dapars2
                for tool_name in ["sierra", "scmapa", "dapars2"]:
                    pd_apa = os.path.join(
                        dapars_de_dir,
                        tool_name,
                        sample,
                        pair_label,
                        f"{g1_clean}_vs_{g2_clean}",
                        "APA_results.txt",
                    )
                    perf_out = os.path.join(
                        dapars_perf_dir,
                        tool_name,
                        sample,
                        pair_label,
                        "dapars2_performance.tsv",
                    )
                    ensure_parent(perf_out)
                    perf_log = os.path.join(dapars_log_dir, tool_name, sample, pair_label, "performance.log")
                    run_cmd(
                        [
                            args.conda_exe,
                            "run",
                            "-n",
                            args.dexseq_env,
                            "Rscript",
                            dapars_perf_script,
                            "--pd_apa",
                            pd_apa,
                            "--gt_apa",
                            gt_apa,
                            "--tool",
                            tool_name,
                            "--sample",
                            sample_pair_id,
                            "--match_type",
                            f"{g1_clean}_vs_{g2_clean}",
                            "--filter_type_1",
                            "scmapa_adPval_0.05",
                            "--filter_type_2",
                            "scmapa_or_0.1",
                            "--output",
                            perf_out,
                        ],
                        log_path=perf_log,
                    )
            except Exception as e:
                dapars_failures.append({"sample": sample, "pair_label": pair_label, "tool": "all", "stage": "dapars_pipeline", "message": str(e)})

    # QC summary
    write_tsv(
        os.path.join(raw_qc_dir, "pas_filter_overview.tsv"),
        qc_rows,
        [
            "sample",
            "pair_label",
            "tool",
            "selected_cells",
            "pas_total",
            "pas_kept_nonzero",
            "pas_dropped_zero",
            "bed_rows_kept",
            "barcode_transform",
        ],
    )

    # Traditional summaries (sierra only in this run)
    match_files = []
    quantify_files = []
    deapa_files = []
    for sample in selected_samples:
        sample_pairs = [r for r in pair_rows if r["sample"] == sample]
        for pair in sample_pairs:
            pair_label = pair["pair_label"]
            sample_pair_id = f"{sample}__{pair_label}"
            base = os.path.join(raw_perf_out_dir, "sierra", f"{sample_pair_id}")
            match_files.append(f"{base}_match_performance.tsv")
            quantify_files.append(f"{base}_pas_quantify_performance.tsv")
            deapa_files.append(f"{base}_de_apa_performance.tsv")

    concat_tsv(
        match_files,
        os.path.join(raw_perf_summary_dir, "raw_performance_match_metrics.tsv"),
        header_fields=["sample", "tool", "protocol", "match_type", "tp", "fp", "fn", "precision", "recall", "f1"],
    )
    concat_tsv(
        quantify_files,
        os.path.join(raw_perf_summary_dir, "raw_performance_pas_quantify_metrics.tsv"),
        header_fields=[
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
        ],
    )
    concat_tsv(
        deapa_files,
        os.path.join(raw_perf_summary_dir, "raw_performance_de_apa_metrics.tsv"),
        header_fields=[
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
        ],
    )

    write_tsv(
        os.path.join(raw_perf_summary_dir, "raw_performance_missing_failed.tsv"),
        perf_failures,
        ["sample", "pair_label", "tool", "stage", "message"],
    )

    # Dapars summary
    dapars_metric_files = []
    for tool_name in ["sierra", "scmapa", "dapars2"]:
        for sample in selected_samples:
            sample_pairs = [r for r in pair_rows if r["sample"] == sample]
            for pair in sample_pairs:
                pair_label = pair["pair_label"]
                dapars_metric_files.append(
                    os.path.join(
                        dapars_perf_dir,
                        tool_name,
                        sample,
                        pair_label,
                        "dapars2_performance.tsv",
                    )
                )

    concat_tsv(
        dapars_metric_files,
        os.path.join(dapars_summary_dir, "raw_dapars2_metrics.tsv"),
        header_fields=[
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
        ],
    )
    write_tsv(
        os.path.join(dapars_summary_dir, "raw_dapars2_missing_failed.tsv"),
        dapars_failures,
        ["sample", "pair_label", "tool", "stage", "message"],
    )

    print("Done.")
    print(f"Run root: {run_root}")
    print(f"Traditional failures: {len(perf_failures)}")
    print(f"DaPars failures: {len(dapars_failures)}")


if __name__ == "__main__":
    main()
