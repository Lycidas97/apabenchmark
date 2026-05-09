#!/usr/bin/env python3
import csv
import os
import re
import subprocess
from typing import List, Dict, Tuple


def ensure_parent(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def read_tsv(path: str) -> List[Dict[str, str]]:
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: str, rows: List[Dict[str, str]], fieldnames: List[str]) -> None:
    ensure_parent(path)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_cmd(cmd: List[str], log_path: str = None) -> None:
    if log_path:
        ensure_parent(log_path)
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
    if log_path:
        with open(log_path, "w", encoding="utf-8") as lf:
            lf.write("COMMAND:\n")
            lf.write(" ".join(cmd) + "\n\n")
            lf.write("STDOUT+STDERR:\n")
            lf.write(proc.stdout or "")
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}")


def touch_done(path: str) -> None:
    ensure_parent(path)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("ok\n")


def parse_cluster_label(pdui_col: str) -> str:
    base = os.path.basename(pdui_col)
    if base.startswith("cluster"):
        label = base[len("cluster"):]
        label = re.sub(r"\.(?:trans|intersect)_PDUI$", "", label)
        label = re.sub(r"_PDUI$", "", label)
        return label
    return base


def clean_cell_type(cell_type: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", cell_type).strip("_")


def score_triplet_match(label: str, pdui_col: str, target_raw: str, target_clean: str) -> int:
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


def extract_dapars_pair_input(src_path: str, group1_raw: str, group2_raw: str, group1_clean: str, group2_clean: str, out_path: str) -> int:
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
            triplets.append({"start": i, "pdui": pdui_col, "label": parse_cluster_label(pdui_col)})

        def pick_triplet(target_raw: str, target_clean: str):
            scored: List[Tuple[int, Dict[str, str]]] = []
            for t in triplets:
                s = score_triplet_match(t["label"], t["pdui"], target_raw, target_clean)
                if s > 0:
                    scored.append((s, t))
            if not scored:
                raise RuntimeError(f"No DaPars2 group match for {target_raw} ({target_clean}) in {src_path}")
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
