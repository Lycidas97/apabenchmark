import argparse
import glob
import os


def has_required_groups(header_columns, group1, group2):
    required = {
        f"{group1}_long_exp",
        f"{group1}_short_exp",
        f"{group2}_long_exp",
        f"{group2}_short_exp",
    }
    return required.issubset(set(header_columns))


def read_header_columns(file_path):
    with open(file_path, "r", encoding="utf-8", errors="ignore") as fh:
        first_line = fh.readline().rstrip("\n")
    if not first_line:
        return []
    return first_line.split("\t")


def main():
    parser = argparse.ArgumentParser(
        description="Generate exclude list for DaPars2 single-group samples."
    )
    parser.add_argument("--data-dir", required=True, help="Benchmark data directory")
    parser.add_argument(
        "--output",
        default=None,
        help="Output exclude file path (default: <data-dir>/dapars2_single_group_exclude.txt)",
    )
    parser.add_argument("--group1", default="Sample_1", help="Expected group1 prefix")
    parser.add_argument("--group2", default="Sample_2", help="Expected group2 prefix")
    args = parser.parse_args()

    data_dir = os.path.abspath(args.data_dir)
    output_path = args.output or os.path.join(data_dir, "dapars2_single_group_exclude.txt")

    input_patterns = [
        # Current benchmark-final layout.
        os.path.join(
            data_dir,
            "result",
            "benchmark",
            "sim_bam",
            "dapars2",
            "*",
            "DaPars2_Result_all_chrs.txt",
        ),
        # Legacy layout kept for migration compatibility.
        os.path.join(data_dir, "sim_bam_result", "dapars2", "*", "DaPars2_Result_all_chrs.txt"),
    ]

    input_files = []
    for pattern in input_patterns:
        input_files.extend(glob.glob(pattern))
    input_files = sorted(set(input_files))

    excluded = []
    for path in input_files:
        sample_name = os.path.basename(os.path.dirname(path))
        header_columns = read_header_columns(path)
        if not has_required_groups(header_columns, args.group1, args.group2):
            excluded.append((f"dapars2/{sample_name}", "single-group"))

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as wf:
        wf.write("# sample_key\treason\n")
        for sample_key, reason in excluded:
            wf.write(f"{sample_key}\t{reason}\n")

    if not input_files:
        print("warning=no_dapars2_inputs_found")
        print("searched_patterns=" + ";".join(input_patterns))
    print(f"total_files={len(input_files)}")
    print(f"excluded_single_group={len(excluded)}")
    print(f"output={output_path}")


if __name__ == "__main__":
    main()
