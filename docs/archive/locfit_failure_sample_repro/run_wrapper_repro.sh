#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
test_root="${repo_root}/docs/archive/locfit_failure_sample_repro"
log_dir="${test_root}/logs"
mkdir -p "${log_dir}"

tool="sierra"
sample="Microwell_human_embryo_SRR9887755_mm10_pas1_gn5000_rep2"

gt_pas="${repo_root}/data/sim_data/sim_pas/mm10_sim_pas_gn5000_rep1.bed"
gt_mtx="${repo_root}/data/sim_data/bam/${sample}.bam.expr.tsv"
pd_pas="${repo_root}/data/result/benchmark/sim_bam/${tool}/${sample}/pas_fix.bed"
pd_mtx="${repo_root}/data/result/benchmark/sim_bam/${tool}/${sample}/pas_counts.tsv"
output_dir="${repo_root}/data/result/performance"
wrapper_script="${repo_root}/scripts/05_benchmark_analysis/scripts/calculate_benchmark_performance.R"

conda_bin="${CONDA_BIN:-${APABENCHMARK_CONDA_EXE:-conda}}"
log_file="${log_dir}/wrapper_repro.log"

echo "repo_root=${repo_root}"
echo "tool=${tool}"
echo "sample=${sample}"
echo "log_file=${log_file}"
echo "conda_bin=${conda_bin}"

for path in "${gt_pas}" "${gt_mtx}" "${pd_pas}" "${pd_mtx}" "${wrapper_script}"; do
  if [[ ! -f "${path}" ]]; then
    echo "Missing input file: ${path}" >&2
    exit 1
  fi
done

if ! command -v "${conda_bin}" >/dev/null 2>&1; then
  echo "Conda binary not found: ${conda_bin}" >&2
  exit 1
fi

set +e
"${conda_bin}" run -n dexseq Rscript "${wrapper_script}" \
  --pd_pas "${pd_pas}" \
  --pd_mtx "${pd_mtx}" \
  --gt_pas "${gt_pas}" \
  --gt_mtx "${gt_mtx}" \
  --tool "${tool}" \
  --sample "${sample}" \
  --output_dir "${output_dir}" \
  --dexseq_split_num 6 2>&1 | tee "${log_file}"
status=${PIPESTATUS[0]}
set -e

echo "exit_code=${status}"
exit "${status}"
