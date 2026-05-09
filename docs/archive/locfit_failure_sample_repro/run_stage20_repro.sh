#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
test_root="${repo_root}/docs/archive/locfit_failure_sample_repro"
log_dir="${test_root}/logs"
mkdir -p "${log_dir}"

tool="sierra"
sample="Microwell_human_embryo_SRR9887755_mm10_pas1_gn5000_rep2"
default_bundle="${repo_root}/data/result/performance/${tool}/${sample}_stage10_bundle.rds"
stage10_bundle="${STAGE10_BUNDLE:-${default_bundle}}"

stage20_script="${repo_root}/scripts/05_benchmark_analysis/scripts/20_stage_apa_and_te.R"
log_file="${log_dir}/stage20_repro.log"
conda_bin="${CONDA_BIN:-${APABENCHMARK_CONDA_EXE:-conda}}"

echo "repo_root=${repo_root}"
echo "stage10_bundle=${stage10_bundle}"
echo "log_file=${log_file}"
echo "conda_bin=${conda_bin}"

if [[ ! -f "${stage10_bundle}" ]]; then
  echo "Missing stage10 bundle: ${stage10_bundle}" >&2
  exit 1
fi

if ! command -v "${conda_bin}" >/dev/null 2>&1; then
  echo "Conda binary not found: ${conda_bin}" >&2
  exit 1
fi

set +e
"${conda_bin}" run -n dexseq Rscript "${stage20_script}" --stage10_bundle "${stage10_bundle}" 2>&1 | tee "${log_file}"
status=${PIPESTATUS[0]}
set -e

echo "exit_code=${status}"
exit "${status}"
