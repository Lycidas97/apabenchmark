#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./scripts/setup_envs.sh [options]

Create or update the conda environments used by the APABenchmark workflow.

Options:
  --conda-exe PATH                  Conda executable. Default: APABENCHMARK_CONDA_EXE or conda on PATH.
  --include-benchmark-performance   Also create/update scripts/envs/benchmark_performance.yml.
  --update-existing                 Update existing environments instead of skipping them.
  -h, --help                        Show this help.
EOF
}

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
conda_exe="${APABENCHMARK_CONDA_EXE:-}"
include_benchmark_performance="0"
update_existing="0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --conda-exe) conda_exe="$2"; shift 2 ;;
    --include-benchmark-performance) include_benchmark_performance="1"; shift ;;
    --update-existing) update_existing="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "${conda_exe}" ]]; then
  if command -v conda >/dev/null 2>&1; then
    conda_exe="$(command -v conda)"
  else
    echo "Cannot find conda. Pass --conda-exe or set APABENCHMARK_CONDA_EXE." >&2
    exit 2
  fi
fi

env_specs=(
  "apasim|scripts/envs/apasim.yml"
  "dexseq|scripts/envs/dexseq.yml"
  "scmapa_r|scripts/envs/scmapa_r.yml"
  "benchmark_plotting|scripts/envs/benchmark_plotting.yml"
)

if [[ "${include_benchmark_performance}" == "1" ]]; then
  env_specs+=("benchmark_performance|scripts/envs/benchmark_performance.yml")
fi

env_exists() {
  local env_name="$1"
  "${conda_exe}" env list | awk '{print $1}' | grep -Fxq "${env_name}"
}

echo "[setup-envs] repository root: ${repo_root}"
echo "[setup-envs] conda executable: ${conda_exe}"

for env_spec in "${env_specs[@]}"; do
  env_name="${env_spec%%|*}"
  rel_env_file="${env_spec#*|}"
  env_file="${repo_root}/${rel_env_file}"
  if [[ ! -f "${env_file}" ]]; then
    echo "Missing environment file: ${env_file}" >&2
    exit 2
  fi

  if env_exists "${env_name}"; then
    if [[ "${update_existing}" == "1" ]]; then
      echo "[setup-envs] updating ${env_name} from ${rel_env_file}"
      "${conda_exe}" env update -n "${env_name}" -f "${env_file}" --prune
    else
      echo "[setup-envs] ${env_name} already exists; skipping"
    fi
  else
    echo "[setup-envs] creating ${env_name} from ${rel_env_file}"
    "${conda_exe}" env create -n "${env_name}" -f "${env_file}"
  fi
done

echo "[setup-envs] done"
