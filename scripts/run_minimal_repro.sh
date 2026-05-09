#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./scripts/run_minimal_repro.sh --data-root /path/to/apabenchmark_final [options]

Default action:
  1. Run scripts/05_benchmark_analysis/Snakefile in sim mode.
  2. Run scripts/06_plotting/sim_data_performance/scripts/run_all.py.

Options:
  --data-root PATH       Project/data root containing data/ or the data dir itself.
  --cores N             Snakemake cores. Default: 8.
  --conda-exe PATH      Conda executable. Default: APABENCHMARK_CONDA_EXE or conda on PATH.
  --driver-env NAME     Conda env used to launch Snakemake. Default: apasim.
  --plot-env NAME       Conda env used to run plotting scripts. Default: benchmark_plotting.
  --pipeline-mode MODE  sim, raw, both, or raw_parallel. Default: sim.
  --raw-run-id ID       Raw run id for raw plot workflow. Default: default.
  --skip-analysis       Do not run 05_benchmark_analysis.
  --plots-only          Skip analysis and require existing plot intermediates.
  --force-prepare       Force regeneration of sim plot intermediate tables.
  --with-raw-plots      Also run raw_data_performance plotting workflow.
  --dry-run             Run only the Snakemake DAG dry run for 05_benchmark_analysis.
  -h, --help            Show this help.
EOF
}

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
data_root="${APABENCHMARK_FINAL_ROOT:-}"
cores="8"
conda_exe="${APABENCHMARK_CONDA_EXE:-}"
driver_env="${APABENCHMARK_DRIVER_ENV:-apasim}"
plot_env="${APABENCHMARK_PLOT_ENV:-benchmark_plotting}"
pipeline_mode="${APABENCHMARK_PIPELINE_MODE:-sim}"
raw_run_id="${APABENCHMARK_RAW_RUN_ID:-default}"
run_analysis="1"
plots_only="0"
force_prepare="0"
with_raw_plots="0"
dry_run="0"

if [[ -z "${data_root}" && -d "${repo_root}/data" ]]; then
  data_root="${repo_root}"
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-root)
      data_root="$2"
      shift 2
      ;;
    --cores)
      cores="$2"
      shift 2
      ;;
    --conda-exe)
      conda_exe="$2"
      shift 2
      ;;
    --driver-env)
      driver_env="$2"
      shift 2
      ;;
    --plot-env)
      plot_env="$2"
      shift 2
      ;;
    --pipeline-mode)
      pipeline_mode="$2"
      shift 2
      ;;
    --raw-run-id)
      raw_run_id="$2"
      shift 2
      ;;
    --skip-analysis)
      run_analysis="0"
      shift
      ;;
    --plots-only)
      run_analysis="0"
      plots_only="1"
      shift
      ;;
    --force-prepare)
      force_prepare="1"
      shift
      ;;
    --with-raw-plots)
      with_raw_plots="1"
      shift
      ;;
    --dry-run)
      dry_run="1"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

case "${pipeline_mode}" in
  sim|raw|both|raw_parallel) ;;
  *)
    echo "Invalid --pipeline-mode: ${pipeline_mode}" >&2
    exit 2
    ;;
esac

if [[ -z "${data_root}" ]]; then
  echo "Missing --data-root. Provide a project root containing data/ or the data directory itself." >&2
  exit 2
fi

if [[ ! -d "${data_root}" ]]; then
  echo "Data root not found: ${data_root}" >&2
  exit 2
fi

data_root="$(cd "${data_root}" && pwd -P)"
if [[ "$(basename "${data_root}")" == "data" ]]; then
  project_root="$(cd "${data_root}/.." && pwd -P)"
  data_dir="${data_root}"
else
  project_root="${data_root}"
  data_dir="${data_root}/data"
fi

if [[ ! -d "${data_dir}" ]]; then
  echo "Data directory not found: ${data_dir}" >&2
  exit 2
fi

if [[ -z "${conda_exe}" ]]; then
  if command -v conda >/dev/null 2>&1; then
    conda_exe="$(command -v conda)"
  else
    echo "Cannot find conda. Pass --conda-exe or set APABENCHMARK_CONDA_EXE." >&2
    exit 2
  fi
fi

export APABENCHMARK_FINAL_ROOT="${project_root}"
export APABENCHMARK_BENCHMARK5_DIR="${repo_root}/scripts/05_benchmark_analysis"
export APABENCHMARK_CONDA_EXE="${conda_exe}"
export APABENCHMARK_PIPELINE_MODE="${pipeline_mode}"
export APABENCHMARK_RAW_RUN_ID="${raw_run_id}"

echo "[minimal-repro] repo root: ${repo_root}"
echo "[minimal-repro] data root: ${project_root}"
echo "[minimal-repro] pipeline mode: ${pipeline_mode}"
echo "[minimal-repro] plot env: ${plot_env}"

if [[ "${run_analysis}" == "1" ]]; then
  snakemake_args=(
    run -n "${driver_env}"
    snakemake
    -s "${repo_root}/scripts/05_benchmark_analysis/Snakefile"
    --cores "${cores}"
  )
  if [[ "${dry_run}" == "1" ]]; then
    snakemake_args+=("-n")
  fi
  echo "[minimal-repro] running 05_benchmark_analysis"
  "${conda_exe}" "${snakemake_args[@]}"
else
  echo "[minimal-repro] skipping 05_benchmark_analysis"
fi

if [[ "${dry_run}" == "1" ]]; then
  echo "[minimal-repro] dry run complete; plotting skipped."
  exit 0
fi

sim_plot_args=(
  "${repo_root}/scripts/06_plotting/sim_data_performance/scripts/run_all.py"
  --data-root "${project_root}"
)
if [[ "${plots_only}" == "1" ]]; then
  sim_plot_args+=(--plots-only)
fi
if [[ "${force_prepare}" == "1" ]]; then
  sim_plot_args+=(--force-prepare)
fi

echo "[minimal-repro] running sim_data_performance plots"
"${conda_exe}" run -n "${plot_env}" python "${sim_plot_args[@]}"

if [[ "${with_raw_plots}" == "1" ]]; then
  echo "[minimal-repro] running raw_data_performance plots"
  "${conda_exe}" run -n "${plot_env}" python \
    "${repo_root}/scripts/06_plotting/raw_data_performance/scripts/run_all.py" \
    --data-root "${project_root}" \
    --raw-run-id "${raw_run_id}"
fi

echo "[minimal-repro] done"
