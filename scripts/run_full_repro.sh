#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./scripts/run_full_repro.sh --data-root /path/to/apabenchmark_final [options]

Runs the connected full APA benchmark flow:
  01_generate_annotation -> 02_peak_analysis -> 03_data_simulation
  -> 04_benchmark -> 05_benchmark_analysis -> connected 06_plotting topics.

Options:
  --data-root PATH       Project root containing data/.
  --cores N             Cores for Snakemake stages. Default: 8.
  --conda-exe PATH      Conda executable. Default: APABENCHMARK_CONDA_EXE or conda on PATH.
  --driver-env NAME     Conda env used to launch Snakemake and stage 01. Default: apasim.
  --plot-env NAME       Conda env used to run plotting scripts. Default: benchmark_plotting.
  --mm10-fasta PATH     Mouse genome FASTA for simulation.
  --pipeline-mode MODE  sim, raw, both, or raw_parallel for stage 05. Default: both.
  --raw-run-id ID       Raw run id for stage 05/plots. Default: default.
  --skip-01             Skip annotation generation.
  --skip-02             Skip peak analysis.
  --skip-03             Skip data simulation.
  --skip-04             Skip benchmark tool runs.
  --skip-05             Skip benchmark analysis.
  --skip-plots          Skip plotting.
  --with-raw-plots      Run raw_data_performance plotting workflow.
  --dry-run             Run Snakemake dry-runs where possible.
  -h, --help            Show this help.
EOF
}

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
data_root="${APABENCHMARK_FINAL_ROOT:-}"
cores="8"
conda_exe="${APABENCHMARK_CONDA_EXE:-}"
driver_env="${APABENCHMARK_DRIVER_ENV:-apasim}"
plot_env="${APABENCHMARK_PLOT_ENV:-benchmark_plotting}"
pipeline_mode="${APABENCHMARK_PIPELINE_MODE:-both}"
raw_run_id="${APABENCHMARK_RAW_RUN_ID:-default}"
mm10_fasta="${APABENCHMARK_MM10_GENOME_FASTA:-}"
run_01="1"
run_02="1"
run_03="1"
run_04="1"
run_05="1"
run_plots="1"
with_raw_plots="0"
dry_run="0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-root) data_root="$2"; shift 2 ;;
    --cores) cores="$2"; shift 2 ;;
    --conda-exe) conda_exe="$2"; shift 2 ;;
    --driver-env) driver_env="$2"; shift 2 ;;
    --plot-env) plot_env="$2"; shift 2 ;;
    --mm10-fasta) mm10_fasta="$2"; shift 2 ;;
    --pipeline-mode) pipeline_mode="$2"; shift 2 ;;
    --raw-run-id) raw_run_id="$2"; shift 2 ;;
    --skip-01) run_01="0"; shift ;;
    --skip-02) run_02="0"; shift ;;
    --skip-03) run_03="0"; shift ;;
    --skip-04) run_04="0"; shift ;;
    --skip-05) run_05="0"; shift ;;
    --skip-plots) run_plots="0"; shift ;;
    --with-raw-plots) with_raw_plots="1"; shift ;;
    --dry-run) dry_run="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

case "${pipeline_mode}" in
  sim|raw|both|raw_parallel) ;;
  *) echo "Invalid --pipeline-mode: ${pipeline_mode}" >&2; exit 2 ;;
esac

if [[ -z "${data_root}" ]]; then
  echo "Missing --data-root. Provide a project root containing data/." >&2
  exit 2
fi
if [[ ! -d "${data_root}" ]]; then
  echo "Data root not found: ${data_root}" >&2
  exit 2
fi
project_root="$(cd "${data_root}" && pwd -P)"
if [[ "$(basename "${project_root}")" == "data" ]]; then
  project_root="$(cd "${project_root}/.." && pwd -P)"
fi
if [[ ! -d "${project_root}/data" ]]; then
  echo "Project data directory not found: ${project_root}/data" >&2
  exit 2
fi

if [[ "${project_root}" != "${repo_root}" && ( "${run_01}" == "1" || "${run_02}" == "1" || "${run_03}" == "1" || "${run_04}" == "1" ) ]]; then
  cat >&2 <<EOF
Stages 01/02/03/04 use ../../data paths relative to this repository checkout.
For a full upstream run, --data-root must be the repository root:
  ${repo_root}
Current resolved project root:
  ${project_root}

Use a checkout with data/ at its root, or skip upstream stages and run only
05/06 with --skip-01 --skip-02 --skip-03 --skip-04.
EOF
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

require_file() {
  local path="$1"
  local label="$2"
  if [[ ! -e "${path}" ]]; then
    echo "Missing ${label}: ${path}" >&2
    exit 2
  fi
}

require_dir() {
  local path="$1"
  local label="$2"
  if [[ ! -d "${path}" ]]; then
    echo "Missing ${label}: ${path}" >&2
    exit 2
  fi
}

if [[ "${run_01}" == "1" ]]; then
  require_file "${project_root}/data/raw_data/annotations/gencode.vM25.polyAs.gtf" "stage 01 mouse Gencode polyA input"
  require_file "${project_root}/data/raw_data/annotations/atlas.clusters.3.0.GRCm38.GENCODE_M25.bed.gz" "stage 01 mouse PolyASite input"
  require_file "${project_root}/data/raw_data/annotations/MousePas/mm10.PAS.main.tsv" "stage 01 mouse PolyA_DB input"
  require_file "${project_root}/data/raw_data/annotations/gencode.vM25.annotation.bed" "stage 01 mouse Gencode BED input"
  require_file "${project_root}/data/raw_data/annotations/gencode.v40.polyAs.gtf" "stage 01 human Gencode polyA input"
  require_file "${project_root}/data/raw_data/annotations/atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz" "stage 01 human PolyASite input"
  require_file "${project_root}/data/raw_data/annotations/HumanPas/hg38.PAS.main.tsv" "stage 01 human PolyA_DB input"
  require_file "${project_root}/data/raw_data/annotations/gencode.v40.annotation.bed" "stage 01 human Gencode BED input"
else
  if [[ "${run_03}" == "1" ]]; then
    require_file "${project_root}/data/int_data/annotations/mouse_integrated_pas.bed" "stage 01 mouse annotation output"
  fi
  if [[ "${run_02}" == "1" ]]; then
    require_file "${project_root}/data/int_data/annotations/mouse_pas_500.bed" "stage 01 filtered mouse PAS output"
  fi
fi
if [[ "${run_02}" == "1" ]]; then
  require_dir "${project_root}/data/raw_data/bam_for_peak_extraction" "stage 02 raw BAM input directory"
fi

if [[ "${run_03}" == "1" ]]; then
  if [[ -z "${mm10_fasta}" ]]; then
    echo "Missing --mm10-fasta or APABENCHMARK_MM10_GENOME_FASTA for stage 03." >&2
    exit 2
  fi
  require_file "${mm10_fasta}" "stage 03 mouse genome FASTA"
  export APABENCHMARK_MM10_GENOME_FASTA="$(cd "$(dirname "${mm10_fasta}")" && pwd -P)/$(basename "${mm10_fasta}")"
fi

export APABENCHMARK_FINAL_ROOT="${project_root}"
export APABENCHMARK_BENCHMARK5_DIR="${repo_root}/scripts/05_benchmark_analysis"
export APABENCHMARK_CONDA_EXE="${conda_exe}"
export APABENCHMARK_PIPELINE_MODE="${pipeline_mode}"
export APABENCHMARK_RAW_RUN_ID="${raw_run_id}"

echo "[full-repro] repo root: ${repo_root}"
echo "[full-repro] project root: ${project_root}"
echo "[full-repro] pipeline mode: ${pipeline_mode}"

snakemake_extra=()
if [[ "${dry_run}" == "1" ]]; then
  snakemake_extra+=("-n")
fi

if [[ "${run_01}" == "1" ]]; then
  if [[ "${dry_run}" == "1" ]]; then
    echo "[full-repro] stage 01 annotation generation preflight complete"
    if [[ ( "${run_02}" == "1" && ! -e "${project_root}/data/int_data/annotations/mouse_pas_500.bed" ) || \
          ( "${run_03}" == "1" && ! -e "${project_root}/data/int_data/annotations/mouse_integrated_pas.bed" ) ]]; then
      cat <<EOF
[full-repro] downstream dry-runs skipped because stage 01 outputs are absent.
Run without --dry-run to generate annotation outputs, or rerun dry-run after
prepared annotation files exist under data/int_data/annotations/.
EOF
      exit 0
    fi
  else
    echo "[full-repro] stage 01 annotation generation"
    "${conda_exe}" run -n "${driver_env}" bash \
      "${repo_root}/scripts/01_generate_annotation/run_annotation.sh"
  fi
fi

if [[ "${run_03}" == "1" ]]; then
  require_file "${project_root}/data/int_data/annotations/mouse_integrated_pas.bed" "stage 01 mouse annotation output"
fi
if [[ "${run_02}" == "1" ]]; then
  require_file "${project_root}/data/int_data/annotations/mouse_pas_500.bed" "stage 01 filtered mouse PAS output"
fi

if [[ "${run_02}" == "1" ]]; then
  echo "[full-repro] stage 02 peak analysis"
  (cd "${repo_root}/scripts/02_peak_analysis" && "${conda_exe}" run -n "${driver_env}" \
    snakemake --cores "${cores}" --use-conda "${snakemake_extra[@]}")
fi

if [[ "${run_03}" == "1" ]]; then
  echo "[full-repro] stage 03 data simulation"
  (cd "${repo_root}/scripts/03_data_simulation" && "${conda_exe}" run -n "${driver_env}" \
    snakemake --cores "${cores}" --use-conda "${snakemake_extra[@]}")
fi

if [[ "${run_04}" == "1" ]]; then
  echo "[full-repro] stage 04 benchmark tool runs"
  (cd "${repo_root}/scripts/04_benchmark" && bash generate_sample_config.sh "${project_root}" ./snakemake_profile/sample.yaml)
  (cd "${repo_root}/scripts/04_benchmark" && "${conda_exe}" run -n "${driver_env}" \
    snakemake --cores "${cores}" --profile snakemake_profile \
    --singularity-args "--bind ${project_root}:${project_root} --cpus ${cores}" \
    "${snakemake_extra[@]}")
fi

if [[ "${run_05}" == "1" ]]; then
  echo "[full-repro] stage 05 benchmark analysis"
  stage05_args=(
    run -n "${driver_env}"
    snakemake
    -s "${repo_root}/scripts/05_benchmark_analysis/Snakefile"
    --cores "${cores}"
  )
  if [[ "${dry_run}" == "1" ]]; then
    stage05_args+=("-n")
  fi
  "${conda_exe}" "${stage05_args[@]}"
fi

if [[ "${dry_run}" == "1" ]]; then
  echo "[full-repro] dry run complete; plotting skipped."
  exit 0
fi

if [[ "${run_plots}" == "1" ]]; then
  echo "[full-repro] stage 06 sim plots"
  "${conda_exe}" run -n "${plot_env}" python \
    "${repo_root}/scripts/06_plotting/sim_data_performance/scripts/run_all.py" \
    --data-root "${project_root}" \
    --force-prepare

  echo "[full-repro] stage 06 peak model plots"
  "${conda_exe}" run -n "${plot_env}" python \
    "${repo_root}/scripts/06_plotting/peak_model_params/scripts/run_all.py" \
    --data-root "${project_root}"

  echo "[full-repro] stage 06 computation resource plots"
  APABENCHMARK_CP_RESOURCE_ROOT="${project_root}/data/result/cp_resource/pf_bam" \
    "${conda_exe}" run -n "${plot_env}" python \
      "${repo_root}/scripts/06_plotting/computation_resource_consumption/scripts/run_all.py"

  if [[ "${with_raw_plots}" == "1" ]]; then
    echo "[full-repro] stage 06 raw plots"
    "${conda_exe}" run -n "${plot_env}" python \
      "${repo_root}/scripts/06_plotting/raw_data_performance/scripts/run_all.py" \
      --data-root "${project_root}" \
      --raw-run-id "${raw_run_id}"
  fi
fi

echo "[full-repro] done"
