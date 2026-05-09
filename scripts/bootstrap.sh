#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./scripts/bootstrap.sh [options]

Prepare the APABenchmark checkout for a first run. By default this creates
conda environments, downloads public annotation inputs, downloads benchmark SIF
images, and builds public tool-reference files.

Options:
  --data-root PATH             Project root containing data/. Default: repository root.
  --conda-exe PATH             Conda executable passed to setup_envs.sh.
  --skip-envs                  Do not create conda environments.
  --skip-public-inputs         Do not download public annotation inputs.
  --skip-sif                   Do not download benchmark SIF images.
  --skip-tool-references       Do not build tool-reference files.
  --skip-scape                 Build tool references without SCAPE UTR/intron references.
  --force                      Recreate downloaded/generated public inputs where supported.
  -h, --help                   Show this help.
EOF
}

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
data_root="${APABENCHMARK_FINAL_ROOT:-${repo_root}}"
conda_exe="${APABENCHMARK_CONDA_EXE:-}"
run_envs="1"
run_public_inputs="1"
run_sif="1"
run_tool_references="1"
skip_scape="0"
force="0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-root) data_root="$2"; shift 2 ;;
    --conda-exe) conda_exe="$2"; shift 2 ;;
    --skip-envs) run_envs="0"; shift ;;
    --skip-public-inputs) run_public_inputs="0"; shift ;;
    --skip-sif) run_sif="0"; shift ;;
    --skip-tool-references) run_tool_references="0"; shift ;;
    --skip-scape) skip_scape="1"; shift ;;
    --force) force="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

project_root="$(cd "${data_root}" && pwd -P)"
mkdir -p \
  "${project_root}/data/raw_data/annotations" \
  "${project_root}/data/raw_data/bam_for_peak_extraction" \
  "${project_root}/data/raw_data/tool_reference" \
  "${project_root}/data/int_data" \
  "${project_root}/data/sim_data" \
  "${project_root}/data/result" \
  "${repo_root}/scripts/04_benchmark/sif"

echo "[bootstrap] repository root: ${repo_root}"
echo "[bootstrap] project root: ${project_root}"

if [[ "${run_envs}" == "1" ]]; then
  env_args=()
  if [[ -n "${conda_exe}" ]]; then
    env_args+=(--conda-exe "${conda_exe}")
  fi
  "${repo_root}/scripts/setup_envs.sh" "${env_args[@]}"
fi

if [[ "${run_public_inputs}" == "1" ]]; then
  public_args=(--data-root "${project_root}")
  if [[ "${force}" == "1" ]]; then
    public_args+=(--force)
  fi
  python3 "${repo_root}/scripts/prepare_public_inputs.py" "${public_args[@]}"
fi

if [[ "${run_sif}" == "1" ]]; then
  sif_args=()
  if [[ "${force}" == "1" ]]; then
    sif_args+=(--force)
  fi
  bash "${repo_root}/scripts/04_benchmark/download_sif.sh" "${sif_args[@]}"
fi

if [[ "${run_tool_references}" == "1" ]]; then
  ref_args=(--output-dir "${project_root}/data/raw_data/tool_reference")
  if [[ "${force}" == "1" ]]; then
    ref_args+=(--force)
  fi
  if [[ "${skip_scape}" == "1" ]]; then
    ref_args+=(--skip-scape)
  else
    ref_args+=(--scape-image "${repo_root}/scripts/04_benchmark/sif/scape.sif")
  fi
  python3 "${repo_root}/scripts/04_benchmark/prepare_tool_references.py" "${ref_args[@]}"
fi

cat <<EOF
[bootstrap] done

Next steps:
  1. Place prepared BAM/BAI files under:
       ${project_root}/data/raw_data/bam_for_peak_extraction/
  2. Run:
       ./scripts/run_full_repro.sh --data-root "${project_root}" --mm10-fasta /path/to/GRCm38.p6.genome.fa --cores 8 --dry-run
EOF
