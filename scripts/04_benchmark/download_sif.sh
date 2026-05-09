#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/04_benchmark/download_sif.sh [options]

Download benchmark Singularity/Apptainer images into scripts/04_benchmark/sif.

Options:
  --force      Re-download images that already exist.
  -h, --help   Show this help.
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
sif_dir="${script_dir}/sif"
force="0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --force) force="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

mkdir -p "${sif_dir}"

download() {
  local name="$1"
  local url="$2"
  local target="${sif_dir}/${name}"
  if [[ -s "${target}" && "${force}" != "1" ]]; then
    echo "[download-sif] exists: ${target}"
    return
  fi
  echo "[download-sif] ${name}"
  if command -v curl >/dev/null 2>&1; then
    curl -L --fail --retry 3 -o "${target}.tmp" "${url}"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "${target}.tmp" "${url}"
  else
    echo "Cannot find curl or wget." >&2
    exit 2
  fi
  mv "${target}.tmp" "${target}"
}

base_url="https://bis.zju.edu.cn/nextcloud/s/apabenchmark_sif/download?path=%2F&files="

download "infernape.sif" "${base_url}infernape.sif"
download "scape.sif" "${base_url}scape.sif"
download "scapa.sif" "${base_url}scapa.sif"
download "scapatrap.sif" "${base_url}scapatrap.sif"
download "sierra.sif" "${base_url}sierra.sif"
download "scapture.sif" "${base_url}scapture.sif"
