#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

SOURCE_DIR="${1:-}"
MAP_FILE="${2:-${SCRIPT_DIR}/snakemake_profile/raw_bam_name_map.tsv}"
TARGET_DIR="${3:-${PROJECT_ROOT}/data/raw_data/bam/raw_bam}"
GT_DIR="${PROJECT_ROOT}/data/raw_data/raw_bam"
LEGACY_DIR="${PROJECT_ROOT}/data/int_data/bam_to_detect_pas"

if [ -z "${SOURCE_DIR}" ]; then
    echo "Usage: bash link_raw_bam_inputs.sh SOURCE_DIR [MAP_FILE] [TARGET_DIR]" >&2
    exit 2
fi

if [ ! -d "${SOURCE_DIR}" ]; then
    echo "Source raw_bam directory does not exist: ${SOURCE_DIR}" >&2
    exit 1
fi

if [ ! -f "${MAP_FILE}" ]; then
    echo "Mapping file does not exist: ${MAP_FILE}" >&2
    exit 1
fi

mkdir -p "${TARGET_DIR}"

safe_link() {
    local src="$1"
    local dst="$2"
    if [ -e "${dst}" ] && [ ! -L "${dst}" ]; then
        echo "Refusing to overwrite non-symlink file: ${dst}" >&2
        return 1
    fi
    ln -sfn "${src}" "${dst}"
}

mapped_count=0
while IFS=$'\t' read -r canonical source || [ -n "${canonical:-}" ]; do
    if [ -z "${canonical:-}" ] || [[ "${canonical}" =~ ^# ]]; then
        continue
    fi
    if [ -z "${source:-}" ]; then
        source="${canonical}"
    fi

    for suffix in "bam" "bam.bai" "barcode_list.txt"; do
        src="${SOURCE_DIR}/${source}.${suffix}"
        dst="${TARGET_DIR}/${canonical}.${suffix}"
        if [ -f "${src}" ]; then
            safe_link "${src}" "${dst}"
        else
            echo "Warning: source file not found, skip link: ${src}" >&2
        fi
    done

    if [ ! -f "${GT_DIR}/${canonical}/bam.expr.tsv" ] || [ ! -f "${GT_DIR}/${canonical}/pas.bed" ]; then
        echo "Warning: missing raw ground truth for ${canonical} under ${GT_DIR}" >&2
    fi
    mapped_count=$((mapped_count + 1))
done < "${MAP_FILE}"

mkdir -p "$(dirname "${LEGACY_DIR}")"
if [ -e "${LEGACY_DIR}" ] && [ ! -L "${LEGACY_DIR}" ]; then
    echo "Warning: legacy raw BAM path exists and is not symlink, skipped: ${LEGACY_DIR}" >&2
else
    ln -sfn "${TARGET_DIR}" "${LEGACY_DIR}"
fi

echo "Mapped samples: ${mapped_count}"
echo "Mapped raw BAM dir: ${TARGET_DIR}"
echo "Legacy compatible path: ${LEGACY_DIR}"
