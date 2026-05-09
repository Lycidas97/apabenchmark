#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Usage:
#   bash generate_sample_config.sh [PROJECT_ROOT] [OUTPUT_CONFIG]
# Example:
#   bash generate_sample_config.sh /path/to/apabenchmark_final \
#     /path/to/apabenchmark_final/scripts/04_benchmark/snakemake_profile/sample.yaml
PROJECT_ROOT="${1:-}"
OUTPUT_CONFIG="${2:-${SCRIPT_DIR}/snakemake_profile/sample.yaml}"
READLEN_SAMPLE_SIZE=2000

has_required_benchmark_inputs() {
    local root="$1"
    [ -d "${root}/data/sim_data/bam" ] || return 1
    [ -d "${root}/data/sim_data/bam_for_cprsb" ] || return 1
    if [ -d "${root}/data/raw_data/bam/raw_bam" ] || [ -d "${root}/data/int_data/bam_to_detect_pas" ]; then
        return 0
    fi
    return 1
}

if [ -z "${PROJECT_ROOT}" ]; then
    DEFAULT_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
    if has_required_benchmark_inputs "${DEFAULT_ROOT}"; then
        PROJECT_ROOT="${DEFAULT_ROOT}"
    else
        PROJECT_ROOT="${DEFAULT_ROOT}"
    fi
fi

RAW_BAM_DIR="${PROJECT_ROOT}/data/raw_data/bam/raw_bam"
if [ ! -d "${RAW_BAM_DIR}" ] && [ -d "${PROJECT_ROOT}/data/int_data/bam_to_detect_pas" ]; then
    RAW_BAM_DIR="${PROJECT_ROOT}/data/int_data/bam_to_detect_pas"
fi
SIM_BAM_DIR="${PROJECT_ROOT}/data/sim_data/bam"
PF_BAM_DIR="${PROJECT_ROOT}/data/sim_data/bam_for_cprsb"

SAMTOOLS_BIN="$(command -v samtools || true)"
if [ -z "${SAMTOOLS_BIN}" ]; then
    echo "samtools not found on PATH" >&2
    exit 1
fi

declare -A READLEN_CACHE
load_cache_file() {
    local cache_file="$1"
    [ -f "${cache_file}" ] || return 0
    while IFS=$'\t' read -r cached_sample cached_len; do
        if [ -n "${cached_sample}" ] && [ -n "${cached_len}" ]; then
            READLEN_CACHE["${cached_sample}"]="${cached_len}"
        fi
    done < <(
        awk '
            /^  [^ ].*:$/ {sample=$1; sub(/:$/, "", sample)}
            /^    read_length:/ && sample != "" {print sample "\t" $2}
        ' "${cache_file}"
    )
}

load_cache_file "${OUTPUT_CONFIG}"
ROOT_DEFAULT_CACHE="${PROJECT_ROOT}/scripts/04_benchmark/snakemake_profile/sample.yaml"
if [ "${ROOT_DEFAULT_CACHE}" != "${OUTPUT_CONFIG}" ]; then
    load_cache_file "${ROOT_DEFAULT_CACHE}"
fi

infer_read_length() {
    local bam_file="$1"
    local sample_name="$2"
    local rl=""

    if [ -n "${READLEN_CACHE[${sample_name}]:-}" ]; then
        echo "${READLEN_CACHE[${sample_name}]}"
        return
    fi

    # Fast path for names carrying rl metadata, e.g. *_rl100_*
    if [[ "${sample_name}" =~ _rl([0-9]+)_ ]]; then
        rl="${BASH_REMATCH[1]}"
    fi

    if [ -n "${rl}" ]; then
        echo "${rl}"
        return
    fi

    local inferred=""
    set +o pipefail
    inferred="$("${SAMTOOLS_BIN}" view "${bam_file}" 2>/dev/null | awk -v sample_size="${READLEN_SAMPLE_SIZE}" '
        $10 != "*" && length($10) > 0 {
            len[length($10)]++
            total++
            if (total >= sample_size) {
                exit
            }
        }
        END {
            max_count = 0
            mode_len = ""
            for (l in len) {
                if (len[l] > max_count || (len[l] == max_count && (mode_len == "" || l + 0 < mode_len + 0))) {
                    max_count = len[l]
                    mode_len = l
                }
            }
            if (mode_len != "") {
                print mode_len
            }
        }
    ')"
    set -o pipefail
    echo "${inferred}"
}

has_sim_genome_token() {
    local sample_name="$1"
    [[ "${sample_name}" =~ (^|_)(mm10|hg38)(_|$) ]]
}

append_samples() {
    local section_name="$1"
    local bam_dir="$2"
    local sample_count=0
    local section_entries=""
    local -a invalid_sim_samples=()

    if [ ! -d "$bam_dir" ]; then
        echo "Warning: BAM directory not found: $bam_dir" >&2
        echo "${section_name}: {}" >> "${OUTPUT_CONFIG}"
        echo "[${section_name}] ${sample_count} samples from ${bam_dir}"
        return
    fi

    shopt -s nullglob
    local bai_files=("${bam_dir}"/*.bam.bai)
    IFS=$'\n' bai_files=($(printf '%s\n' "${bai_files[@]}" | sort))
    local bai_file
    for bai_file in "${bai_files[@]}"; do
        local bam_file sample read_length
        bam_file="${bai_file%.bai}"
        sample="$(basename "$bai_file" .bam.bai)"

        if [ "${section_name}" = "sim_sample" ] && ! has_sim_genome_token "${sample}"; then
            invalid_sim_samples+=("${sample}")
            continue
        fi

        read_length="$(infer_read_length "$bam_file" "$sample")"

        if [ -z "$read_length" ]; then
            echo "Warning: cannot detect read length for sample $sample" >&2
            continue
        fi

        section_entries+="  ${sample}:"$'\n'
        section_entries+="    read_length: ${read_length}"$'\n'
        sample_count=$((sample_count + 1))
    done
    shopt -u nullglob

    if [ "${section_name}" = "sim_sample" ] && [ "${#invalid_sim_samples[@]}" -gt 0 ]; then
        echo "Error: sim sample names must include genome token mm10/hg38 (e.g. *_mm10_pas1_gn5000_rep1)." >&2
        echo "Invalid sim samples:" >&2
        printf '  - %s\n' "${invalid_sim_samples[@]}" >&2
        return 1
    fi

    if [ "${sample_count}" -eq 0 ]; then
        echo "${section_name}: {}" >> "${OUTPUT_CONFIG}"
    else
        echo "${section_name}:" >> "${OUTPUT_CONFIG}"
        printf "%s" "${section_entries}" >> "${OUTPUT_CONFIG}"
    fi
    echo "[${section_name}] ${sample_count} samples from ${bam_dir}"
}

mkdir -p "$(dirname "${OUTPUT_CONFIG}")"
: > "${OUTPUT_CONFIG}"
if ! append_samples "raw_sample" "${RAW_BAM_DIR}" \
    || ! append_samples "sim_sample" "${SIM_BAM_DIR}" \
    || ! append_samples "pf_sample" "${PF_BAM_DIR}"; then
    rm -f "${OUTPUT_CONFIG}"
    exit 1
fi

echo "Config generated: ${OUTPUT_CONFIG}"
echo "Project root used: ${PROJECT_ROOT}"
