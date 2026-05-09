#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

python mouse_pas_integration.py
python human_pas_integration.py
python filter_pas.py
python summarize_annotation_stage_counts.py
