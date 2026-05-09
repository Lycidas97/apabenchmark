#!/usr/bin/env bash
set -euo pipefail

PID_FILE=""
PID=""
LOG_FILE=""
INTERVAL=20
RUN_ROOT="${APABENCHMARK_RAW_RUN_ROOT:-${APABENCHMARK_FINAL_ROOT:-$(pwd)}/data/result/raw/default}"
SAMPLE="Chromium_mouse_hipp_rep2"
PAIRS=(
  "30_Astro_Epen_vs_04_DG_IMN_Glut"
  "04_DG_IMN_Glut_vs_31_OPC_Oligo"
  "30_Astro_Epen_vs_31_OPC_Oligo"
)

usage() {
  cat <<'EOF'
Usage:
  monitor_retry_process.sh --log-file <path> [--pid <n> | --pid-file <path>] [--interval <sec>] [--run-root <path>]
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pid-file)
      PID_FILE="$2"
      shift 2
      ;;
    --pid)
      PID="$2"
      shift 2
      ;;
    --log-file)
      LOG_FILE="$2"
      shift 2
      ;;
    --interval)
      INTERVAL="$2"
      shift 2
      ;;
    --run-root)
      RUN_ROOT="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown arg: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ -z "$LOG_FILE" ]]; then
  echo "--log-file is required" >&2
  exit 2
fi

if [[ -z "$PID_FILE" && -z "$PID" ]]; then
  echo "one of --pid-file or --pid is required" >&2
  exit 2
fi

mkdir -p "$(dirname "$LOG_FILE")"

echo "[monitor-start] $(date '+%F %T %Z')" >>"$LOG_FILE"
echo "[monitor-config] pid_file=${PID_FILE:-<none>} pid=${PID:-<none>} interval=${INTERVAL}s run_root=$RUN_ROOT" >>"$LOG_FILE"

while true; do
  NOW="$(date '+%F %T %Z')"
  CUR_PID="$PID"
  if [[ -n "$PID_FILE" && -s "$PID_FILE" ]]; then
    CUR_PID="$(tr -d ' \t\r\n' < "$PID_FILE")"
  fi

  echo "[$NOW] ----" >>"$LOG_FILE"
  echo "[$NOW] pid=$CUR_PID" >>"$LOG_FILE"

  if [[ -n "$CUR_PID" ]] && kill -0 "$CUR_PID" 2>/dev/null; then
    echo "[$NOW] alive=yes" >>"$LOG_FILE"
    ps -p "$CUR_PID" -o pid,ppid,pcpu,pmem,etime,stat,cmd --no-headers >>"$LOG_FILE" 2>&1 || true
    pgrep -af "raw_performance_eval_one.py|10_stage_pas_match_quantify.R|20_stage_apa_and_te.R|30_stage_export_results.R" >>"$LOG_FILE" 2>&1 || true
  else
    echo "[$NOW] alive=no" >>"$LOG_FILE"
  fi

  for pair in "${PAIRS[@]}"; do
    done_file="$RUN_ROOT/raw_performance/state/eval/scape/$SAMPLE/$pair.done"
    fail_file="$RUN_ROOT/raw_performance/failures/eval/scape/$SAMPLE/$pair/failure.tsv"
    done_status="no"
    fail_status="no"
    [[ -f "$done_file" ]] && done_status="yes"
    [[ -f "$fail_file" ]] && fail_status="yes"
    echo "[$NOW] pair=$pair done=$done_status failure=$fail_status" >>"$LOG_FILE"
  done

  sleep "$INTERVAL"
done
