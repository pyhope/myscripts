#!/usr/bin/bash
# watch_squeue_idle.sh

set -euo pipefail

USER_ID="yp0007"
SQUEUE_FMT="%.18i %.9P %.8j %.2t %.10M %.6D %Z"
WARN_AFTER=${1:-'60'}
KILL_AFTER=${2:-'600'}
KILLED_DIRS_FILE="./rerun_list.txt"

check_once() {
  local now_epoch
  now_epoch=$(date +%s)

  jobs=$(squeue -u "$USER_ID" --noheader --format="$SQUEUE_FMT")
  if [[ -z "$jobs" ]]; then
    echo "No jobs running for $USER_ID, exiting."
    exit 0
  fi

  # jobid | state | workdir
  echo "$jobs" | awk '{print $1"\t"$4"\t"$NF}' \
  | while IFS=$'\t' read -r jobid state workdir; do
      [[ -z "${jobid:-}" || -z "${workdir:-}" ]] && continue

      # only check RUNNING jobs
      if [[ "$state" != "R" ]]; then
        echo "INFO: JOBID=$jobid | STATE=$state | WORK_DIR=$workdir | skip (not running)"
        continue
      fi

      if [[ ! -d "$workdir" ]]; then
        echo "WARN: JOBID=$jobid | STATE=$state | WORK_DIR=$workdir | directory not found"
        continue
      fi

      newest_epoch=$(
        find "$workdir" -type f -printf '%T@\n' 2>/dev/null \
        | awk 'BEGIN{max=0} {if ($1>max) max=$1} END{printf "%d", max}'
      )

      diff=$(( now_epoch - newest_epoch ))

      if (( diff >= KILL_AFTER )); then
        last_str=$(date -d "@$newest_epoch" '+%Y-%m-%d %H:%M:%S%z')
        idle_min=$(( diff / 60 ))
        idle_sec=$(( diff % 60 ))
        echo "STALE: scancel $jobid | WORK_DIR=$workdir | last_update=$last_str | idle=${idle_min}m${idle_sec}s"

        if scancel "$jobid"; then
          printf '%s\n' "$workdir" >> "$KILLED_DIRS_FILE"
        else
          echo "ERROR: failed to scancel $jobid" >&2
        fi

      elif (( diff >= WARN_AFTER )); then
        last_str=$(date -d "@$newest_epoch" '+%Y-%m-%d %H:%M:%S%z')
        idle_min=$(( diff / 60 ))
        idle_sec=$(( diff % 60 ))
        echo "WARN: JOBID=$jobid | WORK_DIR=$workdir | last_update=$last_str | idle=${idle_min}m${idle_sec}s"
      fi
    done
}

while :; do
  echo "===== $(date '+%Y-%m-%d %H:%M:%S%z') ====="
  check_once
  sleep "$WARN_AFTER"
done
