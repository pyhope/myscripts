#!/usr/bin/bash
# watch_squeue_idle.sh

set -euo pipefail

USER_ID="yp0007"

SQUEUE_FMT="%.18i|%.2t|%Z|%L"

WARN_AFTER=${1:-'120'}
KILL_AFTER=${2:-'600'}
KILLED_DIRS_FILE="./rerun_list.txt"

STOPCAR_THRESHOLD_SEC=$((10 * 60))

parse_timelimit_to_seconds() {
  local s="${1:-}"
  s="${s// /}"

  local days=0 hours=0 mins=0 secs=0 left right
  if [[ "$s" == *-* ]]; then
    days="${s%%-*}"
    right="${s#*-}"
  else
    right="$s"
  fi

  IFS=':' read -r a b c <<< "$right"

  if [[ -n "${c:-}" ]]; then
    # a:b:c -> H:M:S
    hours="$a"; mins="$b"; secs="$c"
  elif [[ -n "${b:-}" ]]; then
    # a:b -> M:S
    mins="$a"; secs="$b"
  else
    # a -> S
    secs="$a"
  fi

  [[ "$days" =~ ^[0-9]+$ ]] || return 1
  [[ "$hours" =~ ^[0-9]+$ ]] || return 1
  [[ "$mins" =~ ^[0-9]+$ ]] || return 1
  [[ "$secs" =~ ^[0-9]+$ ]] || return 1

  echo $(( days*86400 + hours*3600 + mins*60 + secs ))
}

maybe_write_stopcar() {
  local jobid="$1"
  local workdir="$2"
  local left_str="$3"

  local left_sec
  if ! left_sec="$(parse_timelimit_to_seconds "$left_str")"; then
    return 0
  fi

  if (( left_sec <= STOPCAR_THRESHOLD_SEC )); then
    local stopcar_path="$workdir/STOPCAR"
    local content='LABORT = .TRUE.'

    if [[ -f "$stopcar_path" ]] && grep -qxF "$content" "$stopcar_path"; then
      echo "INFO: JOBID=$jobid | WORK_DIR=$workdir | time_left=$left_str | STOPCAR already set"
      return 0
    fi

    if ( cd "$workdir" && echo "$content" > STOPCAR ); then
      echo "ACTION: JOBID=$jobid | WORK_DIR=$workdir | time_left=$left_str (<=5m) | wrote STOPCAR: '$content'"
    else
      echo "ERROR: JOBID=$jobid | WORK_DIR=$workdir | failed to write STOPCAR" >&2
    fi
  fi
}

check_once() {
  local now_epoch
  now_epoch=$(date +%s)

  jobs=$(squeue -u "$USER_ID" --noheader --format="$SQUEUE_FMT")
  if [[ -z "$jobs" ]]; then
    echo "No jobs running for $USER_ID, exiting."
    exit 0
  fi

  # jobid|state|workdir|timeleft
  echo "$jobs" \
  | while IFS='|' read -r jobid state workdir timeleft; do
      [[ -z "${jobid:-}" || -z "${workdir:-}" ]] && continue

      # only check RUNNING jobs
      if [[ "$state" != "R" ]]; then
        # echo "INFO: JOBID=$jobid | STATE=$state | WORK_DIR=$workdir | skip (not running)"
        continue
      fi

      if [[ ! -d "$workdir" ]]; then
        echo "WARN: JOBID=$jobid | STATE=$state | WORK_DIR=$workdir | directory not found"
        continue
      fi

      maybe_write_stopcar "$jobid" "$workdir" "${timeleft:-}"

      newest_epoch=$(
        find "$workdir" -type f -printf '%T@\n' 2>/dev/null \
        | awk 'BEGIN{max=0} {if ($1>max) max=$1} END{printf "%d", max}'
      )

      diff=$(( now_epoch - newest_epoch ))

      if (( diff >= KILL_AFTER )); then
        last_str=$(date -d "@$newest_epoch" '+%Y-%m-%d %H:%M:%S%z')
        idle_min=$(( diff / 60 ))
        idle_sec=$(( diff % 60 ))
        echo "STALE: scancel $jobid | WORK_DIR=$workdir | last_update=$last_str | idle=${idle_min}m${idle_sec}s | time_left=${timeleft:-N/A}"

        if scancel "$jobid"; then
          printf '%s\n' "$workdir" >> "$KILLED_DIRS_FILE"
        else
          echo "ERROR: failed to scancel $jobid" >&2
        fi

      elif (( diff >= WARN_AFTER )); then
        last_str=$(date -d "@$newest_epoch" '+%Y-%m-%d %H:%M:%S%z')
        idle_min=$(( diff / 60 ))
        idle_sec=$(( diff % 60 ))
        echo "WARN: JOBID=$jobid | WORK_DIR=$workdir | last_update=$last_str | idle=${idle_min}m${idle_sec}s | time_left=${timeleft:-N/A}"
      fi
    done
}

while :; do
  echo "===== $(date '+%Y-%m-%d %H:%M:%S%z') ====="
  check_once
  sleep "$WARN_AFTER"
done
