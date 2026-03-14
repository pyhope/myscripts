#!/usr/bin/bash
# watch_squeue_idle.sh

set -euo pipefail

USER_ID="yp0007"
SQUEUE_FMT="%.18i %.9P %.8j %.2t %.10M %.6D %Z %L"
WARN_AFTER=${1:-180}
KILL_AFTER=${2:-900}
KILLED_DIRS_FILE="$(pwd)/rerun_list.txt"
LOG_FILE="$(pwd)/sq.log"

STOPCAR_THRESHOLD_SEC=$((20 * 60))
RESTART_SKIP_THRESHOLD_SEC=$((3 * 3600))   # <= 3 hours: do not restart
RESTART_ARG5_THRESHOLD_SEC=$((8 * 3600))   # > 3 and <= 8 hours: use argument 5

# Modes:
# monitor -> CANCEL_ON_STALE=0, RESTART_ON_KILL=0
# cancel  -> CANCEL_ON_STALE=1, RESTART_ON_KILL=0
# restart -> CANCEL_ON_STALE=1, RESTART_ON_KILL=1
CANCEL_ON_STALE=0
RESTART_ON_KILL=0

log_msg() {
  echo "$*" >> "$LOG_FILE"
}

parse_timelimit_to_seconds() {
  local s="${1:-}"
  s="${s// /}"

  if [[ -z "$s" || "$s" == "N/A" || "$s" == "UNLIMITED" || "$s" == "INVALID" ]]; then
    return 1
  fi

  local days=0 hours=0 mins=0 secs=0 right a b c
  if [[ "$s" == *-* ]]; then
    days="${s%%-*}"
    right="${s#*-}"
  else
    right="$s"
  fi

  IFS=':' read -r a b c <<< "$right"

  if [[ -n "${c:-}" ]]; then
    hours="$a"; mins="$b"; secs="$c"
  elif [[ -n "${b:-}" ]]; then
    mins="$a"; secs="$b"
  else
    secs="$a"
  fi

  [[ "$days"  =~ ^[0-9]+$ ]] || return 1
  [[ "$hours" =~ ^[0-9]+$ ]] || return 1
  [[ "$mins"  =~ ^[0-9]+$ ]] || return 1
  [[ "$secs"  =~ ^[0-9]+$ ]] || return 1

  echo $((10#$days*86400 + 10#$hours*3600 + 10#$mins*60 + 10#$secs))
}

print_status() {
  local mode
  if (( CANCEL_ON_STALE == 0 )); then
    mode="monitor_only"
  elif (( RESTART_ON_KILL == 0 )); then
    mode="monitor_cancel"
  else
    mode="monitor_cancel_restart"
  fi

  echo "STATUS: MODE=$mode"
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
    local content='LSTOP = .TRUE.'

    if [[ -f "$stopcar_path" ]] && grep -qxF "$content" "$stopcar_path"; then
      log_msg "INFO: JOBID=$jobid | WORK_DIR=$workdir | time_left=$left_str | STOPCAR already set"
      return 0
    fi

    if ( cd "$workdir" && echo "$content" > STOPCAR ); then
      log_msg "ACTION: JOBID=$jobid | WORK_DIR=$workdir | time_left=$left_str (<=threshold) | wrote STOPCAR"
    else
      log_msg "ERROR: JOBID=$jobid | WORK_DIR=$workdir | failed to write STOPCAR"
    fi
  fi
}

run_restart_script() {
  local jobid="$1"
  local workdir="$2"
  local left_str="$3"

  if (( RESTART_ON_KILL == 0 )); then
    return 0
  fi

  if ! command -v vml_restart >/dev/null 2>&1; then
    log_msg "ERROR: JOBID=$jobid | WORK_DIR=$workdir | vml_restart not found in PATH"
    return 1
  fi

  local left_sec restart_args=()

  if ! left_sec="$(parse_timelimit_to_seconds "$left_str")"; then
    log_msg "INFO: JOBID=$jobid | WORK_DIR=$workdir | time_left=$left_str | cannot parse time left, using default vml_restart"
  else
    if (( left_sec <= RESTART_SKIP_THRESHOLD_SEC )); then
      log_msg "INFO: JOBID=$jobid | WORK_DIR=$workdir | time_left=$left_str (<=3h) | skip vml_restart"
      return 0
    elif (( left_sec <= RESTART_ARG5_THRESHOLD_SEC )); then
      restart_args=(5)
      log_msg "ACTION: JOBID=$jobid | WORK_DIR=$workdir | time_left=$left_str (3h-8h) | starting vml_restart 5"
    else
      log_msg "ACTION: JOBID=$jobid | WORK_DIR=$workdir | time_left=$left_str (>8h) | starting vml_restart"
    fi
  fi

  if (
    cd "$workdir" && \
    vml_restart "${restart_args[@]}" >> "$LOG_FILE" 2>&1
  ); then
    if (( ${#restart_args[@]} > 0 )); then
      log_msg "ACTION: JOBID=$jobid | WORK_DIR=$workdir | executed vml_restart ${restart_args[*]} successfully"
    else
      log_msg "ACTION: JOBID=$jobid | WORK_DIR=$workdir | executed vml_restart successfully"
    fi
  else
    if (( ${#restart_args[@]} > 0 )); then
      log_msg "ERROR: JOBID=$jobid | WORK_DIR=$workdir | vml_restart ${restart_args[*]} failed"
    else
      log_msg "ERROR: JOBID=$jobid | WORK_DIR=$workdir | vml_restart failed"
    fi
    return 1
  fi
}

check_once() {
  local now_epoch
  local jobs
  now_epoch=$(date +%s)

  jobs=$(squeue -u "$USER_ID" --noheader --format="$SQUEUE_FMT")
  if [[ -z "$jobs" ]]; then
    log_msg "No jobs running for $USER_ID, exiting."
    exit 0
  fi

  log_msg "===== $(date '+%Y-%m-%d %H:%M:%S%z') ====="

  # jobid | state | workdir | timeleft
  echo "$jobs" | awk '{print $1"\t"$4"\t"$(NF-1)"\t"$NF}' \
  | while IFS=$'\t' read -r jobid state workdir timeleft; do
      [[ -z "${jobid:-}" || -z "${workdir:-}" ]] && continue

      if [[ "$state" != "R" ]]; then
        continue
      fi

      if [[ ! -d "$workdir" ]]; then
        log_msg "WARN: JOBID=$jobid | STATE=$state | WORK_DIR=$workdir | directory not found"
        continue
      fi

      maybe_write_stopcar "$jobid" "$workdir" "${timeleft:-}"

      newest_epoch=$(
        find "$workdir" -type f -printf '%T@\n' 2>/dev/null \
        | awk 'BEGIN{max=0} {if ($1>max) max=$1} END{printf "%d", max}'
      )

      diff=$(( now_epoch - newest_epoch ))

      if (( diff >= KILL_AFTER )); then
        local last_str idle_min idle_sec
        last_str=$(date -d "@$newest_epoch" '+%Y-%m-%d %H:%M:%S%z')
        idle_min=$(( diff / 60 ))
        idle_sec=$(( diff % 60 ))

        if (( CANCEL_ON_STALE == 0 )); then
          log_msg "STALE: JOBID=$jobid | WORK_DIR=$workdir | last_update=$last_str | idle=${idle_min}m${idle_sec}s | time_left=${timeleft:-N/A} | monitor-only, no scancel"
        else
          log_msg "STALE: scancel $jobid | WORK_DIR=$workdir | last_update=$last_str | idle=${idle_min}m${idle_sec}s | time_left=${timeleft:-N/A}"

          if scancel "$jobid"; then
            printf '%s\n' "$workdir" >> "$KILLED_DIRS_FILE"
            run_restart_script "$jobid" "$workdir" "${timeleft:-}"
          else
            log_msg "ERROR: failed to scancel $jobid"
          fi
        fi

      elif (( diff >= WARN_AFTER )); then
        local last_str idle_min idle_sec
        last_str=$(date -d "@$newest_epoch" '+%Y-%m-%d %H:%M:%S%z')
        idle_min=$(( diff / 60 ))
        idle_sec=$(( diff % 60 ))
        log_msg "WARN: JOBID=$jobid | WORK_DIR=$workdir | last_update=$last_str | idle=${idle_min}m${idle_sec}s | time_left=${timeleft:-N/A}"
      fi
    done
}

wait_with_commands() {
  local remaining="$1"
  local cmd

  while (( remaining > 0 )); do
    if read -r -t 1 cmd; then
      case "$cmd" in
        monitor)
          CANCEL_ON_STALE=0
          RESTART_ON_KILL=0
          echo "COMMAND: monitor-only mode ON"
          print_status
          ;;
        cancel)
          CANCEL_ON_STALE=1
          RESTART_ON_KILL=0
          echo "COMMAND: monitor+cancel mode ON"
          print_status
          ;;
        restart)
          CANCEL_ON_STALE=1
          RESTART_ON_KILL=1
          echo "COMMAND: monitor+cancel+restart mode ON"
          print_status
          ;;
        status)
          print_status
          ;;
        check)
          echo "COMMAND: immediate check"
          check_once
          ;;
        quit|exit)
          echo "COMMAND: exit"
          exit 0
          ;;
        "")
          ;;
        *)
          echo "Unknown command: $cmd"
          echo "Available commands: monitor | cancel | restart | status | check | quit"
          ;;
      esac
    fi
    remaining=$((remaining - 1))
  done
}

echo "Interactive commands enabled:"
echo "  monitor   -> monitor only, do not cancel, do not restart"
echo "  cancel    -> monitor and cancel stale jobs, but do not restart"
echo "  restart   -> monitor, cancel stale jobs, and conditionally run vml_restart"
echo "  status    -> show current mode"
echo "  check     -> run check immediately"
echo "  quit      -> exit script"
echo "Periodic check output and vml_restart output will be appended to: $LOG_FILE"

print_status

while :; do
  check_once
  wait_with_commands "$WARN_AFTER"
done
