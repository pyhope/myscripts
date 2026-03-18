#!/usr/bin/bash
# vasp-watchdog.sh

set -euo pipefail

USER_ID=$(id -u -n)
SQUEUE_FMT="%.18i %.9P %.8j %.2t %.10M %.6D %Z %L"

# Default thresholds (seconds)
WARN_AFTER=$((5 * 60))
KILL_AFTER=$((30 * 60))

KILLED_DIRS_FILE="$(pwd)/rerun_list.txt"
LOG_FILE="$(pwd)/sq.log"

STOPCAR_THRESHOLD_SEC=$((20 * 60))
RESTART_SKIP_THRESHOLD_SEC=$((3 * 3600))    # < 3 hours: do not restart
CLEAN_RESTART_THRESHOLD_SEC=$((21 * 3600))  # > 21 hours: use vml_clean_restart
RESTART_EXTRA_SEC=$((1 * 3600))             # restart with remaining time + 1 h

# Modes:
# monitor -> CANCEL_ON_STALE=0, RESTART_ON_KILL=0, RECOVER_DISAPPEARED=0
# cancel  -> CANCEL_ON_STALE=1, RESTART_ON_KILL=0, RECOVER_DISAPPEARED=0
# restart -> CANCEL_ON_STALE=1, RESTART_ON_KILL=1, RECOVER_DISAPPEARED=0
# recover -> CANCEL_ON_STALE=1, RESTART_ON_KILL=1, RECOVER_DISAPPEARED=1
CANCEL_ON_STALE=0
RESTART_ON_KILL=0
RECOVER_DISAPPEARED=0

declare -A PREV_JOBID_BY_WORKDIR=()
declare -A PREV_TIMELEFT_BY_WORKDIR=()
declare -A CURR_JOBID_BY_WORKDIR=()
declare -A CURR_TIMELEFT_BY_WORKDIR=()

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

calc_restart_hours() {
  local left_sec="$1"
  local total_sec hours remainder

  total_sec=$((left_sec + RESTART_EXTRA_SEC))
  hours=$((total_sec / 3600))
  remainder=$((total_sec % 3600))

  if (( remainder >= 1800 )); then
    hours=$((hours + 1))
  fi

  (( hours < 1 )) && hours=1
  echo "$hours"
}

print_status() {
  local mode warn_min kill_min

  if (( CANCEL_ON_STALE == 0 )); then
    mode="monitor_only"
  elif (( RESTART_ON_KILL == 0 )); then
    mode="monitor_cancel"
  elif (( RECOVER_DISAPPEARED == 0 )); then
    mode="monitor_cancel_restart"
  else
    mode="monitor_cancel_restart_recover"
  fi

  warn_min=$((WARN_AFTER / 60))
  kill_min=$((KILL_AFTER / 60))

  echo "STATUS: MODE=$mode | WARN_AFTER=${warn_min} min | KILL_AFTER=${kill_min} min"
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

run_logged_command() {
  local jobid="$1"
  local workdir="$2"
  local cmd_name="$3"
  shift 3
  local cmd=( "$cmd_name" "$@" )

  if ! command -v "$cmd_name" >/dev/null 2>&1; then
    log_msg "ERROR: JOBID=$jobid | WORK_DIR=$workdir | $cmd_name not found in PATH"
    return 1
  fi

  if (
    cd "$workdir" && \
    "${cmd[@]}" >> "$LOG_FILE" 2>&1
  ); then
    log_msg "ACTION: JOBID=$jobid | WORK_DIR=$workdir | executed ${cmd[*]} successfully"
  else
    log_msg "ERROR: JOBID=$jobid | WORK_DIR=$workdir | ${cmd[*]} failed"
    return 1
  fi
}

run_restart_logic() {
  local jobid="$1"
  local workdir="$2"
  local left_str="$3"
  local left_sec restart_hours

  if ! left_sec="$(parse_timelimit_to_seconds "$left_str")"; then
    log_msg "INFO: JOBID=$jobid | WORK_DIR=$workdir | time_left=$left_str | cannot parse time left, using default vml_restart"
    run_logged_command "$jobid" "$workdir" vml_restart
    return $?
  fi

  if (( left_sec < RESTART_SKIP_THRESHOLD_SEC )); then
    log_msg "INFO: JOBID=$jobid | WORK_DIR=$workdir | time_left=$left_str (<3h) | skip restart"

  elif (( left_sec > CLEAN_RESTART_THRESHOLD_SEC )); then
    log_msg "ACTION: JOBID=$jobid | WORK_DIR=$workdir | time_left=$left_str (>21h) | starting vml_clean_restart"
    run_logged_command "$jobid" "$workdir" vml_clean_restart

  else
    restart_hours="$(calc_restart_hours "$left_sec")"
    log_msg "ACTION: JOBID=$jobid | WORK_DIR=$workdir | time_left=$left_str (3h-21h, +1h then round to nearest hour) | starting vml_restart $restart_hours"
    run_logged_command "$jobid" "$workdir" vml_restart "$restart_hours"
  fi
}

run_restart_script() {
  local jobid="$1"
  local workdir="$2"
  local left_str="$3"

  if (( RESTART_ON_KILL == 0 )); then
    return 0
  fi

  run_restart_logic "$jobid" "$workdir" "$left_str"
}

handle_disappeared_workdirs() {
  local workdir prev_jobid prev_timeleft prev_left_sec

  if (( RECOVER_DISAPPEARED == 0 )); then
    return 0
  fi

  for workdir in "${!PREV_JOBID_BY_WORKDIR[@]}"; do
    if [[ -n "${CURR_JOBID_BY_WORKDIR[$workdir]+x}" ]]; then
      continue
    fi

    prev_jobid="${PREV_JOBID_BY_WORKDIR[$workdir]}"
    prev_timeleft="${PREV_TIMELEFT_BY_WORKDIR[$workdir]}"

    if ! prev_left_sec="$(parse_timelimit_to_seconds "$prev_timeleft")"; then
      log_msg "INFO: PREV_JOBID=$prev_jobid | WORK_DIR=$workdir | prev_time_left=$prev_timeleft | disappeared from squeue but previous time left is unparseable, skip restart logic"
      continue
    fi

    if (( prev_left_sec < RESTART_SKIP_THRESHOLD_SEC )); then
      log_msg "INFO: PREV_JOBID=$prev_jobid | WORK_DIR=$workdir | prev_time_left=$prev_timeleft | disappeared from squeue but previous time left <3h, skip restart logic"
      continue
    fi

    if [[ ! -d "$workdir" ]]; then
      log_msg "WARN: PREV_JOBID=$prev_jobid | WORK_DIR=$workdir | prev_time_left=$prev_timeleft | disappeared from squeue and directory not found"
      continue
    fi

    log_msg "DISAPPEARED: PREV_JOBID=$prev_jobid | WORK_DIR=$workdir | prev_time_left=$prev_timeleft | missing from current squeue output, running restart logic"
    run_restart_logic "$prev_jobid" "$workdir" "$prev_timeleft"
  done
}

sync_previous_snapshot() {
  local workdir
  PREV_JOBID_BY_WORKDIR=()
  PREV_TIMELEFT_BY_WORKDIR=()

  for workdir in "${!CURR_JOBID_BY_WORKDIR[@]}"; do
    PREV_JOBID_BY_WORKDIR["$workdir"]="${CURR_JOBID_BY_WORKDIR[$workdir]}"
    PREV_TIMELEFT_BY_WORKDIR["$workdir"]="${CURR_TIMELEFT_BY_WORKDIR[$workdir]}"
  done
}

check_once() {
  local now_epoch
  local jobs
  now_epoch=$(date +%s)

  jobs=$(squeue -u "$USER_ID" --noheader --format="$SQUEUE_FMT" || true)

  CURR_JOBID_BY_WORKDIR=()
  CURR_TIMELEFT_BY_WORKDIR=()

  log_msg "===== $(date '+%Y-%m-%d %H:%M:%S%z') ====="

  if [[ -n "$jobs" ]]; then
    while IFS=$'\t' read -r jobid state workdir timeleft; do
      [[ -z "${jobid:-}" || -z "${workdir:-}" ]] && continue

      CURR_JOBID_BY_WORKDIR["$workdir"]="$jobid"
      CURR_TIMELEFT_BY_WORKDIR["$workdir"]="$timeleft"

      if [[ "$state" != "R" ]]; then
        continue
      fi

      if [[ ! -d "$workdir" ]]; then
        log_msg "WARN: JOBID=$jobid | STATE=$state | WORK_DIR=$workdir | directory not found"
        continue
      fi

      maybe_write_stopcar "$jobid" "$workdir" "${timeleft:-}"

      local newest_epoch diff last_str idle_min idle_sec
      newest_epoch=$(
        find "$workdir" -type f -printf '%T@\n' 2>/dev/null \
        | awk 'BEGIN{max=0} {if ($1>max) max=$1} END{printf "%d", max}'
      )

      diff=$(( now_epoch - newest_epoch ))

      if (( diff >= KILL_AFTER )); then
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
        last_str=$(date -d "@$newest_epoch" '+%Y-%m-%d %H:%M:%S%z')
        idle_min=$(( diff / 60 ))
        idle_sec=$(( diff % 60 ))
        log_msg "WARN: JOBID=$jobid | WORK_DIR=$workdir | last_update=$last_str | idle=${idle_min}m${idle_sec}s | time_left=${timeleft:-N/A}"
      fi
    done < <(echo "$jobs" | awk '{print $1"\t"$4"\t"$(NF-1)"\t"$NF}')
  fi

  handle_disappeared_workdirs
  sync_previous_snapshot

  if [[ -z "$jobs" ]]; then
    log_msg "No jobs running for $USER_ID, exiting."
    exit 0
  fi
}

set_warn_minutes() {
  local m="$1"
  [[ "$m" =~ ^[0-9]+$ ]] || { echo "ERROR: warn minutes must be a positive integer"; return 1; }
  (( m > 0 )) || { echo "ERROR: warn minutes must be > 0"; return 1; }
  WARN_AFTER=$((m * 60))
  echo "COMMAND: WARN_AFTER set to ${m} min"
  print_status
}

set_kill_minutes() {
  local m="$1"
  [[ "$m" =~ ^[0-9]+$ ]] || { echo "ERROR: kill minutes must be a positive integer"; return 1; }
  (( m > 0 )) || { echo "ERROR: kill minutes must be > 0"; return 1; }
  KILL_AFTER=$((m * 60))
  echo "COMMAND: KILL_AFTER set to ${m} min"
  print_status
}

wait_with_commands() {
  local remaining="$1"
  local cmd arg1

  while (( remaining > 0 )); do
    if read -r -t 1 cmd arg1; then
      case "$cmd" in
        monitor)
          CANCEL_ON_STALE=0
          RESTART_ON_KILL=0
          RECOVER_DISAPPEARED=0
          echo "COMMAND: monitor-only mode ON"
          print_status
          ;;
        cancel)
          CANCEL_ON_STALE=1
          RESTART_ON_KILL=0
          RECOVER_DISAPPEARED=0
          echo "COMMAND: monitor+cancel mode ON"
          print_status
          ;;
        restart)
          CANCEL_ON_STALE=1
          RESTART_ON_KILL=1
          RECOVER_DISAPPEARED=0
          echo "COMMAND: monitor+cancel+restart mode ON"
          print_status
          ;;
        recover)
          CANCEL_ON_STALE=1
          RESTART_ON_KILL=1
          RECOVER_DISAPPEARED=1
          echo "COMMAND: monitor+cancel+restart+recover mode ON"
          print_status
          ;;
        status)
          print_status
          ;;
        check)
          echo "COMMAND: immediate check"
          check_once
          ;;
        warn)
          if [[ -z "${arg1:-}" ]]; then
            echo "Usage: warn <minutes>"
          else
            set_warn_minutes "$arg1"
          fi
          ;;
        kill)
          if [[ -z "${arg1:-}" ]]; then
            echo "Usage: kill <minutes>"
          else
            set_kill_minutes "$arg1"
          fi
          ;;
        quit|exit)
          echo "COMMAND: exit"
          exit 0
          ;;
        "")
          ;;
        *)
          echo "Unknown command: $cmd"
          ;;
      esac
    fi
    remaining=$((remaining - 1))
  done
}

echo "Interactive commands enabled:"
echo "  monitor         -> monitor only, do not cancel, do not restart (default)"
echo "  cancel          -> monitor and cancel stale jobs, but do not restart"
echo "  restart         -> monitor, cancel stale jobs, and conditionally run restart script"
echo "  recover         -> restart mode + recover disappeared workdirs from previous round"
echo "  warn <minutes>  -> set WARN_AFTER in minutes"
echo "  kill <minutes>  -> set KILL_AFTER in minutes"
echo "  status          -> show current mode and thresholds"
echo "  check           -> run check immediately"
echo "  quit            -> exit script"
echo "Periodic check output and restart-script output will be appended to: $LOG_FILE"

print_status

while :; do
  check_once
  wait_with_commands "$WARN_AFTER"
done
