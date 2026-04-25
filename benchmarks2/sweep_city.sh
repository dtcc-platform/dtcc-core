#!/usr/bin/env bash
# Parameter sweep harness for supported Swedish city centers.
#
# Usage:
#   ./sweep_city.sh stockholm              # 300 s per case (default)
#   ./sweep_city.sh stockholm 420          # override per-case timeout
#
# Supported cities:
#   lund stockholm gothenburg malmo uppsala linkoping orebro vasteras helsingborg norrkoping
#
# Requires DTCC_LIDAR_URL and DTCC_GPKG_URL in the calling shell.
#
# Example:
#   export DTCC_LIDAR_URL=http://13.60.69.202:8001
#   export DTCC_GPKG_URL=http://13.60.69.202:8001

set -u

exec 3>&1 4>&2

CITY="${1:-}"
TIMEOUT_SECS="${2:-300}"

case "$CITY" in
  lund|stockholm|gothenburg|malmo|uppsala|linkoping|orebro|vasteras|helsingborg|norrkoping)
    ;;
  "")
    echo "usage: $0 <city> [timeout-seconds]" >&2
    exit 2
    ;;
  *)
    echo "unknown city: $CITY" >&2
    echo "valid cities: lund, stockholm, gothenburg, malmo, uppsala, linkoping, orebro, vasteras, helsingborg, norrkoping" >&2
    exit 2
    ;;
esac

LOG_DIR="./sweep_${CITY}_logs"
mkdir -p "$LOG_DIR"

: "${DTCC_LIDAR_URL:?DTCC_LIDAR_URL must be exported}"
: "${DTCC_GPKG_URL:?DTCC_GPKG_URL must be exported}"

CASES=(
  raster_0.5 raster_1.0 raster_2.0 raster_5.0 raster_10.0
  detail_0.25 detail_0.5 detail_1.0 detail_2.0
  area_1.0 area_5.0 area_15.0 area_25.0 area_100.0
  bbox_50 bbox_100 bbox_200 bbox_350 bbox_500
  max_1 max_2 max_5 max_10 max_20
)

run_with_escalation() {
  local secs=$1; shift
  local log=$1; shift
  local stage_file start end elapsed rc signal
  stage_file=$(mktemp "${TMPDIR:-/tmp}/sweep_${CITY}_stage.XXXXXX")
  start=$(date +%s)

  : >"$log"
  "$@" > >(tee -a "$log" >&3) 2> >(tee -a "$log" >&4) &
  local pid=$!

  ( sleep "$secs"          && { echo SIGINT  >>"$stage_file"; kill -INT  "$pid" 2>/dev/null; } ) </dev/null >/dev/null 2>&1 &
  local w1=$!
  ( sleep "$((secs + 10))" && { echo SIGTERM >>"$stage_file"; kill -TERM "$pid" 2>/dev/null; } ) </dev/null >/dev/null 2>&1 &
  local w2=$!
  ( sleep "$((secs + 20))" && { echo SIGKILL >>"$stage_file"; kill -KILL "$pid" 2>/dev/null; } ) </dev/null >/dev/null 2>&1 &
  local w3=$!

  wait "$pid" 2>/dev/null
  rc=$?
  end=$(date +%s)
  elapsed=$((end - start))

  sleep 0.5
  kill "$w1" "$w2" "$w3" 2>/dev/null
  wait "$w1" "$w2" "$w3" 2>/dev/null

  if [ -s "$stage_file" ]; then
    signal=$(tail -n 1 "$stage_file")
  else
    signal="none"
  fi
  rm -f "$stage_file"

  echo "signal_used=$signal rc=$rc elapsed=$elapsed"
}

declare -a RESULTS

classify_case() {
  local case_id="$1" log="$2" signal="$3" rc="$4" elapsed="$5"
  local line status verts faces bfaces
  line=$(grep "^SWEEP_RESULT " "$log" | tail -n 1 || true)

  if [ -n "$line" ]; then
    status=$(awk 'match($0, /status=[^ ]+/) { print substr($0, RSTART+7, RLENGTH-7) }' <<<"$line")
    verts=$(awk 'match($0, /vertices=[0-9]+/) { print substr($0, RSTART+9, RLENGTH-9) }' <<<"$line")
    faces=$(awk 'match($0, /faces=[0-9]+/) { print substr($0, RSTART+6, RLENGTH-6) }' <<<"$line")
    bfaces=$(awk 'match($0, /building_faces=[0-9]+/) { print substr($0, RSTART+15, RLENGTH-15) }' <<<"$line")
    case "$status" in
      ok)
        RESULTS+=("$case_id|OK|elapsed=${elapsed}s verts=${verts} faces=${faces} building_faces=${bfaces}")
        ;;
      ok_terrain_only)
        RESULTS+=("$case_id|OK_TERRAIN_ONLY|elapsed=${elapsed}s verts=${verts} faces=${faces} building_faces=0")
        ;;
      error)
        RESULTS+=("$case_id|FAIL|elapsed=${elapsed}s (status=error from driver)")
        ;;
      *)
        RESULTS+=("$case_id|FAIL|elapsed=${elapsed}s (unrecognized status=$status)")
        ;;
    esac
    return
  fi

  case "$signal" in
    SIGTERM|SIGKILL)
      RESULTS+=("$case_id|HANG|elapsed=${elapsed}s killed via $signal")
      ;;
    SIGINT|none)
      RESULTS+=("$case_id|FAIL|elapsed=${elapsed}s rc=$rc (no SWEEP_RESULT, died before escalation)")
      ;;
    *)
      RESULTS+=("$case_id|FAIL|elapsed=${elapsed}s rc=$rc (no SWEEP_RESULT, unknown signal=$signal)")
      ;;
  esac
}

repeat_char() {
  local char=$1 count=$2 out="" i
  for ((i = 0; i < count; i++)); do
    out+="$char"
  done
  printf '%s' "$out"
}

print_summary_rule() {
  local left=$1 mid=$2 right=$3
  printf '  %s' "$left"
  repeat_char "─" 16
  printf '%s' "$mid"
  repeat_char "─" 16
  printf '%s' "$mid"
  repeat_char "─" 10
  printf '%s' "$mid"
  repeat_char "─" 10
  printf '%s' "$mid"
  repeat_char "─" 10
  printf '%s' "$mid"
  repeat_char "─" 16
  printf '%s\n' "$right"
}

extract_summary_field() {
  local details=$1 field=$2
  if [[ $details =~ (^|[[:space:]])${field}=([^[:space:]]+) ]]; then
    printf '%s' "${BASH_REMATCH[2]}"
  else
    printf '-'
  fi
}

print_boxed_summary_table() {
  local row case_id rest classification details elapsed verts faces bfaces

  echo
  echo "Sweep results table:"
  print_summary_rule "╭" "┬" "╮"
  printf '  │ %-14s │ %-14s │ %8s │ %8s │ %8s │ %14s │\n' \
    "Case" "Result" "Elapsed" "Vertices" "Faces" "Building Faces"
  print_summary_rule "├" "┼" "┤"

  for row in "${RESULTS[@]}"; do
    case_id=${row%%|*}
    rest=${row#*|}
    classification=${rest%%|*}
    details=${rest#*|}
    elapsed=$(extract_summary_field "$details" "elapsed")
    verts=$(extract_summary_field "$details" "verts")
    faces=$(extract_summary_field "$details" "faces")
    bfaces=$(extract_summary_field "$details" "building_faces")
    printf '  │ %-14s │ %-14s │ %8s │ %8s │ %8s │ %14s │\n' \
      "$case_id" "$classification" "$elapsed" "$verts" "$faces" "$bfaces"
  done

  print_summary_rule "╰" "┴" "╯"
}

print_failed_summary_table() {
  local row case_id rest classification details elapsed verts faces bfaces found=0

  for row in "${RESULTS[@]}"; do
    rest=${row#*|}
    classification=${rest%%|*}
    if [ "$classification" != "OK" ]; then
      found=1
      break
    fi
  done

  echo
  echo "Failed tests table:"
  if [ "$found" -eq 0 ]; then
    echo "  No failed tests."
    return
  fi

  print_summary_rule "╭" "┬" "╮"
  printf '  │ %-14s │ %-14s │ %8s │ %8s │ %8s │ %14s │\n' \
    "Case" "Result" "Elapsed" "Vertices" "Faces" "Building Faces"
  print_summary_rule "├" "┼" "┤"

  for row in "${RESULTS[@]}"; do
    case_id=${row%%|*}
    rest=${row#*|}
    classification=${rest%%|*}
    if [ "$classification" = "OK" ]; then
      continue
    fi
    details=${rest#*|}
    elapsed=$(extract_summary_field "$details" "elapsed")
    verts=$(extract_summary_field "$details" "verts")
    faces=$(extract_summary_field "$details" "faces")
    bfaces=$(extract_summary_field "$details" "building_faces")
    printf '  │ %-14s │ %-14s │ %8s │ %8s │ %8s │ %14s │\n' \
      "$case_id" "$classification" "$elapsed" "$verts" "$faces" "$bfaces"
  done

  print_summary_rule "╰" "┴" "╯"
}

for case_id in "${CASES[@]}"; do
  log="${LOG_DIR}/${case_id}.log"
  echo "===> ${CITY}/${case_id} (timeout=${TIMEOUT_SECS}s, log=${log})"
  watchdog_out=$(run_with_escalation "$TIMEOUT_SECS" "$log" \
    uv run --frozen --python 3.11 python sweep_city.py "$CITY" "$case_id")
  signal=$(awk 'match($0, /signal_used=[A-Za-z]+/) { print substr($0, RSTART+12, RLENGTH-12) }' <<<"$watchdog_out")
  rc=$(awk 'match($0, /rc=-?[0-9]+/) { print substr($0, RSTART+3, RLENGTH-3) }' <<<"$watchdog_out")
  elapsed=$(awk 'match($0, /elapsed=[0-9]+/) { print substr($0, RSTART+8, RLENGTH-8) }' <<<"$watchdog_out")
  classify_case "$case_id" "$log" "${signal:-none}" "${rc:-0}" "${elapsed:-0}"
done

echo
echo "===================== ${CITY} SWEEP SUMMARY ====================="
prev_axis=""
for row in "${RESULTS[@]}"; do
  case_id=${row%%|*}
  rest=${row#*|}
  classification=${rest%%|*}
  details=${rest#*|}
  axis=${case_id%%_*}
  if [ "$axis" != "$prev_axis" ]; then
    case "$axis" in
      raster) header="raster_cell_size sweep:" ;;
      detail) header="min_building_detail sweep:" ;;
      area)   header="min_building_area sweep:" ;;
      bbox)   header="bbox_size_m sweep:" ;;
      max)    header="max_mesh_size sweep:" ;;
      *)      header="$axis sweep:" ;;
    esac
    [ -n "$prev_axis" ] && echo
    echo "$header"
    prev_axis="$axis"
  fi
  printf '  %-20s %-18s (%s)\n' "$case_id:" "$classification" "$details"
done
echo
print_boxed_summary_table
print_failed_summary_table
echo
echo "=========================================================="
echo "logs in ${LOG_DIR}/"
