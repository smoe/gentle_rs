#!/usr/bin/env bash
# Linux/X11 semantic-GUI smoke example. This script contains no synthetic input
# API: xdotool converts read-only semantic rectangles into normal OS events.
set -euo pipefail

for tool in jq xdotool; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "SKIP: semantic GUI smoke requires $tool" >&2
    exit 0
  fi
done
if [[ -z "${DISPLAY:-}" ]]; then
  echo "SKIP: semantic GUI smoke requires DISPLAY (run under xvfb-run)" >&2
  exit 0
fi
if [[ $# -lt 2 ]]; then
  echo "usage: $0 GENTLE_BINARY TP73_PROJECT [SNAPSHOT.json]" >&2
  exit 2
fi

binary=$1
project=$2
snapshot=${3:-/tmp/gentle-gui-semantic.json}
repeat_scope=${GENTLE_REPEAT_SCOPE:-}
array_scope=${GENTLE_ARRAY_SCOPE:-}
sequence_scope=${GENTLE_SEQUENCE_SCOPE:-}
rm -f "$snapshot"
GENTLE_GUI_TEST_SNAPSHOT="$snapshot" "$binary" --project "$project" &
gentle_pid=$!
trap 'kill "$gentle_pid" 2>/dev/null || true' EXIT

click_item() {
  local semantic_id=$1 scope=${2:-}
  local selector='.items[] | select(.semantic_id == $id)'
  if [[ -n "$scope" ]]; then
    selector+=' | select(.subject_scope == $scope)'
  fi
  local point
  point=$(jq -r --arg id "$semantic_id" --arg scope "$scope" \
    "$selector | select(.state.visible and .state.enabled) | [((.rect_logical_points.min_x + .rect_logical_points.max_x) / 2 * .pixels_per_point | floor), ((.rect_logical_points.min_y + .rect_logical_points.max_y) / 2 * .pixels_per_point | floor)] | @tsv" \
    "$snapshot" | head -1)
  [[ -n "$point" ]] || { echo "FAIL: missing enabled semantic item $semantic_id $scope" >&2; exit 1; }
  read -r x y <<<"$point"
  xdotool mousemove --sync "$x" "$y" click 1
}

double_click_item() {
  local semantic_id=$1 scope=${2:-}
  local selector='.items[] | select(.semantic_id == $id)'
  if [[ -n "$scope" ]]; then
    selector+=' | select(.subject_scope == $scope)'
  fi
  local point
  point=$(jq -r --arg id "$semantic_id" --arg scope "$scope" \
    "$selector | select(.state.visible and .state.enabled) | [((.rect_logical_points.min_x + .rect_logical_points.max_x) / 2 * .pixels_per_point | floor), ((.rect_logical_points.min_y + .rect_logical_points.max_y) / 2 * .pixels_per_point | floor)] | @tsv" \
    "$snapshot" | head -1)
  [[ -n "$point" ]] || { echo "FAIL: missing enabled semantic item $semantic_id $scope" >&2; exit 1; }
  read -r x y <<<"$point"
  xdotool mousemove --sync "$x" "$y" click --repeat 2 --delay 100 1
}

deadline=$((SECONDS + 30))
while (( SECONDS < deadline )); do
  if [[ -s "$snapshot" ]] && jq -e '
    .settled
    and (.generation > 0)
    and any(.items[]; .semantic_id == "main.project.sequence.open" and .state.visible and .state.enabled)
  ' "$snapshot" >/dev/null; then
    break
  fi
  sleep 0.1
done
if ! [[ -s "$snapshot" ]] || ! jq -e '
  .settled
  and any(.items[]; .semantic_id == "main.project.sequence.open" and .state.visible and .state.enabled)
' "$snapshot" >/dev/null; then
  echo "FAIL: semantic GUI snapshot did not settle with an enabled project sequence control" >&2
  exit 1
fi

double_click_item main.project.sequence.open "$sequence_scope"

deadline=$((SECONDS + 30))
while (( SECONDS < deadline )); do
  if [[ -s "$snapshot" ]] && jq -e '
    .settled
    and any(.items[]; .semantic_id == "window.dna_viewer")
    and any(.items[]; .semantic_id == "dna.splitter.info_width")
  ' "$snapshot" >/dev/null; then
    break
  fi
  sleep 0.1
done
if ! [[ -s "$snapshot" ]] || ! jq -e '
  .settled
  and any(.items[]; .semantic_id == "window.dna_viewer")
  and any(.items[]; .semantic_id == "dna.splitter.info_width")
' "$snapshot" >/dev/null; then
  echo "FAIL: semantic GUI snapshot did not settle with the DNA viewer and splitter ready" >&2
  exit 1
fi

jq -e '.items[] | select(.semantic_id == "window.dna_viewer")' "$snapshot" >/dev/null
click_item dna.splitter.info_width
if [[ -n "$repeat_scope" ]]; then click_item dna.feature_tree.row "$repeat_scope"; fi
if [[ -n "$array_scope" ]]; then click_item dna.feature_tree.row "$array_scope"; fi

if command -v scrot >/dev/null 2>&1; then
  scrot "${snapshot%.json}.png"
fi
if [[ -n "${GENTLE_VERIFY_COMMAND:-}" ]]; then
  bash -lc "$GENTLE_VERIFY_COMMAND"
else
  echo "NOTE: set GENTLE_VERIFY_COMMAND to assert the typed project/report result" >&2
fi
jq '{schema, generation, settled, items}' "$snapshot"
