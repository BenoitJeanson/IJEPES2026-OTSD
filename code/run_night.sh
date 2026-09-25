#!/bin/bash
# Regenerate both campaigns from the current code, in sequence by construction --
# no polling, no pgrep.
#
# `set -e` matters: a campaign that dies mid-stage used to leave this script printing
# "all done" and exiting 0 over six cells that never ran. Each stage exits nonzero when
# any cell fails, and this stops on the first one.
#
# results/ is cleared first, on purpose. A campaign skips cells already on disk, so
# without this it would silently keep records produced by superseded code -- which is
# exactly the trap that made IEEE-118 report 79 LP solves for code that now takes 101.
# The old cells are moved aside rather than deleted, in case a number is needed before
# the new run lands.
set -euo pipefail
cd /Users/benoitjeanson/vsCode/TUD/IJEPES2026-OTSD/code
ts() { date +%Y%m%d_%H%M%S; }
STAMP=$(ts)

ARCHIVE="/tmp/otsd_results_superseded_$STAMP"
if compgen -G "results/*/result.json" > /dev/null; then
  mkdir -p "$ARCHIVE"
  mv results/*/ "$ARCHIVE"/ 2>/dev/null || true
  echo "=== $STAMP archived $(ls "$ARCHIVE" | wc -l | tr -d ' ') superseded cells to $ARCHIVE ===" >> logs/stages.log
fi

stage() {          # stage <backend>
  local backend=$1 log
  log="logs/$(ts)_${backend}_campaign.log"
  echo "=== $backend start $(date) -> $log ===" >> logs/stages.log
  if julia --project=. experiments/run_ablation.jl "$backend" > "$log" 2>&1; then
    echo "=== $backend done $(date) ===" >> logs/stages.log
  else
    local code=$?
    echo "=== $backend FAILED ($code) $(date); see $log ===" >> logs/stages.log
    return $code
  fi
}

stage gurobi
stage scip
echo "=== all done $(date) ===" >> logs/stages.log
