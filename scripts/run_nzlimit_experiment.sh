#!/usr/bin/env bash
# run_nzlimit_experiment.sh
#
# Runs a 5-condition experiment for perRoundNzCutLimitFactor tuning on v0.3.3.
# Conditions run sequentially; each uses --parallel to exploit all available cores.
#
# Usage: ./scripts/run_nzlimit_experiment.sh [--bin BINARY] [--outdir BASEDIR]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_EXP="$SCRIPT_DIR/run_experiments.sh"

BIN="/usr/bin/mipster"
INSTANCES="/home/haroldo/inst/miplib/2017+spp/"
OUTBASE="/home/haroldo/experiments/cbc/v0.3.3_nzlimit"
PARALLEL=125
TIMELIMIT=10800
OVERTIME_GRACE=600

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bin)     BIN="$2";     shift 2 ;;
    --outdir)  OUTBASE="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

COMMON=(
  --bin "$BIN"
  --instances "$INSTANCES"
  --parallel "$PARALLEL"
  --timelimit "$TIMELIMIT"
  --threads 1
  --overtime-grace "$OVERTIME_GRACE"
  --opts "-maxNodes 1000000"
)

declare -A CONDITIONS=(
  [C0_baseline]=""
  [C1_factor_0.5]="-perRoundNzCutLimitFactor 0.5"
  [C2_factor_1.0]="-perRoundNzCutLimitFactor 1.0"
  [C3_factor_2.0]="-perRoundNzCutLimitFactor 2.0"
  [C4_factor_4.0]="-perRoundNzCutLimitFactor 4.0"
)

ORDER=(C0_baseline C1_factor_0.5 C2_factor_1.0 C3_factor_2.0 C4_factor_4.0)

total=${#ORDER[@]}
idx=0
for cond in "${ORDER[@]}"; do
  idx=$((idx + 1))
  opts="${CONDITIONS[$cond]}"
  outdir="$OUTBASE/$cond"
  echo ""
  echo "══════════════════════════════════════════════════════════"
  echo "  Condition $idx/$total: $cond"
  [[ -n "$opts" ]] && echo "  Extra opts: $opts"
  echo "  Output: $outdir"
  echo "══════════════════════════════════════════════════════════"
  if [[ -n "$opts" ]]; then
    "$RUN_EXP" "${COMMON[@]}" --opts "$opts" --outdir "$outdir"
  else
    "$RUN_EXP" "${COMMON[@]}" --outdir "$outdir"
  fi
done

echo ""
echo "All $total conditions complete. Results in: $OUTBASE"
