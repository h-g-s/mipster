#!/usr/bin/env bash
# run_parity_pathaggr_ablation.sh
#
# Ablation study for two recently added features:
#
#   1. parityPresolve — GF(2) parity reduction in OsiPresolve (CglPreProcess).
#      Controlled by:  -parityPresolve off
#      Note: ClpPresolve also runs GF(2) parity unconditionally (no param yet).
#
#   2. pathAggrCuts — Path-aggregation MIR cuts (CglPathAggregation).
#      Default mode:   root (root node only)
#      Options:        off | root | ifmove | on
#
# Conditions
# ──────────────────────────────────────────────────────────────────────────────
#  C0_baseline          : defaults (parityPresolve on, pathAggrCuts root)
#  C1_no_parity         : -parityPresolve off
#  C2_pathaggr_off      : -pathAggrCuts off
#  C3_pathaggr_ifmove   : -pathAggrCuts ifmove   (tree-wide when improving)
#  C4_pathaggr_on       : -pathAggrCuts on        (always generate)
#  C5_no_parity_no_pathaggr : -parityPresolve off -pathAggrCuts off
# ──────────────────────────────────────────────────────────────────────────────
#
# Usage:
#   ./scripts/run_parity_pathaggr_ablation.sh [--bin PATH] [--timelimit T]
#                                              [--parallel N] [--dry-run]
#
# Defaults:
#   --bin        $MIPSTER_PREFIX/bin/mipster
#   --timelimit  300
#   --parallel   $(nproc)
#   --instances  $MIPSTER_INSTANCES/miplib/2017+spp
#
# After the run, compare conditions with (once compare_experiments.py is available):
#   python3 scripts/compare_multi_experiments.py --dir $OUTDIR

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
RUN_EXP="$SCRIPT_DIR/run_experiments.sh"

# ── Defaults ──────────────────────────────────────────────────────────────────
BIN="${MIPSTER_PREFIX:-$HOME/prog/cbc}/bin/mipster"
TIMELIMIT=300
PARALLEL=$(nproc)
OVERTIME_GRACE=600
INSTANCES="${MIPSTER_INSTANCES:-$HOME/inst}/miplib/2017+spp"
DATE=$(date +%Y_%m_%d)
OUTDIR="${MIPSTER_EXPERIMENTS:-$HOME/experiments/cbc}/parity_pathaggr_${DATE}"
DRY_RUN=0

# ── Argument parsing ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bin)            BIN="$2";            shift 2 ;;
    --timelimit)      TIMELIMIT="$2";      shift 2 ;;
    --parallel)       PARALLEL="$2";       shift 2 ;;
    --overtime-grace) OVERTIME_GRACE="$2"; shift 2 ;;
    --instances)      INSTANCES="$2";      shift 2 ;;
    --outdir)         OUTDIR="$2";         shift 2 ;;
    --dry-run)        DRY_RUN=1;           shift   ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# ── Validate ──────────────────────────────────────────────────────────────────
if [[ ! -x "$BIN" ]]; then
  echo "Error: mipster binary not found or not executable: $BIN" >&2
  echo "Set \$MIPSTER_PREFIX or pass --bin PATH" >&2
  exit 1
fi

if [[ ! -d "$INSTANCES" ]]; then
  echo "Error: instances directory not found: $INSTANCES" >&2
  exit 1
fi

mkdir -p "$OUTDIR"

echo "════════════════════════════════════════════════════════════════"
echo "  Parity + PathAggr ablation study"
echo "  Binary:         $BIN"
echo "  Instances:      $INSTANCES"
echo "  Timelimit:      ${TIMELIMIT}s  (+${OVERTIME_GRACE}s overtime grace)"
echo "  Parallel:       $PARALLEL concurrent instances per condition"
echo "  Output:         $OUTDIR"
echo "  Dry-run:        $DRY_RUN"
echo "════════════════════════════════════════════════════════════════"
echo ""

# ── Condition table ────────────────────────────────────────────────────────────
# Format: "CONDITION_NAME|extra opts..."
CONDITIONS=(
  "C0_baseline|"
  "C1_no_parity|-parityPresolve off"
  "C2_pathaggr_off|-pathAggrCuts off"
  "C3_pathaggr_ifmove|-pathAggrCuts ifmove"
  "C4_pathaggr_on|-pathAggrCuts on"
  "C5_no_parity_no_pathaggr|-parityPresolve off -pathAggrCuts off"
)

# ── Runner ────────────────────────────────────────────────────────────────────
TOTAL=${#CONDITIONS[@]}
IDX=0

for entry in "${CONDITIONS[@]}"; do
  IDX=$((IDX + 1))
  CNAME="${entry%%|*}"
  COPTS="${entry#*|}"

  CDIR="$OUTDIR/$CNAME"

  echo "────────────────────────────────────────────────────────────────"
  echo "  [$IDX/$TOTAL] $CNAME"
  if [[ -n "$COPTS" ]]; then
    echo "  Extra opts: $COPTS"
  else
    echo "  Extra opts: (none — default settings)"
  fi
  echo "────────────────────────────────────────────────────────────────"

  # Build --opts arguments (one per flag token)
  OPTS_ARGS=()
  if [[ -n "$COPTS" ]]; then
    # Split COPTS into space-separated tokens and pair flags with values
    read -ra TOKENS <<< "$COPTS"
    i=0
    while [[ $i -lt ${#TOKENS[@]} ]]; do
      tok="${TOKENS[$i]}"
      # If next token exists and doesn't start with '-', it's the value
      if [[ $((i+1)) -lt ${#TOKENS[@]} && "${TOKENS[$((i+1))]}" != -* ]]; then
        OPTS_ARGS+=(--opts "${tok} ${TOKENS[$((i+1))]}")
        i=$((i + 2))
      else
        OPTS_ARGS+=(--opts "$tok")
        i=$((i + 1))
      fi
    done
  fi

  CMD=(
    bash "$RUN_EXP"
    --bin "$BIN"
    --timelimit "$TIMELIMIT"
    --overtime-grace "$OVERTIME_GRACE"
    --parallel "$PARALLEL"
    --threads 1
    --instances "$INSTANCES"
    --outdir "$CDIR"
    "${OPTS_ARGS[@]}"
  )

  echo "  Command: ${CMD[*]}"
  echo ""

  if [[ $DRY_RUN -eq 1 ]]; then
    echo "  [dry-run] skipping execution"
  else
    "${CMD[@]}"
  fi

  echo ""
done

echo "════════════════════════════════════════════════════════════════"
echo "  All conditions complete!"
echo "  Results in: $OUTDIR"
echo ""
echo "  To compare (once compare scripts are available):"
echo "    python3 scripts/compare_multi_experiments.py --dir $OUTDIR"
echo "════════════════════════════════════════════════════════════════"
