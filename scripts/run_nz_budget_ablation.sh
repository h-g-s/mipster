#!/usr/bin/env bash
# run_nz_budget_ablation.sh
#
# Ablation study for the perRoundNzCutLimitFactor parameter, which limits
# the total non-zeros (NZ) added to the LP per root cut-generation round.
#
# The budget formula is:  budget = factor × max(currentNZ/10, 2×ncols + 100)
#
# When the budget is exceeded after running cut generators, row cuts are scored
# by the active ranking metric (fitness or violation) and greedily selected from
# best to worst until the budget is consumed.  Column cuts are always preserved.
# When stored (global) cuts already fill the budget, cut generators are skipped.
#
# A factor of -1.0 (default) disables the budget entirely (no trimming).
#
# Conditions
# ──────────────────────────────────────────────────────────────────────────────
#  C0_baseline       : -perRoundNzCutLimitFactor -1  (disabled, default)
#  C1_factor_2.0     : -perRoundNzCutLimitFactor 2.0 (permissive: 2× formula)
#  C2_factor_1.0     : -perRoundNzCutLimitFactor 1.0 (same formula as old global-cuts
#                       maximumAdd, but NOW also trims generator cuts — more restrictive
#                       than old behaviour, which only limited stored global cuts)
#  C3_factor_0.5     : -perRoundNzCutLimitFactor 0.5 (moderate trimming)
#  C4_factor_0.2     : -perRoundNzCutLimitFactor 0.2 (aggressive trimming)
# ──────────────────────────────────────────────────────────────────────────────
#
# Key metrics to compare (all found in stats.csv and summary.tsv):
#   cut_time         — total root cut-generation wall time
#   lp_seconds       — LP reoptimisation time during cut generation
#   continuous       — LP bound achieved after cut generation
#   time             — total solve time
#   result           — SOLVED / TIMEOUT / etc.
#
# Usage:
#   ./scripts/run_nz_budget_ablation.sh [--bin PATH] [--timelimit T]
#                                        [--parallel N] [--dry-run]
#
# Defaults:
#   --bin        $MIPSTER_PREFIX/bin/mipster
#   --timelimit  600
#   --parallel   $(nproc)
#   --instances  $MIPSTER_INSTANCES/miplib/2017+spp
#
# After the run, compare with:
#   python3 scripts/compare_multi_experiments.py --dir $OUTDIR

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
RUN_EXP="$SCRIPT_DIR/run_experiments.sh"

# ── Defaults ──────────────────────────────────────────────────────────────────
BIN="${MIPSTER_PREFIX:-$HOME/prog/mipster}/bin/mipster"
TIMELIMIT=600
PARALLEL=$(nproc)
OVERTIME_GRACE=600
INSTANCES="${MIPSTER_INSTANCES:-$HOME/inst}/miplib/2017+spp"
DATE=$(date +%Y_%m_%d)
OUTDIR="${MIPSTER_EXPERIMENTS:-$HOME/experiments/cbc}/nz_budget_${DATE}"
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
echo "  perRoundNzCutLimitFactor ablation study"
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
  "C0_baseline|-perRoundNzCutLimitFactor -1"
  "C1_factor_2.0|-perRoundNzCutLimitFactor 2.0"
  "C2_factor_1.0|-perRoundNzCutLimitFactor 1.0"
  "C3_factor_0.5|-perRoundNzCutLimitFactor 0.5"
  "C4_factor_0.2|-perRoundNzCutLimitFactor 0.2"
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
  echo "  Extra opts: $COPTS"
  echo "────────────────────────────────────────────────────────────────"

  # Build --opts arguments (one per flag-value pair).
  # A token is treated as a value (not a new flag) when it starts with a digit
  # or with '-' followed by a digit (negative numbers like -1, -0.5).
  OPTS_ARGS=()
  read -ra TOKENS <<< "$COPTS"
  i=0
  while [[ $i -lt ${#TOKENS[@]} ]]; do
    tok="${TOKENS[$i]}"
    next_is_value=0
    if [[ $((i+1)) -lt ${#TOKENS[@]} ]]; then
      nxt="${TOKENS[$((i+1))]}"
      # Value: doesn't start with '-', OR starts with '-' followed by a digit
      if [[ "$nxt" != -* ]] || [[ "$nxt" =~ ^-[0-9] ]]; then
        next_is_value=1
      fi
    fi
    if [[ $next_is_value -eq 1 ]]; then
      OPTS_ARGS+=(--opts "${tok} ${TOKENS[$((i+1))]}")
      i=$((i + 2))
    else
      OPTS_ARGS+=(--opts "$tok")
      i=$((i + 1))
    fi
  done

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
echo "  To compare conditions:"
echo "    python3 scripts/compare_multi_experiments.py --dir $OUTDIR"
echo "════════════════════════════════════════════════════════════════"
