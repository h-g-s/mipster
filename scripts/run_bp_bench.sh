#!/usr/bin/env bash
# run_bp_bench.sh — run bound propagation benchmark over an instance set
#
# Usage:
#   ./scripts/run_bp_bench.sh [OPTIONS] [INSTANCE_DIR]
#
# Arguments:
#   INSTANCE_DIR   Directory containing .mps.gz instances
#                  (default: $MIPSTER_INSTANCES/miplib/2017+spp)
#
# Options:
#   --level <singletons|milpbt|fixpoint>   BP level (default: fixpoint)
#   --max-rounds <N>                        Max rounds for milpbt (default: 100)
#   --out <file.csv>                        Output CSV file (default: bp_bench.csv in cwd)
#   --jobs <N>                              Parallel jobs (default: nproc)
#   --bin <path>                            Path to mipster_bp_bench binary
#                                           (default: looks in PATH then build tree)
#
# Example:
#   ./scripts/run_bp_bench.sh --level fixpoint --out /tmp/bp_results.csv
#   ./scripts/run_bp_bench.sh --level milpbt --max-rounds 50
#
# To compare levels, run three times with different --level flags and
# combine the CSV files (the 'level' column distinguishes them).

set -euo pipefail

# ── Defaults ────────────────────────────────────────────────────────────────
LEVEL="fixpoint"
MAX_ROUNDS=100
OUT="bp_bench.csv"
JOBS=$(nproc)
INSTANCE_DIR="${MIPSTER_INSTANCES:-$HOME/inst}/miplib/2017+spp"
BP_BIN=""

# ── Parse arguments ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --level)       LEVEL="$2";       shift 2 ;;
    --max-rounds)  MAX_ROUNDS="$2";  shift 2 ;;
    --out)         OUT="$2";         shift 2 ;;
    --jobs)        JOBS="$2";        shift 2 ;;
    --bin)         BP_BIN="$2";      shift 2 ;;
    -*)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
    *)
      INSTANCE_DIR="$1"
      shift
      ;;
  esac
done

# ── Locate binary ────────────────────────────────────────────────────────────
if [[ -z "$BP_BIN" ]]; then
  if command -v mipster_bp_bench &>/dev/null; then
    BP_BIN="mipster_bp_bench"
  elif [[ -x "$(dirname "$0")/../src/mipster_bp_bench" ]]; then
    BP_BIN="$(dirname "$0")/../src/mipster_bp_bench"
  else
    echo "Error: mipster_bp_bench not found in PATH or build tree." >&2
    echo "Build with: cd src && make -j\$(nproc) mipster_bp_bench" >&2
    exit 1
  fi
fi

# ── Validate inputs ──────────────────────────────────────────────────────────
if [[ ! -d "$INSTANCE_DIR" ]]; then
  echo "Error: instance directory not found: $INSTANCE_DIR" >&2
  exit 1
fi

INSTANCES=( "$INSTANCE_DIR"/*.mps.gz )
if [[ ${#INSTANCES[@]} -eq 0 || ! -f "${INSTANCES[0]}" ]]; then
  echo "Error: no .mps.gz files found in $INSTANCE_DIR" >&2
  exit 1
fi

echo "Instance dir : $INSTANCE_DIR"
echo "Instances    : ${#INSTANCES[@]}"
echo "Level        : $LEVEL"
[[ "$LEVEL" == "milpbt" ]] && echo "Max rounds   : $MAX_ROUNDS"
echo "Jobs         : $JOBS"
echo "Output       : $OUT"
echo ""

# ── Write header ─────────────────────────────────────────────────────────────
"$BP_BIN" --header-only > "$OUT"

# ── Run in parallel ──────────────────────────────────────────────────────────
EXTRA_ARGS=()
[[ "$LEVEL" == "milpbt" ]] && EXTRA_ARGS+=("--max-rounds" "$MAX_ROUNDS")

parallel -j"$JOBS" \
  "$BP_BIN" --no-header --level "$LEVEL" "${EXTRA_ARGS[@]}" {} \
  >> "$OUT" \
  ::: "${INSTANCES[@]}"

echo "Done. Results written to: $OUT"
echo "Rows: $(( $(wc -l < "$OUT") - 1 )) (excluding header)"
