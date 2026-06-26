#!/bin/bash
# compare_probing_filter.sh — compare mipster with/without probing row cut filter
# Runs instances in parallel (12 at a time), 1-hour time limit each.
# Usage: ./scripts/compare_probing_filter.sh [OUTDIR]

set -euo pipefail

OUTDIR="${1:-$HOME/experiments/cbc/probing_filter_$(date +%Y%m%d_%H%M%S)}"
INST_DIR="${MIPSTER_INSTANCES}/miplib/2017+spp"
TIMELIMIT=3600
JOBS=12

BIN_FILTER="$(dirname "$0")/../src/mipster_probfilter"
BIN_NOFILTER="$(dirname "$0")/../src/mipster_nofilter"

# Instances with slow preprocessing (server > 1000s, proxy for locally slow)
INSTANCES=(
  ex9 supportcase10 cdma sing326 rail01 usafa physiciansched3-3 bab6
  chromaticindex1024-7 s100 rocII-5-11 netdiversion ivu06-big splice1k1
  ns1116954 rocII-4-11 rocII-3-11 rocII-2-11 ns1208400 neos-4722843-widden
)

mkdir -p "$OUTDIR/filter" "$OUTDIR/nofilter"
echo "Output: $OUTDIR"
echo "Jobs: $JOBS  Timelimit: ${TIMELIMIT}s"

run_instance() {
  local bin="$1" tag="$2" inst="$3"
  local mps="$INST_DIR/${inst}.mps.gz"
  local log="$OUTDIR/$tag/${inst}.log"
  local res="$OUTDIR/$tag/${inst}.result"

  [ -f "$mps" ] || { echo "MISSING: $mps"; return; }

  timeout $((TIMELIMIT + 60)) "$bin" "$mps" -seconds "$TIMELIMIT" -solve \
    > "$log" 2>&1 || true

  # Extract key metrics
  local preptime obj bound gap nodes
  preptime=$(grep -oP 'Preprocessing complete.*Time:\s*\K[\d.]+' "$log" || echo "-")
  obj=$(grep -oP 'Objective value:\s*\K[-\d.eE+]+' "$log" | tail -1 || echo "-")
  bound=$(grep -oP 'Lower bound:\s*\K[-\d.eE+]+' "$log" | tail -1 || echo "-")
  gap=$(grep -oP 'Gap:\s*\K[\d.]+' "$log" | tail -1 || echo "-")
  nodes=$(grep -oP 'Enumerated nodes:\s*\K\d+' "$log" | tail -1 || echo "-")
  status=$(grep -oP 'Result - \K[^\n]+' "$log" | tail -1 || echo "?")

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$inst" "$preptime" "$obj" "$bound" "$gap" "$nodes" "$status" > "$res"
  echo "  [$tag] $inst: prep=${preptime}s  gap=${gap}%  nodes=$nodes  $status"
}

export -f run_instance
export OUTDIR TIMELIMIT BIN_FILTER BIN_NOFILTER INST_DIR

echo "=== Launching filter runs ==="
printf '%s\n' "${INSTANCES[@]}" | \
  parallel -j"$JOBS" run_instance "$BIN_FILTER" filter {}

echo "=== Launching no-filter runs ==="
printf '%s\n' "${INSTANCES[@]}" | \
  parallel -j"$JOBS" run_instance "$BIN_NOFILTER" nofilter {}

echo ""
echo "=== Results comparison ==="
printf "%-30s  %8s %8s  %8s %8s  %8s %8s  %8s\n" \
  "Instance" "prep_F" "prep_NF" "gap_F%" "gap_NF%" "nodes_F" "nodes_NF" "status_F"
echo "$(printf '%.0s-' {1..100})"

for inst in "${INSTANCES[@]}"; do
  rf="$OUTDIR/filter/${inst}.result"
  rn="$OUTDIR/nofilter/${inst}.result"
  [ -f "$rf" ] || continue
  read -r _ pf of bf gf nf sf < "$rf" 2>/dev/null || continue
  read -r _ pn on bn gn nn sn < "$rn" 2>/dev/null || pn="-"; gn="-"; nn="-"; sn="-"
  printf "%-30s  %8s %8s  %8s %8s  %8s %8s  %s\n" \
    "$inst" "$pf" "$pn" "$gf" "$gn" "$nf" "$nn" "$sf"
done

echo ""
echo "Done. Full logs in $OUTDIR"
