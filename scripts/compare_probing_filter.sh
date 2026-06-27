#!/bin/bash
# compare_probing_filter.sh — compare mipster with/without probing row cut filter
# Usage: ./scripts/compare_probing_filter.sh [OUTDIR]

OUTDIR="${1:-$HOME/experiments/cbc/probing_filter_$(date +%Y%m%d_%H%M%S)}"
INST_DIR="${MIPSTER_INSTANCES}/miplib/2017+spp"
TIMELIMIT=3600
OVERTIME=600   # hard-kill 10 min after time limit
REPO="$(cd "$(dirname "$0")/.." && pwd)"

BIN_FILTER="$REPO/src/mipster_probfilter"
BIN_NOFILTER="$REPO/src/mipster_nofilter"

INSTANCES=(
  ex9 supportcase10 cdma sing326 rail01 usafa physiciansched3-3 bab6
  chromaticindex1024-7 s100 rocII-5-11 netdiversion ivu06-big splice1k1
  ns1116954 ns1208400 neos-4722843-widden
)

mkdir -p "$OUTDIR/filter" "$OUTDIR/nofilter"
echo "Output:    $OUTDIR"
echo "Timelimit: ${TIMELIMIT}s   Overtime kill: +${OVERTIME}s"
echo "Instances: ${#INSTANCES[@]}"

run_one() {
    local bin="$1" tagdir="$2" inst="$3"
    local mps="$INST_DIR/${inst}.mps.gz"
    local log="$tagdir/${inst}.log"
    local res="$tagdir/${inst}.result"

    if [ ! -f "$mps" ]; then
        printf "MISSING\t-\t-\t-\t-\t-\tMISSING_MPS\n" > "$res"
        echo "  MISSING: $inst"; return
    fi

    timeout --kill-after=30 $((TIMELIMIT + OVERTIME)) \
        "$bin" "$mps" -seconds "$TIMELIMIT" -solve > "$log" 2>&1 || true

    local preptime obj bound gap nodes walltime status
    preptime=$(grep -oP 'Preprocess.*Time:\s*\K[\d.]+' "$log" | tail -1); [ -z "$preptime" ] && preptime="-"
    obj=$(grep -oP 'Objective value:\s*\K[-\d.eE+]+' "$log" | tail -1); [ -z "$obj" ] && obj="-"
    bound=$(grep -oP 'Lower bound:\s*\K[-\d.eE+]+' "$log" | tail -1); [ -z "$bound" ] && bound="-"
    gap=$(grep -oP 'Gap:\s*\K[\d.eE+]+' "$log" | tail -1); [ -z "$gap" ] && gap="-"
    nodes=$(grep -oP 'Enumerated nodes:\s*\K\d+' "$log" | tail -1); [ -z "$nodes" ] && nodes="-"
    walltime=$(grep -oP 'Total time \(Wallclock seconds\):\s*\K[\d.]+' "$log" | tail -1); [ -z "$walltime" ] && walltime="-"
    status=$(grep -oP 'Result - \K.+' "$log" | tail -1); [ -z "$status" ] && status="?"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$inst" "$preptime" "$obj" "$bound" "$gap" "$nodes" "$walltime" "$status" > "$res"
    echo "  DONE [$(basename "$tagdir")] $inst: prep=${preptime}s gap=${gap}% nodes=$nodes wall=${walltime}s"
}

# Launch runs in batches capped at JOBS to avoid OOM
PIDS=()
slot() {
    # Wait until fewer than JOBS background jobs are running
    while [ "${#PIDS[@]}" -ge "$JOBS" ]; do
        NEWPIDS=()
        for p in "${PIDS[@]}"; do kill -0 "$p" 2>/dev/null && NEWPIDS+=("$p"); done
        PIDS=("${NEWPIDS[@]}")
        [ "${#PIDS[@]}" -ge "$JOBS" ] && sleep 5
    done
}

echo "=== Launching filter runs ==="
for inst in "${INSTANCES[@]}"; do
    slot; ( run_one "$BIN_FILTER" "$OUTDIR/filter" "$inst" ) & PIDS+=($!)
done

echo "=== Launching nofilter runs ==="
for inst in "${INSTANCES[@]}"; do
    slot; ( run_one "$BIN_NOFILTER" "$OUTDIR/nofilter" "$inst" ) & PIDS+=($!)
done

echo "Waiting for all ${#PIDS[@]} active runs to finish..."
wait
echo ""
echo "=== Results comparison ==="
printf "%-30s  %7s %7s  %7s %7s  %8s %8s  %7s %7s\n" \
    "Instance" "prep_F" "prep_NF" "gap_F" "gap_NF" "nodes_F" "nodes_NF" "wall_F" "wall_NF"
printf '%0.s-' {1..110}; echo

for inst in "${INSTANCES[@]}"; do
    rf="$OUTDIR/filter/${inst}.result"
    rn="$OUTDIR/nofilter/${inst}.result"
    [ -f "$rf" ] || { echo "MISSING result: $inst filter"; continue; }
    IFS=$'\t' read -r _ pf of bf gf nf wf sf < "$rf" || true
    IFS=$'\t' read -r _ pn on bn gn nn wn sn < "$rn" 2>/dev/null || pn="-"; gn="-"; nn="-"; wn="-"
    printf "%-30s  %7s %7s  %7s %7s  %8s %8s  %7s %7s\n" \
        "$inst" "$pf" "$pn" "$gf" "$gn" "$nf" "$nn" "$wf" "$wn"
done

echo ""
echo "Full logs: $OUTDIR"
