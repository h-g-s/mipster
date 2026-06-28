#!/bin/bash
# compare_probing_filter.sh
OUTDIR="${1:-$HOME/experiments/cbc/probing_filter_$(date +%Y%m%d_%H%M%S)}"
INST_DIR="${MIPSTER_INSTANCES}/miplib/2017+spp"
TIMELIMIT=3600
JOBS=7
REPO="$(cd "$(dirname "$0")/.." && pwd)"
RUNNER="$REPO/scripts/run_probing_one.sh"
BIN_FILTER="$REPO/src/mipster_probfilter"
BIN_NOFILTER="$REPO/src/mipster_nofilter"

INSTANCES=(
  ex9 supportcase10 cdma sing326 rail01 usafa physiciansched3-3 bab6
  chromaticindex1024-7 s100 rocII-5-11 netdiversion ivu06-big splice1k1
  ns1116954 ns1208400 neos-4722843-widden
)

# If not inside a tmux session, re-launch self inside one so the experiment
# survives terminal closure (parallel reinstalls SIGHUP handler, defeating nohup).
if [ -z "$TMUX" ]; then
    SESSION="probing_$(date +%H%M%S)"
    mkdir -p "$(dirname "$OUTDIR")"
    LOG="$HOME/experiments/cbc/probing_filter_launch.log"
    tmux new-session -d -s "$SESSION" \
        "bash '$0' '$OUTDIR' 2>&1 | tee '$LOG'"
    echo "Launched in tmux session '$SESSION'"
    echo "Monitor:  tmux attach -t $SESSION"
    echo "Log:      $LOG"
    echo "Output:   $OUTDIR"
    exit 0
fi

mkdir -p "$OUTDIR/filter" "$OUTDIR/nofilter"
echo "Output:    $OUTDIR"
echo "Jobs:      $JOBS   Timelimit: ${TIMELIMIT}s"
echo "Total runs: $((${#INSTANCES[@]}*2))"

# Build jobs list: one command per line, piped to parallel
JOBS_FILE=$(mktemp)
for inst in "${INSTANCES[@]}"; do
    echo "$RUNNER $BIN_FILTER   $OUTDIR/filter   $INST_DIR $inst $TIMELIMIT"
    echo "$RUNNER $BIN_NOFILTER $OUTDIR/nofilter $INST_DIR $inst $TIMELIMIT"
done > "$JOBS_FILE"

echo "Jobs file: $JOBS_FILE ($(wc -l < "$JOBS_FILE") lines)"
echo "=== Launching ==="
parallel --jobs "$JOBS" --no-run-if-empty < "$JOBS_FILE"
rm -f "$JOBS_FILE"

echo ""
echo "=== Results ==="
printf "%-30s  %7s %7s  %8s %8s  %7s %7s\n" \
    "Instance" "prep_F" "prep_NF" "nodes_F" "nodes_NF" "wall_F" "wall_NF"
printf '%0.s-' {1..90}; echo
for inst in "${INSTANCES[@]}"; do
    rf="$OUTDIR/filter/${inst}.result"
    rn="$OUTDIR/nofilter/${inst}.result"
    [ -f "$rf" ] || { printf "%-30s  MISSING\n" "$inst"; continue; }
    IFS=$'\t' read -r _ pf of bf gf nf wf sf < "$rf" || true
    [ -f "$rn" ] && { IFS=$'\t' read -r _ pn on bn gn nn wn sn < "$rn" || true; } || { pn="-"; nn="-"; wn="-"; }
    printf "%-30s  %7s %7s  %8s %8s  %7s %7s\n" \
        "$inst" "$pf" "$pn" "$nf" "$nn" "$wf" "$wn"
done
echo "Logs: $OUTDIR"
