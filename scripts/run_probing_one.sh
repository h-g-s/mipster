#!/bin/bash
# run_probing_one.sh <bin> <tagdir> <inst_dir> <inst> <timelimit>
# Runs one solver instance and writes a .result file.
BIN="$1"
TAGDIR="$2"
INST_DIR="$3"
INST="$4"
TIMELIMIT="${5:-3600}"
OVERTIME=600

MPS="$INST_DIR/${INST}.mps.gz"
LOG="$TAGDIR/${INST}.log"
RES="$TAGDIR/${INST}.result"

if [ ! -f "$MPS" ]; then
    printf "%s\t-\t-\t-\t-\t-\t-\tMISSING_MPS\n" "$INST" > "$RES"
    echo "MISSING [$( basename "$TAGDIR")] $INST"
    exit 0
fi

timeout --kill-after=30 $((TIMELIMIT + OVERTIME)) \
    "$BIN" "$MPS" -seconds "$TIMELIMIT" -solve > "$LOG" 2>&1 || true

preptime=$(grep -oP 'Preprocess.*Time:\s*\K[\d.]+' "$LOG" | tail -1); [ -z "$preptime" ] && preptime="-"
obj=$(grep -oP 'Objective value:\s*\K[-\d.eE+]+' "$LOG" | tail -1);   [ -z "$obj" ]     && obj="-"
bound=$(grep -oP 'Lower bound:\s*\K[-\d.eE+]+' "$LOG" | tail -1);     [ -z "$bound" ]   && bound="-"
gap=$(grep -oP 'Gap:\s*\K[\d.eE+]+' "$LOG" | tail -1);                [ -z "$gap" ]     && gap="-"
nodes=$(grep -oP 'Enumerated nodes:\s*\K\d+' "$LOG" | tail -1);       [ -z "$nodes" ]   && nodes="-"
walltime=$(grep -oP 'Total time \(Wallclock seconds\):\s*\K[\d.]+' "$LOG" | tail -1); [ -z "$walltime" ] && walltime="-"
status=$(grep -oP 'Result - \K.+' "$LOG" | tail -1);                  [ -z "$status" ]  && status="?"

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$INST" "$preptime" "$obj" "$bound" "$gap" "$nodes" "$walltime" "$status" > "$RES"
echo "DONE [$(basename "$TAGDIR")] $INST: prep=${preptime}s nodes=$nodes wall=${walltime}s"
