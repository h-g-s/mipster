#!/usr/bin/env bash
#
# run_fbbt_cmp.sh — overnight FBBT correctness + performance comparison
#
# Two conditions per instance:
#   C0_baseline  — mipster without FBBT
#   C1_fbbt      — mipster with MIPSTER_FBBT=1 (enableFBBT at root)
#
# Usage: ./run_fbbt_cmp.sh [--instances FILE] [--out DIR] [--sec N] [--jobs N]
#
set -euo pipefail

INST_DIR="${MIPSTER_INSTANCES:-$HOME/inst/miplib/2017+spp}"
OUT_DIR=""
TIME_LIMIT=10800
JOBS=$(nproc)
INSTANCES_FILE=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --instances) INSTANCES_FILE="$2"; shift 2;;
    --out)       OUT_DIR="$2"; shift 2;;
    --sec)       TIME_LIMIT="$2"; shift 2;;
    --jobs)      JOBS="$2"; shift 2;;
    *) echo "Unknown option: $1"; exit 1;;
  esac
done

[[ -z "$OUT_DIR" ]] && OUT_DIR="${MIPSTER_EXPERIMENTS:-$HOME/experiments/cbc}/fbbt_cmp_$(date +%Y_%m_%d)"
mkdir -p "$OUT_DIR/C0_baseline" "$OUT_DIR/C1_fbbt"

MIPSTER_BIN="${MIPSTER_BIN:-mipster}"

# ── Instance list ──────────────────────────────────────────────────────────

if [[ -n "$INSTANCES_FILE" ]]; then
  mapfile -t INSTANCES < "$INSTANCES_FILE"
else
  mapfile -t INSTANCES < <(ls "$INST_DIR"/*.mps.gz 2>/dev/null)
fi

echo "Instances  : ${#INSTANCES[@]}"
echo "Time limit : ${TIME_LIMIT}s"
echo "Jobs       : ${JOBS}"
echo "Output     : $OUT_DIR"
echo "Started    : $(date)"
echo ""

# ── Single-instance runner ─────────────────────────────────────────────────
# Called as:  run_one <inst_path>  <C0|C1>  <out_dir>  <time_limit>  <mipster>
#
# Written as a shell script wrapper so GNU parallel can invoke it cleanly,
# avoiding env-var and argument-splitting issues.

run_one() {
  local inst_path="$1"
  local cond="$2"       # C0 or C1
  local out_dir="$3"
  local time_limit="$4"
  local mipster="$5"

  local name; name=$(basename "$inst_path" .mps.gz); name=$(basename "$name" .lp.gz)
  local log="$out_dir/$name.log"
  local sol="$out_dir/$name.sol"
  local result="$out_dir/$name.result"

  local t0=$SECONDS
  if [[ "$cond" == "C1" ]]; then
    env MIPSTER_FBBT=1 "$mipster" "$inst_path" -sec "$time_limit" \
        -solu "$sol" -solve > "$log" 2>&1 || true
  else
    "$mipster" "$inst_path" -sec "$time_limit" \
        -solu "$sol" -solve > "$log" 2>&1 || true
  fi
  local elapsed=$(( SECONDS - t0 ))

  # Parse log
  local status="UNKNOWN"
  local obj="NA"
  local bound="NA"
  local nodes="NA"

  if grep -q "Result - Optimal solution found" "$log" 2>/dev/null; then
    status="OPTIMAL"
  elif grep -q "Result - Stopped on time" "$log" 2>/dev/null; then
    status="TIMELIMIT"
  elif grep -q "Result - Infeasible" "$log" 2>/dev/null; then
    status="INFEASIBLE"
  elif grep -q "infeasibility proved" "$log" 2>/dev/null; then
    status="INFEASIBLE"
  fi

  obj=$(grep "Objective value:" "$log" 2>/dev/null | awk '{print $NF}' | head -1 || true)
  [[ -z "$obj" ]] && obj="NA"

  local last_tick
  last_tick=$(grep "^✔" "$log" 2>/dev/null | tail -1 || true)
  if [[ -n "$last_tick" ]]; then
    bound=$(echo "$last_tick" | grep -oP 'Bound:\s*\K[0-9.eE+\-]+' || true)
    nodes=$(echo "$last_tick" | grep -oP 'Nodes:\s*\K[0-9.KM]+' || true)
  fi
  [[ -z "$bound" ]] && bound="NA"
  [[ -z "$nodes" ]] && nodes="NA"

  printf "%s\t%s\t%s\t%d\t%s\n" "$status" "$obj" "$bound" "$elapsed" "$nodes" > "$result"

  printf "[%s] %-4s %-45s  %-12s  obj=%-18s  t=%ds  nodes=%s\n" \
    "$(date +%H:%M:%S)" "$cond" "$name" "$status" "$obj" "$elapsed" "$nodes"
}
export -f run_one
export MIPSTER_BIN

# ── Generate job wrapper scripts ───────────────────────────────────────────
# Write one small shell script per job so parallel just executes them.
# This avoids all arg-splitting issues.

JOB_DIR=$(mktemp -d)
for inst in "${INSTANCES[@]}"; do
  for cond in C0 C1; do
    out_subdir="$OUT_DIR/C0_baseline"
    [[ "$cond" == "C1" ]] && out_subdir="$OUT_DIR/C1_fbbt"
    job="$JOB_DIR/$(basename "$inst" .mps.gz)_${cond}.sh"
    printf '#!/usr/bin/env bash\nrun_one %q %s %q %d %q\n' \
      "$inst" "$cond" "$out_subdir" "$TIME_LIMIT" "$MIPSTER_BIN" > "$job"
    chmod +x "$job"
  done
done

echo "Total jobs : $(ls "$JOB_DIR"/*.sh | wc -l)"

parallel --no-notice -j "$JOBS" bash {} ::: "$JOB_DIR"/*.sh \
  2>"$OUT_DIR/parallel.log" || true

rm -rf "$JOB_DIR"

echo ""
echo "Finished: $(date)"
echo ""

# ── Comparison summary ─────────────────────────────────────────────────────

python3 - "$OUT_DIR" << 'PYEOF'
import sys, os, glob

out_dir = sys.argv[1]

def read_dir(path):
    d = {}
    for f in glob.glob(f"{path}/*.result"):
        name = os.path.basename(f).replace('.result','')
        with open(f) as fh:
            parts = fh.read().strip().split('\t')
        if len(parts) >= 4:
            d[name] = {'status': parts[0], 'obj': parts[1],
                       'bound': parts[2], 'time': parts[3],
                       'nodes': parts[4] if len(parts) > 4 else 'NA'}
    return d

c0 = read_dir(f"{out_dir}/C0_baseline")
c1 = read_dir(f"{out_dir}/C1_fbbt")
common = sorted(c0.keys() & c1.keys())

w = 45
hdr = (f"{'Instance':<{w}} {'C0_status':>10} {'C0_obj':>16} {'C0_t(s)':>7} {'C0_nodes':>8}"
       f"  {'C1_status':>10} {'C1_obj':>16} {'C1_t(s)':>7} {'C1_nodes':>8}  {'match':>14}  {'speedup':>7}")
print(hdr)
print('-' * len(hdr))

wrong = []
for name in common:
    r0, r1 = c0[name], c1[name]
    check = "OK"
    speedup = "NA"
    both_optimal = (r0['status'] == 'OPTIMAL' and r1['status'] == 'OPTIMAL')
    try:
        o0, o1 = float(r0['obj']), float(r1['obj'])
        diff = abs(o0 - o1)
        rel = diff / max(1.0, abs(o0))
        if both_optimal and rel > 1e-5:
            # Both proved optimal but with different values — real bug
            check = f"OBJ_DIFF({diff:.3g})"
            wrong.append(name)
        elif not both_optimal and rel > 1e-5:
            # Different incumbents at time limit — performance note, not a bug
            check = "C1_better" if o1 < o0 else "C0_better"
    except:
        if r1['status'] == 'INFEASIBLE' and r0['status'] != 'INFEASIBLE':
            check = "INFEAS_WRONG"
            wrong.append(name)
        elif r0['status'] == 'INFEASIBLE' and r1['status'] != 'INFEASIBLE':
            check = "INFEAS_WRONG_C0"
            wrong.append(name)
        elif r0['obj'] == 'NA' and r1['obj'] == 'NA':
            check = "both_no_sol"
        else:
            check = "parse_err"
    try:
        t0v, t1v = int(r0['time']), int(r1['time'])
        if t1v > 0:
            speedup = f"{t0v/t1v:.2f}x"
    except:
        pass
    print(f"{name:<{w}} {r0['status']:>10} {r0['obj']:>16} {r0['time']:>7} {r0['nodes']:>8}"
          f"  {r1['status']:>10} {r1['obj']:>16} {r1['time']:>7} {r1['nodes']:>8}  {check:>14}  {speedup:>7}")

print('-' * len(hdr))
s0 = sum(1 for r in c0.values() if r['status']=='OPTIMAL')
s1 = sum(1 for r in c1.values() if r['status']=='OPTIMAL')
print(f"\nTotal: {len(common)}   C0_optimal={s0}   C1_optimal={s1}   wrong={len(wrong)}")
if wrong:
    print(f"WRONG: {', '.join(wrong)}")
PYEOF
