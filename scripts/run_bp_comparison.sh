#!/usr/bin/env bash
# run_bp_comparison.sh — compare binary-only vs full-FBBT bound propagation
# on the full 2017+spp instance set.
#
# Usage:
#   ./scripts/run_bp_comparison.sh [--jobs N] [--out DIR] [--instances DIR]
#
# Output:
#   <out>/results.csv   — one row per (instance × condition)
#   <out>/summary.txt   — printed after all runs finish

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

BENCH="$REPO_DIR/src/mipster_bp_bench"
JOBS=$(nproc)
OUT_DIR="$REPO_DIR/experiments/bp_comparison_$(date +%Y%m%d_%H%M%S)"
INST_DIR="${MIPSTER_INSTANCES:-$HOME/inst}/miplib/2017+spp"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --jobs) JOBS="$2"; shift 2 ;;
    --out)  OUT_DIR="$2"; shift 2 ;;
    --instances) INST_DIR="$2"; shift 2 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

if [[ ! -x "$BENCH" ]]; then
  echo "ERROR: bench binary not found at $BENCH — run 'cd src && make' first." >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
RESULTS="$OUT_DIR/results.csv"

echo "Instance set : $INST_DIR"
echo "Output       : $OUT_DIR"
echo "Jobs         : $JOBS"

INSTANCES=( "$INST_DIR"/*.mps.gz )
N=${#INSTANCES[@]}
echo "Instances    : $N"
echo ""

# Write CSV header once.
"$BENCH" --header-only > "$RESULTS"

# Run both conditions in parallel. Each job appends one row.
# Two jobs per instance (full_fbbt + binary_only) → 2*N total tasks.
run_one() {
  local f="$1"
  local mode="$2"
  if [[ "$mode" == "binary_only" ]]; then
    "$BENCH" --no-header --binary-only "$f" 2>/dev/null
  else
    "$BENCH" --no-header "$f" 2>/dev/null
  fi
}
export -f run_one
export BENCH

# Build task list: "path mode" pairs
TASKS=()
for f in "${INSTANCES[@]}"; do
  TASKS+=("$f full_fbbt")
  TASKS+=("$f binary_only")
done

printf '%s\n' "${TASKS[@]}" | \
  parallel -j"$JOBS" --colsep ' ' run_one {1} {2} >> "$RESULTS"

echo "Done. Results -> $RESULTS"
echo ""

# Print summary using Python.
python3 - "$RESULTS" "$OUT_DIR" << 'PYEOF'
import sys, csv, os
from collections import defaultdict

results_file, out_dir = sys.argv[1], sys.argv[2]

rows = []
with open(results_file) as f:
    reader = csv.DictReader(f)
    for row in reader:
        rows.append(row)

by_inst = defaultdict(dict)
for r in rows:
    by_inst[r['instance']][r['condition']] = r

# Instances present in both conditions.
both = {k: v for k, v in by_inst.items()
        if 'full_fbbt' in v and 'binary_only' in v}

n = len(both)
total_fbbt   = sum(int(v['full_fbbt']['fbbt_tightened']) for v in both.values())
total_bpfix  = sum(int(v['full_fbbt']['bp_fixed'])       for v in both.values())
total_singfix= sum(int(v['full_fbbt']['singleton_fixed']) for v in both.values())

time_fbbt    = [float(v['full_fbbt']['time_sec'])   for v in both.values()]
time_binary  = [float(v['binary_only']['time_sec']) for v in both.values()]

slowdowns = [(tf - tb, tb, tf, inst)
             for inst, v in both.items()
             for tb, tf in [(float(v['binary_only']['time_sec']),
                             float(v['full_fbbt']['time_sec']))]]
slowdowns.sort(reverse=True)

# Instances where FBBT actually does something.
fbbt_active = {k: v for k, v in both.items()
               if int(v['full_fbbt']['fbbt_tightened']) > 0}

print(f"{'='*60}")
print(f"Bound Propagation: binary-only vs full-FBBT")
print(f"{'='*60}")
print(f"Instances compared    : {n}")
print(f"Instances with FBBT>0 : {len(fbbt_active)}")
print()
print(f"Tightenings (full_fbbt condition):")
print(f"  singleton_fixed : {total_singfix:,}")
print(f"  bp_fixed (bin)  : {total_bpfix:,}")
print(f"  fbbt_tightened  : {total_fbbt:,}")
print()
print(f"Time (seconds, all instances):")
print(f"  binary_only  total={sum(time_binary):.4f}  mean={sum(time_binary)/n:.6f}  "
      f"max={max(time_binary):.6f}")
print(f"  full_fbbt    total={sum(time_fbbt):.4f}  mean={sum(time_fbbt)/n:.6f}  "
      f"max={max(time_fbbt):.6f}")
overhead = sum(time_fbbt) - sum(time_binary)
print(f"  overhead     total={overhead:.4f}  "
      f"relative={100*overhead/max(sum(time_binary),1e-9):.1f}%")
print()

# Subset: only instances where FBBT fires
if fbbt_active:
    tf_a = [float(v['full_fbbt']['time_sec'])   for v in fbbt_active.values()]
    tb_a = [float(v['binary_only']['time_sec'])  for v in fbbt_active.values()]
    oh_a = sum(tf_a) - sum(tb_a)
    print(f"Time (instances where fbbt_tightened>0, n={len(fbbt_active)}):")
    print(f"  binary_only  total={sum(tb_a):.4f}  mean={sum(tb_a)/len(tb_a):.6f}")
    print(f"  full_fbbt    total={sum(tf_a):.4f}  mean={sum(tf_a)/len(tf_a):.6f}")
    print(f"  overhead     total={oh_a:.4f}  "
          f"relative={100*oh_a/max(sum(tb_a),1e-9):.1f}%")
    print()

print(f"Top 20 slowest (by absolute overhead  full_fbbt - binary_only):")
print(f"  {'instance':<40} {'binary_only':>12} {'full_fbbt':>10} {'overhead':>10} {'fbbt_tight':>10}")
print(f"  {'-'*40} {'-'*12} {'-'*10} {'-'*10} {'-'*10}")
for overhead_i, tb, tf, inst in slowdowns[:20]:
    ft = int(both[inst]['full_fbbt']['fbbt_tightened'])
    print(f"  {inst:<40} {tb:>12.6f} {tf:>10.6f} {overhead_i:>+10.6f} {ft:>10}")

# Save slow instances list for profiling.
slow_file = os.path.join(out_dir, 'slowest_instances.txt')
with open(slow_file, 'w') as f:
    for _, _, _, inst in slowdowns[:50]:
        if inst in both:
            f.write(both[inst]['full_fbbt'].get('instance', inst) + '\n')
print()
print(f"Top-50 slowest instance names -> {slow_file}")
PYEOF
