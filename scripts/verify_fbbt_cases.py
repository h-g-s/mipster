#!/usr/bin/env python3
"""
verify_fbbt_cases.py  --  LP-verify bound-tightening cases from mipster_bp_bench.

Each line in the input JSON-lines file represents one bound-tightening event:

  { "beff":..., "isUB":0|1, "isBinaryFix":0|1,
    "claimedBound":..., "oldLB":..., "oldUB":...,
    "tightenedIdx":N, "rowIdx":N,
    "vars": [[coef, lb, ub, type], ...] }

Verification:
  - Build the LP: max/min x[tightenedIdx]
    s.t.  sum(coef[k] * x[k]) <= beff
          lb[k] <= x[k] <= ub[k]   for all k
  - Integer variables: apply floor (UB tightening) or ceil (LB tightening)
    to the LP optimum before comparison.
  - Tolerance: 1e-6 relative to max(1, |claimedBound|).

Usage:
  python3 verify_fbbt_cases.py cases.jsonl [--sample N] [--verbose] [--fail-fast]

Options:
  --sample N    Randomly sample at most N cases per category (default: all)
  --verbose     Print details of every case
  --fail-fast   Stop on first failure
"""

import json
import math
import random
import sys
from collections import defaultdict

try:
    from scipy.optimize import linprog
except ImportError:
    print("ERROR: scipy not installed. Run: pip install scipy", file=sys.stderr)
    sys.exit(2)

TOL = 1e-6
PRIMAL_TOL = 1e-7

# Variable types matching CoinColumnType enum in CoinPackedMatrix.hpp
TYPE_CONTINUOUS = 0
TYPE_BINARY = 1
TYPE_GENERAL_INT = 2
TYPE_SEMI_CONT = 3
TYPE_SEMI_INT = 4


def category(case: dict) -> str:
    """Return a readable category label for a case."""
    if case["isBinaryFix"]:
        return "a_binary_fix"
    lb, ub = case["oldLB"], case["oldUB"]
    var_type = case["vars"][case["tightenedIdx"]][3]
    kind = "int" if var_type == TYPE_GENERAL_INT else "cont"
    if lb >= -PRIMAL_TOL and ub > PRIMAL_TOL:
        return f"b_[0,b]_{kind}"
    if lb < -PRIMAL_TOL and ub <= PRIMAL_TOL:
        return f"c_[b,0]_{kind}"
    if lb < -PRIMAL_TOL and ub > PRIMAL_TOL:
        return f"d_cross_{kind}"
    return f"other_{kind}"


def verify_case(case: dict) -> tuple[bool, str]:
    """
    Returns (ok, message).  ok=True means claimed bound is verified.
    """
    vars_ = case["vars"]
    n = len(vars_)
    if n == 0:
        return False, "empty vars list"

    ti = case["tightenedIdx"]
    if ti < 0 or ti >= n:
        return False, f"tightenedIdx={ti} out of range [0,{n})"

    coefs = [v[0] for v in vars_]
    bounds = [(v[1], v[2]) for v in vars_]
    types = [v[3] for v in vars_]
    beff = case["beff"]
    is_ub = bool(case["isUB"])
    claimed = case["claimedBound"]

    # Objective: max x[ti]  (isUB)  or  min x[ti]  (!isUB)
    c_obj = [0.0] * n
    c_obj[ti] = -1.0 if is_ub else 1.0  # linprog always minimises

    res = linprog(c_obj,
                  A_ub=[coefs], b_ub=[beff],
                  bounds=bounds,
                  method="highs")

    if res.status == 2:
        return False, "LP infeasible — beff/bounds inconsistent"
    if res.status == 3:
        return False, "LP unbounded"
    if res.status != 0:
        return False, f"LP solver status={res.status}: {res.message}"

    lp_opt = -res.fun if is_ub else res.fun

    # Apply integrality rounding.
    if types[ti] in (TYPE_GENERAL_INT, TYPE_BINARY):
        lp_opt = math.floor(lp_opt + PRIMAL_TOL) if is_ub else math.ceil(lp_opt - PRIMAL_TOL)

    tol = TOL * max(1.0, abs(claimed))
    ok = abs(lp_opt - claimed) <= tol
    if ok:
        msg = f"OK  lp={lp_opt:.10g}  claimed={claimed:.10g}"
    else:
        msg = (f"FAIL  lp={lp_opt:.10g}  claimed={claimed:.10g}  "
               f"diff={abs(lp_opt - claimed):.3e}  tol={tol:.3e}  "
               f"rowIdx={case['rowIdx']}  ti={ti}  isUB={is_ub}")
    return ok, msg


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cases_file", help="JSON-lines file from mipster_bp_bench --collect-cases")
    ap.add_argument("--sample", type=int, default=0,
                    help="max cases per category to verify (0 = all)")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--fail-fast", action="store_true")
    args = ap.parse_args()

    # Load all cases.
    all_cases: list[dict] = []
    with open(args.cases_file) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                all_cases.append(json.loads(line))
            except json.JSONDecodeError as e:
                print(f"Warning: bad JSON on line {lineno}: {e}", file=sys.stderr)

    print(f"Loaded {len(all_cases)} cases from {args.cases_file}")

    # Group by category.
    by_cat: dict[str, list] = defaultdict(list)
    for c in all_cases:
        by_cat[category(c)].append(c)

    print(f"\nCategory distribution:")
    for cat, cases in sorted(by_cat.items()):
        print(f"  {cat:30s}: {len(cases):6d} cases")

    # Sample per category if requested.
    to_verify: list[dict] = []
    for cat, cases in sorted(by_cat.items()):
        if args.sample > 0 and len(cases) > args.sample:
            sampled = random.sample(cases, args.sample)
            print(f"  Sampling {args.sample} from {cat}")
        else:
            sampled = cases
        to_verify.extend(sampled)

    print(f"\nVerifying {len(to_verify)} cases...\n")

    passed = 0
    failed = 0
    errors: list[str] = []

    for i, case in enumerate(to_verify):
        ok, msg = verify_case(case)
        if ok:
            passed += 1
            if args.verbose:
                print(f"[{i+1:5d}] {category(case):30s}  {msg}")
        else:
            failed += 1
            errors.append(f"[{i+1:5d}] {category(case):30s}  {msg}")
            print(f"[{i+1:5d}] {category(case):30s}  {msg}")
            if args.fail_fast:
                break

    print(f"\n{'='*60}")
    print(f"Results: {passed} PASSED, {failed} FAILED out of {passed+failed} verified")

    if failed == 0:
        print("✓ All bound-tightening cases verified correct.")
        return 0
    else:
        print(f"\nFirst failures:")
        for e in errors[:10]:
            print(f"  {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
