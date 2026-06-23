#!/usr/bin/env python3
"""
stats_analysis.py — Detailed time/efficiency breakdown from MIPster stats.csv files.

Reads the per-instance stats CSV files produced by run_experiments.sh
(via -writeStatistics / -csvStatistics) and reports:

  1. Time budget breakdown  — lp_seconds, cut_time, cgraph_time,
                              total heuristic time, strong-branching, other
  2. Cut generator ranking  — by total time, cuts produced, avg NZ
  3. Heuristic ranking      — by total time and solutions found
  4. Worst instances        — top-N slowest overall and by phase

Usage:
    python3 stats_analysis.py --outdir /path/to/exp/
    python3 stats_analysis.py --outdir /path/to/exp/ --top 20 --plain
    python3 stats_analysis.py --statsfile /path/to/stats.csv --top 15
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from rich.box import ROUNDED
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

# ── Known cut generators (prefix in column names) ────────────────────────────

CUT_GENERATORS = [
    ("probing",             "Probing"),
    ("gomory",              "Gomory (root)"),
    ("gomoryl1",            "Gomory L1"),
    ("gomoryl2",            "Gomory L2"),
    ("gomory_2",            "Gomory (tree)"),
    ("knapsack",            "Knapsack"),
    ("clique",              "Clique (BK)"),
    ("oddwheel",            "OddWheel"),
    ("mixedintegerrounding2","MIR2"),
    ("flowcover",           "FlowCover"),
    ("twomircuts",          "TwoMIR"),
    ("twomircutsl1",        "TwoMIR L1"),
    ("twomircutsl2",        "TwoMIR L2"),
    ("liftandproject",      "LiftAndProject"),
    ("residualcapacity",    "ResidualCap"),
    ("zerohalf",            "ZeroHalf"),
    ("reduce_and_split",    "ReduceSplit"),
    ("reduce_and_split_2",  "ReduceSplit2"),
    ("stored",              "Stored cuts"),
]

# ── Known heuristics (prefix in column names) ─────────────────────────────────

HEURISTICS = [
    ("feasibility_pump",            "FeasibilityPump"),
    ("feasibilityjump",             "FeasibilityJump"),
    ("rins",                        "RINS"),
    ("rens",                        "RENS"),
    ("rensdj",                      "RENS-DJ"),
    ("rensub",                      "RENS-UB"),
    ("rounding",                    "Rounding"),
    ("random_rounding",             "RandomRounding"),
    ("combine_solutions",           "CombineSol"),
    ("greedy_cover",                "GreedyCover"),
    ("greedy_equality",             "GreedyEquality"),
    ("divecoefficient",             "DiveCoefficient"),
    ("divefractional",              "DiveFractional"),
    ("diveguided",                  "DiveGuided"),
    ("divelinesearch",              "DiveLineSearch"),
    ("divepseudocost",              "DivePseudoCost"),
    ("divevectorlength",            "DiveVectorLength"),
    ("diveany",                     "DiveAny"),
    ("vnd",                         "VND"),
    ("naive",                       "Naive"),
    ("multiple_root_solvers",       "MultipleRootSolvers"),
    ("dynamic_pass_thru",           "DynPassThru"),
    ("linked",                      "Linked"),
    ("partial_solution_given",      "PartialSolGiven"),
    ("dantzig_wolfe_expansion",     "DantzigWolfe"),
]


def _safe_sum(df: pd.DataFrame, col: str, default: float = 0.0) -> float:
    if col in df.columns:
        return pd.to_numeric(df[col], errors="coerce").fillna(0.0).sum()
    return default


def _safe_col(df: pd.DataFrame, col: str) -> pd.Series:
    if col in df.columns:
        return pd.to_numeric(df[col], errors="coerce").fillna(0.0)
    return pd.Series([0.0] * len(df), index=df.index)


def load_stats(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    # Rename first column to 'instance' for consistency
    if df.columns[0] != "instance":
        df = df.rename(columns={df.columns[0]: "instance"})
    return df


def fmt_time(seconds: float) -> str:
    if seconds < 0.001:
        return "0"
    if seconds >= 3600:
        return f"{seconds/3600:.2f}h"
    if seconds >= 60:
        return f"{seconds/60:.1f}min"
    return f"{seconds:.2f}s"


def fmt_pct(frac: float) -> str:
    return f"{frac * 100:.1f}%"


def bar(frac: float, width: int = 20) -> str:
    filled = round(frac * width)
    return "█" * filled + "░" * (width - filled)


# ─────────────────────────────────────────────────────────────────────────────
# Section 1: Time budget breakdown
# ─────────────────────────────────────────────────────────────────────────────

def section_time_budget(df: pd.DataFrame, console: Console) -> None:
    total_time  = _safe_sum(df, "time")
    lp_time     = _safe_sum(df, "lp_seconds")
    cut_time    = _safe_sum(df, "cut_time")
    cgraph_time = _safe_sum(df, "cgraph_time")

    # Sum all heuristic times
    heur_time = sum(
        _safe_sum(df, f"heur_{pfx}_time") for pfx, _ in HEURISTICS
    )

    # Strong branching / other LP time is captured in lp_seconds
    # Remaining = everything not accounted for
    accounted = lp_time + cut_time + cgraph_time + heur_time
    other_time = max(total_time - accounted, 0.0)

    if total_time <= 0:
        console.print("[dim]No timing data available.[/dim]")
        return

    phases = [
        ("LP solve",        lp_time,     "blue"),
        ("Cut generation",  cut_time,    "yellow"),
        ("Heuristics",      heur_time,   "green"),
        ("Conflict graph",  cgraph_time, "magenta"),
        ("Other",           other_time,  "dim"),
    ]

    table = Table(
        title="⏱  Time budget (summed across all instances)",
        box=ROUNDED, header_style="bold cyan", show_lines=False,
    )
    table.add_column("Phase",         style="bold",    width=22)
    table.add_column("Total time",    justify="right", width=12)
    table.add_column("Share",         justify="right", width=7)
    table.add_column("Distribution",  width=24)

    for name, t, color in phases:
        frac = t / total_time if total_time > 0 else 0.0
        table.add_row(
            name,
            fmt_time(t),
            fmt_pct(frac),
            f"[{color}]{bar(frac)}[/{color}]",
        )

    # Footer: total
    table.add_section()
    table.add_row("[bold]Total[/bold]", fmt_time(total_time), "100%", "")

    console.print(table)


# ─────────────────────────────────────────────────────────────────────────────
# Section 2: Cut generator ranking
# ─────────────────────────────────────────────────────────────────────────────

def section_cuts(df: pd.DataFrame, console: Console, top_n: int) -> None:
    rows = []
    for pfx, label in CUT_GENERATORS:
        time_col  = f"cut_{pfx}_time"
        cuts_col  = f"cut_{pfx}_cuts"
        calls_col = f"cut_{pfx}_calls"
        avgnz_col = f"cut_{pfx}_avgnz"

        t     = _safe_sum(df, time_col)
        cuts  = _safe_sum(df, cuts_col)
        calls = _safe_sum(df, calls_col)
        avgnz = (
            pd.to_numeric(df[avgnz_col], errors="coerce").mean()
            if avgnz_col in df.columns else float("nan")
        )

        if t > 0 or cuts > 0:
            rows.append({
                "label": label,
                "time":  t,
                "cuts":  int(cuts),
                "calls": int(calls),
                "avgnz": avgnz,
                "cuts_per_call": cuts / calls if calls > 0 else 0.0,
            })

    if not rows:
        console.print("[dim]No cut generator data found.[/dim]")
        return

    rows.sort(key=lambda r: r["time"], reverse=True)
    max_time = rows[0]["time"] if rows else 1.0

    table = Table(
        title=f"✂  Cut generators — top {min(top_n, len(rows))} by time",
        box=ROUNDED, header_style="bold cyan", show_lines=False,
    )
    table.add_column("Generator",    style="bold",    width=22)
    table.add_column("Total time",   justify="right", width=12)
    table.add_column("Bar",          width=20)
    table.add_column("Total cuts",   justify="right", width=11)
    table.add_column("Calls",        justify="right", width=8)
    table.add_column("Cuts/call",    justify="right", width=10)
    table.add_column("Avg NZ",       justify="right", width=9)

    for r in rows[:top_n]:
        frac = r["time"] / max_time if max_time > 0 else 0.0
        avgnz_str = f"{r['avgnz']:.1f}" if r["avgnz"] == r["avgnz"] else "-"
        table.add_row(
            r["label"],
            fmt_time(r["time"]),
            f"[yellow]{bar(frac)}[/yellow]",
            f"{r['cuts']:,}",
            f"{r['calls']:,}",
            f"{r['cuts_per_call']:.2f}",
            avgnz_str,
        )

    console.print(table)


# ─────────────────────────────────────────────────────────────────────────────
# Section 3: Heuristic ranking
# ─────────────────────────────────────────────────────────────────────────────

def section_heuristics(df: pd.DataFrame, console: Console, top_n: int) -> None:
    rows = []
    for pfx, label in HEURISTICS:
        time_col  = f"heur_{pfx}_time"
        sols_col  = f"heur_{pfx}_sols"
        execs_col = f"heur_{pfx}_execs"

        t    = _safe_sum(df, time_col)
        sols = _safe_sum(df, sols_col)
        execs = _safe_sum(df, execs_col)

        if t > 0 or sols > 0:
            rows.append({
                "label": label,
                "time":  t,
                "sols":  int(sols),
                "execs": int(execs),
                "sols_per_exec": sols / execs if execs > 0 else 0.0,
            })

    if not rows:
        console.print("[dim]No heuristic data found.[/dim]")
        return

    rows.sort(key=lambda r: r["time"], reverse=True)
    max_time = rows[0]["time"] if rows else 1.0

    table = Table(
        title=f"💡 Heuristics — top {min(top_n, len(rows))} by time",
        box=ROUNDED, header_style="bold cyan", show_lines=False,
    )
    table.add_column("Heuristic",    style="bold",    width=22)
    table.add_column("Total time",   justify="right", width=12)
    table.add_column("Bar",          width=20)
    table.add_column("Solutions",    justify="right", width=10)
    table.add_column("Executions",   justify="right", width=11)
    table.add_column("Sol/exec",     justify="right", width=9)

    for r in rows[:top_n]:
        frac = r["time"] / max_time if max_time > 0 else 0.0
        sol_color = "green" if r["sols"] > 0 else "dim"
        table.add_row(
            r["label"],
            fmt_time(r["time"]),
            f"[green]{bar(frac)}[/green]",
            f"[{sol_color}]{r['sols']:,}[/{sol_color}]",
            f"{r['execs']:,}",
            f"{r['sols_per_exec']:.3f}",
        )

    console.print(table)


# ─────────────────────────────────────────────────────────────────────────────
# Section 4: Worst instances
# ─────────────────────────────────────────────────────────────────────────────

def section_worst_instances(df: pd.DataFrame, console: Console, top_n: int) -> None:
    if "time" not in df.columns:
        return

    df = df.copy()
    df["_time"] = pd.to_numeric(df.get("time", pd.Series()), errors="coerce").fillna(0.0)
    df["_cut_time"]   = _safe_col(df, "cut_time")
    df["_lp_time"]    = _safe_col(df, "lp_seconds")
    df["_heur_time"]  = sum(_safe_col(df, f"heur_{pfx}_time") for pfx, _ in HEURISTICS)
    df["_cgraph_time"] = _safe_col(df, "cgraph_time")

    worst = df.nlargest(top_n, "_time")

    table = Table(
        title=f"🐢 Top-{top_n} slowest instances",
        box=ROUNDED, header_style="bold cyan", show_lines=False,
    )
    table.add_column("Instance",      style="bold", no_wrap=True)
    table.add_column("Total",         justify="right", width=9)
    table.add_column("LP",            justify="right", width=9)
    table.add_column("Cuts",          justify="right", width=9)
    table.add_column("Heuristics",    justify="right", width=11)
    table.add_column("CGraph",        justify="right", width=9)
    table.add_column("Nodes",         justify="right", width=10)

    for _, r in worst.iterrows():
        nodes_val = r.get("nodes", "")
        try:
            nodes_str = f"{int(float(nodes_val)):,}"
        except (ValueError, TypeError):
            nodes_str = str(nodes_val)

        table.add_row(
            str(r.get("instance", "?")),
            fmt_time(r["_time"]),
            fmt_time(r["_lp_time"]),
            fmt_time(r["_cut_time"]),
            fmt_time(r["_heur_time"]),
            fmt_time(r["_cgraph_time"]),
            nodes_str,
        )

    console.print(table)


# ─────────────────────────────────────────────────────────────────────────────
# Section 5: Wasted cut time
# ─────────────────────────────────────────────────────────────────────────────

def section_wasted_cuts(df: pd.DataFrame, console: Console, top_n: int) -> None:
    """
    Two sub-tables:
      A) Cut generators ranked by time-per-cut (lots of time, few cuts produced).
         Only generators that actually ran (calls > 0) are included.
      B) Per-instance top-N spotlight: for each generator with a high time-per-cut,
         the worst individual instances (most time spent, 0 cuts produced).
    """
    # ── A: generator-level waste ranking ─────────────────────────────────────
    gen_rows = []
    for pfx, label in CUT_GENERATORS:
        t     = _safe_sum(df, f"cut_{pfx}_time")
        cuts  = _safe_sum(df, f"cut_{pfx}_cuts")
        calls = _safe_sum(df, f"cut_{pfx}_calls")
        if t <= 0 or calls <= 0:
            continue
        time_per_cut = t / cuts if cuts > 0 else float("inf")
        zero_cut_instances = int(
            ((_safe_col(df, f"cut_{pfx}_calls") > 0) &
             (_safe_col(df, f"cut_{pfx}_cuts") == 0)).sum()
        )
        zero_cut_time = (
            _safe_col(df, f"cut_{pfx}_time")[
                (_safe_col(df, f"cut_{pfx}_calls") > 0) &
                (_safe_col(df, f"cut_{pfx}_cuts") == 0)
            ].sum()
        )
        gen_rows.append({
            "label":             label,
            "pfx":               pfx,
            "time":              t,
            "cuts":              int(cuts),
            "calls":             int(calls),
            "time_per_cut":      time_per_cut,
            "zero_cut_insts":    zero_cut_instances,
            "zero_cut_time":     zero_cut_time,
        })

    if not gen_rows:
        return

    # Sort by wasted time (time on instances where cuts=0)
    gen_rows.sort(key=lambda r: r["zero_cut_time"], reverse=True)
    max_wasted = gen_rows[0]["zero_cut_time"] if gen_rows else 1.0

    table = Table(
        title="⚠  Cut generators: time spent where no cuts were produced",
        box=ROUNDED, header_style="bold cyan", show_lines=False,
    )
    table.add_column("Generator",      style="bold",    width=22)
    table.add_column("Wasted time",    justify="right", width=13)
    table.add_column("Bar",            width=20)
    table.add_column("# Insts(0 cuts)", justify="right", width=14)
    table.add_column("Total time",     justify="right", width=12)
    table.add_column("Time/cut (s)",   justify="right", width=13)

    for r in gen_rows[:top_n]:
        wasted_frac = r["zero_cut_time"] / max_wasted if max_wasted > 0 else 0.0
        tpc = r["time_per_cut"]
        tpc_str = f"{tpc:.4f}" if tpc < float("inf") else "[bold red]∞[/bold red]"
        waste_color = "red" if r["zero_cut_time"] > 60 else "yellow" if r["zero_cut_time"] > 5 else "dim"
        table.add_row(
            r["label"],
            f"[{waste_color}]{fmt_time(r['zero_cut_time'])}[/{waste_color}]",
            f"[red]{bar(wasted_frac)}[/red]",
            f"{r['zero_cut_insts']:,}",
            fmt_time(r["time"]),
            tpc_str,
        )

    console.print(table)

    # ── B: per-instance spotlight for top wasted generators ──────────────────
    # Show top-3 generators by wasted time, each with their worst instances
    spotlight_gens = [r for r in gen_rows if r["zero_cut_time"] > 0][:3]
    for gen in spotlight_gens:
        pfx   = gen["pfx"]
        label = gen["label"]
        t_col = f"cut_{pfx}_time"
        c_col = f"cut_{pfx}_cuts"
        k_col = f"cut_{pfx}_calls"

        mask = (_safe_col(df, k_col) > 0) & (_safe_col(df, c_col) == 0)
        sub = df[mask].copy()
        if sub.empty:
            continue
        sub["_wt"] = _safe_col(sub, t_col)
        worst = sub.nlargest(min(top_n, 8), "_wt")

        inst_table = Table(
            title=f"  ↳ {label}: top instances with 0 cuts, most time spent",
            box=ROUNDED, header_style="bold red", show_lines=False,
        )
        inst_table.add_column("Instance",  style="bold", no_wrap=True)
        inst_table.add_column("Time (wasted)", justify="right", width=15)
        inst_table.add_column("Calls",        justify="right", width=8)
        inst_table.add_column("Total time",   justify="right", width=12)

        for _, r in worst.iterrows():
            total_t = pd.to_numeric(r.get("time", 0), errors="coerce") or 0.0
            inst_table.add_row(
                str(r.get("instance", "?")),
                f"[red]{fmt_time(r['_wt'])}[/red]",
                f"{int(_safe_col(df.loc[[r.name]], k_col).iloc[0]):,}",
                fmt_time(total_t),
            )
        console.print(inst_table)


# ─────────────────────────────────────────────────────────────────────────────
# Section 6: Wasted heuristic time
# ─────────────────────────────────────────────────────────────────────────────

def section_wasted_heuristics(df: pd.DataFrame, console: Console, top_n: int) -> None:
    """
    Two sub-tables:
      A) Heuristics ranked by time spent on instances where they found 0 solutions.
      B) Per-instance spotlight for the top wasted heuristics.
    """
    heur_rows = []
    for pfx, label in HEURISTICS:
        t     = _safe_sum(df, f"heur_{pfx}_time")
        sols  = _safe_sum(df, f"heur_{pfx}_sols")
        execs = _safe_sum(df, f"heur_{pfx}_execs")
        if t <= 0 or execs <= 0:
            continue

        zero_sol_mask = (
            (_safe_col(df, f"heur_{pfx}_execs") > 0) &
            (_safe_col(df, f"heur_{pfx}_sols") == 0)
        )
        zero_sol_time = _safe_col(df, f"heur_{pfx}_time")[zero_sol_mask].sum()
        zero_sol_insts = int(zero_sol_mask.sum())

        heur_rows.append({
            "label":          label,
            "pfx":            pfx,
            "time":           t,
            "sols":           int(sols),
            "execs":          int(execs),
            "zero_sol_time":  zero_sol_time,
            "zero_sol_insts": zero_sol_insts,
            "time_per_sol":   t / sols if sols > 0 else float("inf"),
        })

    if not heur_rows:
        return

    heur_rows.sort(key=lambda r: r["zero_sol_time"], reverse=True)
    max_wasted = heur_rows[0]["zero_sol_time"] if heur_rows else 1.0

    table = Table(
        title="⚠  Heuristics: time spent where no solution was found",
        box=ROUNDED, header_style="bold cyan", show_lines=False,
    )
    table.add_column("Heuristic",       style="bold",    width=22)
    table.add_column("Wasted time",     justify="right", width=13)
    table.add_column("Bar",             width=20)
    table.add_column("# Insts(0 sols)", justify="right", width=14)
    table.add_column("Total time",      justify="right", width=12)
    table.add_column("Time/sol (s)",    justify="right", width=13)

    for r in heur_rows[:top_n]:
        wasted_frac = r["zero_sol_time"] / max_wasted if max_wasted > 0 else 0.0
        tps = r["time_per_sol"]
        tps_str = fmt_time(tps) if tps < float("inf") else "[bold red]∞[/bold red]"
        waste_color = "red" if r["zero_sol_time"] > 60 else "yellow" if r["zero_sol_time"] > 5 else "dim"
        table.add_row(
            r["label"],
            f"[{waste_color}]{fmt_time(r['zero_sol_time'])}[/{waste_color}]",
            f"[red]{bar(wasted_frac)}[/red]",
            f"{r['zero_sol_insts']:,}",
            fmt_time(r["time"]),
            tps_str,
        )

    console.print(table)

    # ── Per-instance spotlight for top wasted heuristics ─────────────────────
    spotlight = [r for r in heur_rows if r["zero_sol_time"] > 0][:3]
    for h in spotlight:
        pfx   = h["pfx"]
        label = h["label"]
        t_col = f"heur_{pfx}_time"
        s_col = f"heur_{pfx}_sols"
        e_col = f"heur_{pfx}_execs"

        mask = (_safe_col(df, e_col) > 0) & (_safe_col(df, s_col) == 0)
        sub = df[mask].copy()
        if sub.empty:
            continue
        sub["_wt"] = _safe_col(sub, t_col)
        worst = sub.nlargest(min(top_n, 8), "_wt")

        inst_table = Table(
            title=f"  ↳ {label}: top instances with 0 solutions, most time spent",
            box=ROUNDED, header_style="bold red", show_lines=False,
        )
        inst_table.add_column("Instance",      style="bold", no_wrap=True)
        inst_table.add_column("Time (wasted)", justify="right", width=15)
        inst_table.add_column("Executions",    justify="right", width=12)
        inst_table.add_column("Total time",    justify="right", width=12)

        for _, r in worst.iterrows():
            total_t = pd.to_numeric(r.get("time", 0), errors="coerce") or 0.0
            execs_val = int(_safe_col(df.loc[[r.name]], e_col).iloc[0])
            inst_table.add_row(
                str(r.get("instance", "?")),
                f"[red]{fmt_time(r['_wt'])}[/red]",
                f"{execs_val:,}",
                fmt_time(total_t),
            )
        console.print(inst_table)


# ─────────────────────────────────────────────────────────────────────────────
# Section 7: Per-instance cut breakdown (optional, verbose)
# ─────────────────────────────────────────────────────────────────────────────

def section_cut_breakdown_by_instance(df: pd.DataFrame, console: Console, top_n: int) -> None:
    """Top-N instances by total cut_time with cut-type breakdown."""
    if "cut_time" not in df.columns:
        return

    df = df.copy()
    df["_cut_time"] = _safe_col(df, "cut_time")
    worst = df.nlargest(top_n, "_cut_time")

    # Find cut generators that actually have data
    active_gens = [
        (pfx, label) for pfx, label in CUT_GENERATORS
        if _safe_sum(df, f"cut_{pfx}_time") > 0
    ]
    if not active_gens:
        return

    table = Table(
        title=f"✂  Cut time breakdown — top-{top_n} instances by cut time",
        box=ROUNDED, header_style="bold cyan", show_lines=False,
    )
    table.add_column("Instance", style="bold", no_wrap=True, width=22)
    table.add_column("Total",    justify="right", width=9)
    for _, label in active_gens:
        table.add_column(label[:12], justify="right", width=10)

    for _, r in worst.iterrows():
        total_ct = r["_cut_time"]
        if total_ct <= 0:
            break
        row_data = [str(r.get("instance", "?")), fmt_time(total_ct)]
        for pfx, _ in active_gens:
            t = pd.to_numeric(r.get(f"cut_{pfx}_time", 0), errors="coerce")
            t = 0.0 if t != t else float(t)
            pct = t / total_ct if total_ct > 0 else 0.0
            if t > 0:
                row_data.append(f"{fmt_time(t)} ({pct*100:.0f}%)")
            else:
                row_data.append("[dim]–[/dim]")
        table.add_row(*row_data)

    console.print(table)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Detailed time/efficiency analysis from MIPster stats CSV files."
    )
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("--outdir",    help="Experiment output directory (looks for stats.csv inside)")
    grp.add_argument("--statsfile", help="Direct path to a stats.csv file")
    parser.add_argument("--top",    type=int, default=15, metavar="N",
                        help="Number of top entries to show per table [default: 15]")
    parser.add_argument("--plain",  action="store_true", help="Suppress colour output")
    parser.add_argument("--cut-breakdown", action="store_true",
                        help="Also show per-instance cut-type breakdown table")
    args = parser.parse_args()

    if args.statsfile:
        stats_path = Path(args.statsfile)
    else:
        outdir = Path(args.outdir)
        # Try both stats.csv (combined) and individual *.stats.csv files
        stats_path = outdir / "stats.csv"
        if not stats_path.exists():
            # Combine per-instance stats files
            parts = sorted(outdir.glob("*.stats.csv"))
            if not parts:
                print(f"Error: no stats.csv or *.stats.csv files found in {outdir}",
                      file=sys.stderr)
                sys.exit(1)
            combined = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
            stats_path = outdir / "stats.csv"
            combined.to_csv(stats_path, index=False)
            print(f"Combined {len(parts)} stats files → {stats_path}")

    if not stats_path.exists():
        print(f"Error: {stats_path} not found", file=sys.stderr)
        sys.exit(1)

    df = load_stats(stats_path)
    n = len(df)

    console = Console(record=True, force_terminal=not args.plain, width=160)

    label = Path(args.outdir).name if args.outdir else stats_path.name
    console.print()
    console.print(Panel(
        f"[bold]Stats analysis[/bold] — [cyan]{label}[/cyan]\n"
        f"[dim]{n} instances   {stats_path}[/dim]",
        border_style="cyan",
    ))

    section_time_budget(df, console)
    console.print()
    section_cuts(df, console, args.top)
    console.print()
    section_heuristics(df, console, args.top)
    console.print()
    section_wasted_cuts(df, console, args.top)
    console.print()
    section_wasted_heuristics(df, console, args.top)
    console.print()
    section_worst_instances(df, console, args.top)

    if args.cut_breakdown:
        console.print()
        section_cut_breakdown_by_instance(df, console, args.top)

    # Save plain-text copy alongside stats.csv
    if args.outdir:
        out_txt = Path(args.outdir) / "stats_analysis.txt"
        plain_text = console.export_text()
        with open(out_txt, "w") as fh:
            fh.write(plain_text)
        console.print(f"\n[dim]Analysis saved to: {out_txt}[/dim]")


if __name__ == "__main__":
    main()
