#!/usr/bin/env python3
"""Render a rich summary table from a run_experiments.sh output directory.

Includes cut and heuristic analysis from per-instance .stats.csv files.

Usage:
    python3 scripts/summarize_results.py --outdir exp_results/myrun/
    python3 scripts/summarize_results.py --outdir exp_results/myrun/ --plain
    python3 scripts/summarize_results.py --outdir exp_results/myrun/ --no-stats

Cost system (lower is better):
    Solved to optimality / proven infeasible  →    0 pts
    Timeout with feasible solution            →  gap% (0–100)
    No feasible solution found                →  200 pts  (sentinel)
    Overtime (killed past wallclock)          →  300 pts  (sentinel)
    Wrong result                              →  200 pts
"""
import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

from rich.box import ROUNDED, SIMPLE_HEAD
from rich.console import Console
from rich.panel import Panel
from rich.rule import Rule
from rich.table import Table
from rich.text import Text


# ── Cost / gap helpers ────────────────────────────────────────────────────────

def _is_infeasible(status: str) -> bool:
    return ("INFEASIBLE" in status or "(inf)" in status) and "WRONG" not in status


def parse_gap_field(text: str):
    """Parse a gap string like '7.9%', '>100%', '-' → float or None."""
    text = str(text).strip()
    if not text or text == "-":
        return None
    if text.startswith(">"):
        return 100.0
    if text.endswith("%"):
        text = text[:-1]
    try:
        return min(float(text), 100.0)
    except ValueError:
        return None


def compute_cost(row: dict) -> float:
    """Return unified cost score (0=optimal, gap%=timeout, 200=no sol, 300=overtime)."""
    status = row.get("status", "")
    if status == "OVERTIME":
        return 300.0
    if "WRONG" in status:
        return 200.0
    if _is_infeasible(status):
        return 0.0
    if status.startswith("SOLVED"):
        return 0.0

    gf = parse_gap_field(row.get("gap_field", ""))
    if gf is not None:
        return gf
    m = re.search(r"gap=([0-9.]+)%", status)
    if m:
        return min(float(m.group(1)), 100.0)
    if "no_sol" in status:
        return 200.0

    try:
        obj  = float(row.get("objective", "nan"))
        dual = float(row.get("dual_bound", "nan"))
        if obj == obj and dual == dual:
            if abs(dual) > 1e-10:
                return min(abs(obj - dual) / abs(dual) * 100, 100.0)
            return 0.0 if abs(obj) <= 1e-10 else 100.0
        if obj == obj:
            return 100.0
    except (ValueError, TypeError):
        pass

    return 200.0


# ── Status display helpers ────────────────────────────────────────────────────

def status_style(status: str) -> str:
    if status.startswith("SOLVED") or _is_infeasible(status):
        return "green"
    if re.match(r"^TIMEOUT\(gap=", status):
        return "yellow"
    if status in ("TIMEOUT", "TIMEOUT(no_sol)"):
        return "dark_orange"
    if status in ("NO_OBJ", "OVERTIME", "MISSING"):
        return "dim"
    return "bold red"


def status_icon(status: str) -> str:
    if status.startswith("SOLVED") or _is_infeasible(status):
        return "✓"
    if re.match(r"^TIMEOUT", status):
        return "⏱"
    if status in ("NO_OBJ", "OVERTIME", "MISSING"):
        return "–"
    return "✗"


def result_category(status: str) -> str:
    if status.startswith("SOLVED") or _is_infeasible(status):
        return "pass"
    if re.match(r"^TIMEOUT\(gap=", status) or status == "TIMEOUT":
        return "timeout_sol"
    if status == "TIMEOUT(no_sol)":
        return "timeout_nosol"
    if status == "OVERTIME":
        return "overtime"
    return "fail"


_CAT_ORDER = {"pass": 0, "timeout_sol": 1, "timeout_nosol": 2, "overtime": 3, "fail": 4}


def fmt_num(s: str) -> str:
    if not s or s == "-":
        return "-"
    try:
        v = float(s)
        return f"{v:,.8g}"
    except ValueError:
        return s


def fmt_cost(cost: float) -> str:
    if cost == 0.0:
        return "[green]0[/green]"
    if cost >= 300.0:
        return "[dim]OVT[/dim]"
    if cost >= 200.0:
        return "[bold red]NO SOL[/bold red]"
    return f"[yellow]{cost:.1f}%[/yellow]"


def clean_status_label(status: str) -> str:
    if re.match(r"^TIMEOUT\(gap=", status):
        return "TIMEOUT"
    return status


# ── Stats loading ─────────────────────────────────────────────────────────────

def load_stats_rows(outdir: Path, instances: list) -> list:
    """Load first row of each .stats.csv; return list of dicts."""
    rows = []
    for inst in instances:
        path = outdir / f"{inst}.stats.csv"
        if path.exists():
            with open(path, newline="") as fh:
                data = list(csv.DictReader(fh))
                if data:
                    rows.append(data[0])
    return rows


def _flt(val, default=0.0) -> float:
    try:
        return float(val or 0)
    except (ValueError, TypeError):
        return default


def _int(val, default=0) -> int:
    try:
        return int(float(val or 0))
    except (ValueError, TypeError):
        return default


# ── Stats aggregation ─────────────────────────────────────────────────────────

def aggregate_cuts(stats_rows: list) -> dict:
    """Return {name: {cuts, time, calls, active}} for all cut types."""
    agg = defaultdict(lambda: dict(cuts=0, time=0.0, calls=0, active=0))
    for row in stats_rows:
        for key in row:
            if key.startswith("cut_") and key.endswith("_cuts"):
                name = key[4:-5]  # strip 'cut_' prefix and '_cuts' suffix
                n = _int(row[key])
                agg[name]["cuts"] += n
                if n > 0:
                    agg[name]["active"] += 1
                agg[name]["time"]  += _flt(row.get(f"cut_{name}_time"))
                agg[name]["calls"] += _int(row.get(f"cut_{name}_calls"))
    return dict(agg)


def aggregate_heuristics(stats_rows: list) -> dict:
    """Return {name: {execs, time, sols, active}} for all heuristic types."""
    agg = defaultdict(lambda: dict(execs=0, time=0.0, sols=0, active=0))
    for row in stats_rows:
        for key in row:
            if key.startswith("heur_") and key.endswith("_execs"):
                name = key[5:-6]  # strip 'heur_' and '_execs'
                n = _int(row[key])
                agg[name]["execs"] += n
                if n > 0:
                    agg[name]["active"] += 1
                agg[name]["time"] += _flt(row.get(f"heur_{name}_time"))
                agg[name]["sols"] += _int(row.get(f"heur_{name}_sols"))
    return dict(agg)


def aggregate_time_budget(stats_rows: list) -> dict:
    """Sum key time fields across all stats rows."""
    budget = dict(elapsed=0.0, lp=0.0, cuts=0.0, cgraph=0.0, heur=0.0,
                  nodes=0, iterations=0)
    for row in stats_rows:
        budget["elapsed"]    += _flt(row.get("elapsed"))
        budget["lp"]         += _flt(row.get("lp_seconds"))
        budget["cuts"]       += _flt(row.get("cut_time"))
        budget["cgraph"]     += _flt(row.get("cgraph_time"))
        budget["nodes"]      += _int(row.get("nodes"))
        budget["iterations"] += _int(row.get("iterations"))
        for key in row:
            if key.startswith("heur_") and key.endswith("_time"):
                budget["heur"] += _flt(row[key])
    return budget


# ── Pretty name maps ──────────────────────────────────────────────────────────

_CUT_NAMES = {
    "mixedintegerrounding2": "MIR2",
    "twomircuts":            "TwoMIR",
    "twomircutsl1":          "TwoMIR-L1",
    "twomircutsl2":          "TwoMIR-L2",
    "gomoryl1":              "Gomory-L1",
    "gomoryl2":              "Gomory-L2",
    "gomory_2":              "Gomory-root",
    "gomory":                "Gomory",
    "probing":               "Probing",
    "knapsack":              "Knapsack",
    "reduce_and_split":      "Reduce&Split",
    "reduce_and_split_2":    "Reduce&Split-root",
    "clique":                "Clique",
    "oddwheel":              "OddWheel",
    "flowcover":             "FlowCover",
    "liftandproject":        "Lift&Project",
    "residualcapacity":      "ResidualCap",
    "zerohalf":              "ZeroHalf",
    "stored":                "Stored",
}

_HEUR_NAMES = {
    "feasibility_pump":        "FeasPump",
    "feasibilityjump":         "FeasJump",
    "rounding":                "Rounding",
    "random_rounding":         "RandRound",
    "combine_solutions":       "Combine",
    "greedy_cover":            "Greedy-Cover",
    "greedy_equality":         "Greedy-Eq",
    "dynamic_pass_thru":       "DynPassThru",
    "linked":                  "Linked",
    "partial_solution_given":  "PartialSol",
    "dantzig_wolfe_expansion": "DW-Expand",
    "rins":                    "RINS",
    "rens":                    "RENS",
    "rensdj":                  "RENS-DJ",
    "rensub":                  "RENS-UB",
    "vnd":                     "VND",
    "naive":                   "Naive",
    "diveany":                 "Dive-Any",
    "divecoefficient":         "Dive-Coeff",
    "divefractional":          "Dive-Frac",
    "diveguided":              "Dive-Guided",
    "divelinesearch":          "Dive-LineSearch",
    "divepseudocost":          "Dive-Pseudo",
    "divevectorlength":        "Dive-VecLen",
    "multiple_root_solvers":   "MultiRoot",
}


def pretty_cut(name: str) -> str:
    return _CUT_NAMES.get(name, name)


def pretty_heur(name: str) -> str:
    return _HEUR_NAMES.get(name, name)


# ── Sparkbar helper ───────────────────────────────────────────────────────────

def spark_bar(pct: float, width: int = 20) -> str:
    """Return a UTF-8 block-element progress bar for a percentage 0–100."""
    filled = int(round(pct / 100 * width))
    filled = max(0, min(filled, width))
    return "█" * filled + "░" * (width - filled)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render a rich summary for a run_experiments.sh result directory."
    )
    parser.add_argument("--outdir", required=True, help="Experiment output directory")
    parser.add_argument("--summary", help="Override path to summary.tsv")
    parser.add_argument("--report",  help="Override path to report.txt (plain-text append)")
    parser.add_argument("--plain",    action="store_true", help="Suppress colour output")
    parser.add_argument("--no-stats", action="store_true", help="Skip .stats.csv analysis")
    args = parser.parse_args()

    outdir       = Path(args.outdir)
    summary_path = Path(args.summary) if args.summary else outdir / "summary.tsv"
    report_path  = Path(args.report)  if args.report  else outdir / "report.txt"

    if not summary_path.exists():
        print(f"Error: {summary_path} not found", file=sys.stderr)
        sys.exit(1)

    # ── Read run metadata from existing report.txt header ────────────────────
    metadata: dict = {}
    if report_path.exists():
        with open(report_path) as fh:
            for line in fh:
                line = line.rstrip()
                if not line or line.startswith("═") or line.startswith("─"):
                    break
                if ":" in line:
                    key, _, val = line.partition(":")
                    metadata[key.strip()] = val.strip()

    # ── Read summary.tsv ──────────────────────────────────────────────────────
    rows: list = []
    with open(summary_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append({k.strip(): v.strip() for k, v in row.items()})

    if not rows:
        print("No data in summary.tsv", file=sys.stderr)
        sys.exit(1)

    has_dual_col = "dual_bound" in rows[0]

    for r in rows:
        r["_cost"] = compute_cost(r)

    rows.sort(key=lambda r: (_CAT_ORDER.get(result_category(r["status"]), 9),
                              r["instance"]))

    # ── Load stats CSVs ───────────────────────────────────────────────────────
    stats_rows: list = []
    if not args.no_stats:
        instances = [r["instance"] for r in rows]
        stats_rows = load_stats_rows(outdir, instances)

    # ── Console setup ─────────────────────────────────────────────────────────
    console = Console(record=True, force_terminal=not args.plain, width=160)

    # ── Per-instance results table ────────────────────────────────────────────
    table = Table(
        title=f"[bold]MIPster Experiment Results[/bold] — {outdir.name}",
        box=ROUNDED,
        header_style="bold cyan",
        show_lines=False,
        expand=True,
    )
    table.add_column("#",         justify="right", style="dim", no_wrap=True, width=4)
    table.add_column("Instance",  style="bold",    no_wrap=True)
    table.add_column("Status",    no_wrap=True,    width=24)
    table.add_column("Objective", justify="right", width=18)
    if has_dual_col:
        table.add_column("Dual bound", justify="right", width=18)
    table.add_column("Reference", justify="right", width=18)
    table.add_column("Cost",      justify="right", width=10)
    table.add_column("Time (s)",  justify="right", width=9)

    counts_by_status: dict = {}
    pass_count = timeout_sol = timeout_nosol = overtime_count = fail_count = 0
    total_cost = 0.0

    for idx, row in enumerate(rows, start=1):
        inst    = row["instance"]
        status  = row["status"]
        obj     = row.get("objective", "-")
        ref     = row.get("expected", "-")
        elapsed = row.get("elapsed_s", "-")
        dual    = row.get("dual_bound", "-") if has_dual_col else None
        cost    = row["_cost"]
        cat     = result_category(status)
        style   = status_style(status)
        icon    = status_icon(status)

        counts_by_status[status] = counts_by_status.get(status, 0) + 1
        if   cat == "pass":          pass_count     += 1
        elif cat == "timeout_sol":   timeout_sol    += 1
        elif cat == "timeout_nosol": timeout_nosol  += 1
        elif cat == "overtime":      overtime_count += 1
        else:                        fail_count     += 1
        total_cost += cost

        status_cell = Text(f"{icon} {clean_status_label(status)}", style=style)
        row_style   = "dim" if cat in ("timeout_sol", "timeout_nosol", "overtime") else ""

        row_data = [str(idx), inst, status_cell, fmt_num(obj)]
        if has_dual_col:
            row_data.append(fmt_num(dual))
        row_data += [fmt_num(ref), fmt_cost(cost), elapsed]
        table.add_row(*row_data, style=row_style)

    total = len(rows)
    console.print()
    console.print(table)

    # ── Summary panel ─────────────────────────────────────────────────────────
    lines = []
    lines.append(f"[bold]Total:[/bold]              {total} instances")
    lines.append(f"[green]✓  Solved:[/green]           {pass_count}  ({pass_count / total * 100:.1f}%)")
    if timeout_sol:
        lines.append(f"[yellow]⏱  Timeout+sol:[/yellow]     {timeout_sol}  ({timeout_sol / total * 100:.1f}%)")
    if timeout_nosol:
        lines.append(f"[dark_orange]⏱  Timeout(no sol):[/dark_orange]  {timeout_nosol}  ({timeout_nosol / total * 100:.1f}%)")
    if overtime_count:
        lines.append(f"[dim]—  Overtime:[/dim]          {overtime_count}  ({overtime_count / total * 100:.1f}%)")
    if fail_count:
        lines.append(f"[bold red]✗  Failed:[/bold red]          {fail_count}  ({fail_count / total * 100:.1f}%)")
    lines.append("")
    lines.append(f"[bold]Total cost:[/bold]         {total_cost:.1f}  (avg {total_cost / total:.2f} per instance)")
    lines.append(f"[dim]  0 = optimal · gap% = timeout with sol · 200 = no sol · 300 = overtime[/dim]")
    lines.append("")
    lines.append("[bold]By status:[/bold]")
    for st, cnt in sorted(counts_by_status.items(),
                          key=lambda x: (_CAT_ORDER.get(result_category(x[0]), 9), -x[1])):
        color = status_style(st)
        bar   = "█" * min(cnt, 40)
        lines.append(f"  [{color}]{st:<34}[/{color}]  {cnt:>4}  [dim]{bar}[/dim]")

    console.print(Panel("\n".join(lines), title="Summary", border_style="cyan"))

    # ── Stats sections ────────────────────────────────────────────────────────
    if stats_rows:
        n_stats = len(stats_rows)
        budget  = aggregate_time_budget(stats_rows)
        total_elapsed = budget["elapsed"]

        # ── Time budget panel ─────────────────────────────────────────────────
        tb_lines = [
            f"[bold]Instances with stats:[/bold]  {n_stats} / {total}"
            + (f"  [dim](missing: {total - n_stats})[/dim]" if n_stats < total else ""),
            f"[bold]Total elapsed:[/bold]         {total_elapsed:,.1f} s"
            + f"  [dim](avg {total_elapsed / n_stats:.1f} s/inst)[/dim]",
            f"[bold]Total nodes:[/bold]           {budget['nodes']:,}"
            + f"  [dim](avg {budget['nodes'] / n_stats:,.0f}/inst)[/dim]",
            f"[bold]Total LP iters:[/bold]        {budget['iterations']:,}"
            + f"  [dim](avg {budget['iterations'] / n_stats:,.0f}/inst)[/dim]",
            "",
            "[bold]Time breakdown (summed across all instances):[/bold]",
        ]

        budget_items = [
            ("LP solves",          budget["lp"],     "cyan"),
            ("Cut generation",     budget["cuts"],   "yellow"),
            ("Conflict graph",     budget["cgraph"], "blue"),
            ("Heuristics",         budget["heur"],   "magenta"),
        ]
        accounted = sum(v for _, v, _ in budget_items)
        other     = max(0.0, total_elapsed - accounted)
        budget_items.append(("B&B / other overhead", other, "dim"))

        for label, t, color in budget_items:
            pct = t / total_elapsed * 100 if total_elapsed > 0 else 0
            bar = spark_bar(pct, 24)
            tb_lines.append(
                f"  [{color}]{label:<22}[/{color}]"
                f"  {t:>10,.1f} s"
                f"  {pct:>5.1f}%"
                f"  [dim]{bar}[/dim]"
            )

        console.print(Panel("\n".join(tb_lines), title="⏱  Time Budget", border_style="blue"))

        # ── Cut statistics table ──────────────────────────────────────────────
        cut_agg      = aggregate_cuts(stats_rows)
        active_cuts  = {k: v for k, v in cut_agg.items() if v["cuts"] > 0}
        total_cuts   = sum(v["cuts"]  for v in active_cuts.values())
        total_ctime  = sum(v["time"]  for v in active_cuts.values())

        if active_cuts:
            cut_table = Table(
                title=(
                    f"[bold yellow]✂  Cut Statistics[/bold yellow]"
                    f"  [dim]({len(active_cuts)} active types"
                    f" · {total_cuts:,} cuts total"
                    f" · {total_ctime:,.1f} s total)[/dim]"
                ),
                box=SIMPLE_HEAD,
                header_style="bold yellow",
                show_lines=False,
                expand=True,
            )
            cut_table.add_column("Cut Type",         no_wrap=True, min_width=18)
            cut_table.add_column("Cuts",              justify="right", width=10)
            cut_table.add_column("Cuts %",            justify="right", width=8)
            cut_table.add_column("Time (s)",          justify="right", width=10)
            cut_table.add_column("Time %",            justify="right", width=8)
            cut_table.add_column("Calls",             justify="right", width=10)
            cut_table.add_column("Active\ninsts",     justify="right", width=8)
            cut_table.add_column("Avg cuts\n/active", justify="right", width=10)
            cut_table.add_column("Cut/call\nrate",    justify="right", width=10)

            for name, data in sorted(active_cuts.items(), key=lambda x: -x[1]["time"]):
                time_pct  = data["time"]  / total_ctime * 100 if total_ctime  > 0 else 0
                cuts_pct  = data["cuts"]  / total_cuts  * 100 if total_cuts   > 0 else 0
                avg_cuts  = data["cuts"]  / data["active"] if data["active"] > 0 else 0
                call_rate = data["cuts"]  / data["calls"]  * 100 if data["calls"] > 0 else 0

                if time_pct > 15:
                    sty = "bold yellow"
                elif time_pct > 3:
                    sty = "white"
                else:
                    sty = "dim"

                cut_table.add_row(
                    f"[{sty}]{pretty_cut(name)}[/{sty}]",
                    f"[{sty}]{data['cuts']:,}[/{sty}]",
                    f"[{sty}]{cuts_pct:.1f}%[/{sty}]",
                    f"[{sty}]{data['time']:,.1f}[/{sty}]",
                    f"[{sty}]{time_pct:.1f}%[/{sty}]",
                    f"[{sty}]{data['calls']:,}[/{sty}]",
                    f"[{sty}]{data['active']}[/{sty}]",
                    f"[{sty}]{avg_cuts:,.0f}[/{sty}]",
                    f"[{sty}]{call_rate:.1f}%[/{sty}]",
                )

            console.print(cut_table)
        else:
            console.print(Panel(
                "[dim]No cuts found in .stats.csv files.[/dim]",
                title="✂  Cut Statistics", border_style="yellow"
            ))

        # ── Heuristic statistics table ────────────────────────────────────────
        heur_agg       = aggregate_heuristics(stats_rows)
        active_heurs   = {k: v for k, v in heur_agg.items() if v["execs"] > 0}
        total_htime    = sum(v["time"] for v in active_heurs.values())
        total_heur_sols = sum(v["sols"] for v in active_heurs.values())
        total_execs    = sum(v["execs"] for v in active_heurs.values())

        if active_heurs:
            heur_table = Table(
                title=(
                    f"[bold magenta]🎲  Heuristic Statistics[/bold magenta]"
                    f"  [dim]({len(active_heurs)} active"
                    f" · {total_execs:,} execs"
                    f" · {total_heur_sols:,} solutions found"
                    f" · {total_htime:,.1f} s total)[/dim]"
                ),
                box=SIMPLE_HEAD,
                header_style="bold magenta",
                show_lines=False,
                expand=True,
            )
            heur_table.add_column("Heuristic",         no_wrap=True, min_width=18)
            heur_table.add_column("Execs",              justify="right", width=9)
            heur_table.add_column("Time (s)",           justify="right", width=10)
            heur_table.add_column("Time %",             justify="right", width=8)
            heur_table.add_column("Solutions",          justify="right", width=10)
            heur_table.add_column("Sol %\nof total",    justify="right", width=10)
            heur_table.add_column("Active\ninsts",      justify="right", width=8)
            heur_table.add_column("Sol/exec\nrate",     justify="right", width=10)

            for name, data in sorted(active_heurs.items(), key=lambda x: -x[1]["time"]):
                time_pct  = data["time"] / total_htime     * 100 if total_htime     > 0 else 0
                sol_pct   = data["sols"] / total_heur_sols * 100 if total_heur_sols > 0 else 0
                sol_rate  = data["sols"] / data["execs"]   * 100 if data["execs"]   > 0 else 0

                if data["sols"] > 0 and time_pct > 5:
                    sty = "bold magenta"
                elif data["sols"] > 0:
                    sty = "magenta"
                elif time_pct > 5:
                    sty = "white"
                else:
                    sty = "dim"

                heur_table.add_row(
                    f"[{sty}]{pretty_heur(name)}[/{sty}]",
                    f"[{sty}]{data['execs']:,}[/{sty}]",
                    f"[{sty}]{data['time']:,.2f}[/{sty}]",
                    f"[{sty}]{time_pct:.1f}%[/{sty}]",
                    f"[{sty}]{data['sols']:,}[/{sty}]",
                    f"[{sty}]{sol_pct:.1f}%[/{sty}]",
                    f"[{sty}]{data['active']}[/{sty}]",
                    f"[{sty}]{sol_rate:.2f}%[/{sty}]",
                )

            console.print(heur_table)
        else:
            console.print(Panel(
                "[dim]No heuristic activity found in .stats.csv files.[/dim]",
                title="🎲  Heuristic Statistics", border_style="magenta"
            ))

    elif not args.no_stats:
        console.print(Panel(
            f"[dim]No .stats.csv files found in {outdir}[/dim]",
            title="Stats", border_style="dim"
        ))

    # ── Run configuration panel ───────────────────────────────────────────────
    if metadata:
        meta_lines = [f"[bold]{k}:[/bold] {v}" for k, v in metadata.items()]
        console.print(Panel("\n".join(meta_lines), title="Run configuration",
                            border_style="blue"))

    # ── Failed / wrong instances panel ───────────────────────────────────────
    failed = [r for r in rows if result_category(r["status"]) == "fail"]
    if failed:
        err_lines = []
        for r in failed:
            line = f"[bold red]{r['instance']}[/bold red]  →  {r['status']}"
            err_path = outdir / f"{r['instance']}.err"
            if err_path.exists():
                line += f"  [dim](see {err_path})[/dim]"
            err_lines.append(line)
        console.print(Panel("\n".join(err_lines), title="Failed instances",
                            border_style="red"))

    # ── Footer ────────────────────────────────────────────────────────────────
    console.print(f"\n[dim]Logs:        {outdir}/[/dim]")
    console.print(f"[dim]Summary TSV: {summary_path}[/dim]")
    if stats_rows:
        console.print(f"[dim]Stats CSVs:  {len(stats_rows)} files loaded[/dim]")

    # ── Append plain-text version to report.txt ───────────────────────────────
    plain_text = console.export_text()
    with open(report_path, "a") as fh:
        fh.write("\n" + plain_text)


if __name__ == "__main__":
    main()
