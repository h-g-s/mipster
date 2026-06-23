#!/usr/bin/env python3
"""Render a rich summary table from a run_experiments.sh output directory.

Usage:
    python3 summarize_results.py --outdir exp_results/myrun/
    python3 summarize_results.py --outdir exp_results/myrun/ --plain

Cost system (lower is better, matches compare_multi_experiments.py):
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
from pathlib import Path

from rich.box import ROUNDED
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text


# ── Cost / gap helpers (aligned with compare_multi_experiments.py) ────────────

def _is_infeasible(status: str) -> bool:
    return ("INFEASIBLE" in status or "(inf)" in status) and "WRONG" not in status


def parse_gap_field(text: str):
    """Parse a gap_field string like '7.9%', '>100%', '-' → float or None."""
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
    """Return unified cost score (0=optimal, 200=no sol, 300=overtime)."""
    status = row.get("status", "")
    if status == "OVERTIME":
        return 300.0
    if "WRONG" in status:
        return 200.0
    if _is_infeasible(status):
        return 0.0
    if status.startswith("SOLVED"):
        return 0.0

    # Timeout with solution: use gap_field, then status-embedded gap
    gf = parse_gap_field(row.get("gap_field", ""))
    if gf is not None:
        return gf
    m = re.search(r"gap=([0-9.]+)%", status)
    if m:
        return min(float(m.group(1)), 100.0)

    # No solution found
    if "no_sol" in status:
        return 200.0

    # Try to compute from dual/primal bounds
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
    if status.startswith("SOLVED"):
        return "green"
    if _is_infeasible(status):
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


def extract_gap_from_status(status: str) -> str:
    m = re.match(r"TIMEOUT\(gap=(.+)\)", status)
    return m.group(1) if m else ""


def clean_status_label(status: str) -> str:
    if re.match(r"^TIMEOUT\(gap=", status):
        return "TIMEOUT"
    return status


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render rich summary table for a run_experiments.sh result directory."
    )
    parser.add_argument("--outdir", required=True, help="Experiment output directory")
    parser.add_argument("--summary", help="Override path to summary.tsv")
    parser.add_argument("--report", help="Override path to report.txt (plain-text append)")
    parser.add_argument("--plain", action="store_true", help="Suppress colour output")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    summary_path = Path(args.summary) if args.summary else outdir / "summary.tsv"
    report_path  = Path(args.report)  if args.report  else outdir / "report.txt"

    if not summary_path.exists():
        print(f"Error: {summary_path} not found", file=sys.stderr)
        sys.exit(1)

    # ── Read run metadata from existing report.txt header ────────────────────
    metadata: dict[str, str] = {}
    if report_path.exists():
        with open(report_path) as fh:
            for line in fh:
                line = line.rstrip()
                if not line or line.startswith("═") or line.startswith("─"):
                    break
                if ":" in line:
                    key, _, val = line.partition(":")
                    metadata[key.strip()] = val.strip()

    # ── Read summary.tsv ─────────────────────────────────────────────────────
    rows: list[dict] = []
    with open(summary_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append({k.strip(): v.strip() for k, v in row.items()})

    if not rows:
        print("No data in summary.tsv", file=sys.stderr)
        sys.exit(1)

    has_gap_col    = "gap"        in rows[0]
    has_dual_col   = "dual_bound" in rows[0]
    has_gapfield   = "gap_field"  in rows[0]

    # Compute cost for each row
    for r in rows:
        r["_cost"] = compute_cost(r)

    # Sort: pass → timeout with solution → timeout no sol → overtime → fail; then by name
    rows.sort(key=lambda r: (_CAT_ORDER.get(result_category(r["status"]), 9),
                             r["instance"]))

    # ── Build rich table ─────────────────────────────────────────────────────
    console = Console(record=True, force_terminal=not args.plain, width=160)

    table = Table(
        title=f"[bold]MIPster Experiment Results[/bold] — {outdir.name}",
        box=ROUNDED,
        header_style="bold cyan",
        show_lines=False,
        expand=True,
    )
    table.add_column("#",          justify="right",  style="dim",  no_wrap=True, width=4)
    table.add_column("Instance",   style="bold",     no_wrap=True)
    table.add_column("Status",     no_wrap=True,     width=24)
    table.add_column("Objective",  justify="right",  width=18)
    if has_dual_col:
        table.add_column("Dual bound", justify="right", width=18)
    table.add_column("Reference",  justify="right",  width=18)
    table.add_column("Cost",       justify="right",  width=10)
    table.add_column("Time (s)",   justify="right",  width=9)

    counts_by_status: dict[str, int] = {}
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

        cat   = result_category(status)
        style = status_style(status)
        icon  = status_icon(status)

        counts_by_status[status] = counts_by_status.get(status, 0) + 1
        if   cat == "pass":          pass_count      += 1
        elif cat == "timeout_sol":   timeout_sol     += 1
        elif cat == "timeout_nosol": timeout_nosol   += 1
        elif cat == "overtime":      overtime_count  += 1
        else:                        fail_count      += 1
        total_cost += cost

        status_cell = Text(f"{icon} {clean_status_label(status)}", style=style)

        row_style = "dim" if cat in ("timeout_sol", "timeout_nosol", "overtime") else ""

        row_data = [str(idx), inst, status_cell, fmt_num(obj)]
        if has_dual_col:
            row_data.append(fmt_num(dual))
        row_data += [fmt_num(ref), fmt_cost(cost), elapsed]

        table.add_row(*row_data, style=row_style)

    total = len(rows)
    console.print()
    console.print(table)

    # ── Summary panel ─────────────────────────────────────────────────────────
    lines: list[str] = []
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
    lines.append(f"[bold]  0=optimal  gap%=timeout  200=no sol  300=overtime[/bold]")
    lines.append("")
    lines.append("[bold]By status:[/bold]")
    for st, cnt in sorted(counts_by_status.items(),
                          key=lambda x: (_CAT_ORDER.get(result_category(x[0]), 9), -x[1])):
        color = status_style(st)
        bar = "█" * min(cnt, 40)
        lines.append(f"  [{color}]{st:<32}[/{color}]  {cnt:>4}  [dim]{bar}[/dim]")

    console.print(Panel("\n".join(lines), title="Summary", border_style="cyan"))

    # ── Run configuration panel ───────────────────────────────────────────────
    if metadata:
        meta_lines = [f"[bold]{k}:[/bold] {v}" for k, v in metadata.items()]
        console.print(Panel("\n".join(meta_lines), title="Run configuration",
                            border_style="blue"))

    # ── Failed / wrong instances panel ───────────────────────────────────────
    failed = [r for r in rows if result_category(r["status"]) == "fail"]
    if failed:
        err_lines: list[str] = []
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

    # ── Append plain-text version to report.txt ───────────────────────────────
    plain_text = console.export_text()
    with open(report_path, "a") as fh:
        fh.write("\n" + plain_text)


if __name__ == "__main__":
    main()
