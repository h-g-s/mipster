/*
 * mipster_bp_bench — bound propagation benchmark for a single MIP instance
 *
 * Loads a problem, runs CbcBoundPropagation at the requested aggression level,
 * and emits one CSV row with timing and effectiveness statistics.
 *
 * Designed for parallel execution driven by an external script:
 *
 *   # Print header once
 *   mipster_bp_bench --header-only > results.csv
 *
 *   # Solve all instances in parallel (GNU parallel)
 *   ls $MIPSTER_INSTANCES/miplib/2017+spp/*.mps.gz | \
 *     parallel -j$(nproc) 'mipster_bp_bench --no-header {} >> results.csv'
 *
 * Usage:
 *   mipster_bp_bench [OPTIONS] <problem.{mps,lp}[.gz]>
 *
 * Options:
 *   --level <singletons|milpbt|fixpoint>   Aggression level (default: fixpoint)
 *   --max-rounds <N>                        Max BP rounds for milpbt (default: 100)
 *   --no-header                             Suppress CSV header line
 *   --header-only                           Print CSV header and exit (no instance needed)
 *
 * Output columns:
 *   instance         Basename of the problem file (no path, no extension)
 *   ncols            Total number of structural columns
 *   nrows            Total number of constraints
 *   nbin             Number of binary variables
 *   nint             Number of general-integer (non-binary) variables
 *   singleton_fixed  Variables fully fixed by the singleton-row pass
 *   singleton_tight  Variables tightened (not fixed) by the singleton-row pass
 *   bp_fixed         Variables fixed by CoinBoundPropagation rounds
 *   total_fixed      singleton_fixed + bp_fixed
 *   total_tight      singleton_tight (non-fixed tightenings from singleton pass)
 *   rounds           Number of CoinBoundPropagation rounds executed
 *   stop_reason      fixpoint | maxrounds | timelimit | infeasible | not_run
 *   infeasible       1 if infeasibility was proved, 0 otherwise
 *   time_sec         Wall-clock time (seconds) spent in bound propagation
 *   level            Level used: singletons | milpbt | fixpoint
 *
 * Exit codes:
 *   0  success
 *   1  problem file error
 *   2  usage error
 */

#include "Cbc_C_Interface.h"
#include "CbcBoundPropagation.hpp"
#include "CoinTime.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiSolverInterface.hpp"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <string>
#include <unistd.h>

/* ── Silence Cbc/Clp banner output during model loading ─────────────── */

static int devNull_ = -1;
static int savedOut_ = -1;
static int savedErr_ = -1;

static void silenceBegin()
{
  devNull_ = open("/dev/null", O_WRONLY);
  if (devNull_ < 0)
    return;
  savedOut_ = dup(STDOUT_FILENO);
  savedErr_ = dup(STDERR_FILENO);
  dup2(devNull_, STDOUT_FILENO);
  dup2(devNull_, STDERR_FILENO);
}

static void silenceEnd()
{
  fflush(stdout);
  fflush(stderr);
  if (savedOut_ >= 0) {
    dup2(savedOut_, STDOUT_FILENO);
    close(savedOut_);
    savedOut_ = -1;
  }
  if (savedErr_ >= 0) {
    dup2(savedErr_, STDERR_FILENO);
    close(savedErr_);
    savedErr_ = -1;
  }
  if (devNull_ >= 0) {
    close(devNull_);
    devNull_ = -1;
  }
}

/* ── Helpers ─────────────────────────────────────────────────────────── */

static const char *CSV_HEADER =
  "instance,ncols,nrows,nbin,nint,"
  "singleton_fixed,singleton_tight,bp_fixed,"
  "total_fixed,total_tight,rounds,"
  "stop_reason,infeasible,time_sec,level";

static const char *stopReasonStr(CbcBoundPropagation::StopReason r)
{
  switch (r) {
  case CbcBoundPropagation::NotRun:
    return "not_run";
  case CbcBoundPropagation::ReachedFixpoint:
    return "fixpoint";
  case CbcBoundPropagation::HitMaxRounds:
    return "maxrounds";
  case CbcBoundPropagation::HitTimeLimit:
    return "timelimit";
  case CbcBoundPropagation::InfeasibleDetected:
    return "infeasible";
  }
  return "unknown";
}

/* Strip directory and up to two extensions (.mps.gz, .lp.gz, .mps, .lp). */
static std::string instanceBasename(const char *path)
{
  std::string s(path);
  // strip directory
  size_t slash = s.rfind('/');
  if (slash != std::string::npos)
    s = s.substr(slash + 1);
  // strip extensions
  for (int i = 0; i < 2; ++i) {
    size_t dot = s.rfind('.');
    if (dot == std::string::npos)
      break;
    s = s.substr(0, dot);
  }
  return s;
}

static void printUsage(const char *prog)
{
  fprintf(stderr,
    "Usage: %s [OPTIONS] <problem.{mps,lp}[.gz]>\n"
    "\n"
    "Options:\n"
    "  --level <singletons|milpbt|fixpoint>  Aggression level (default: fixpoint)\n"
    "  --max-rounds <N>                       Max rounds for milpbt (default: 100)\n"
    "  --no-header                            Suppress CSV header\n"
    "  --header-only                          Print CSV header and exit\n",
    prog);
}

/* ── main ─────────────────────────────────────────────────────────────── */

int main(int argc, char *argv[])
{
  CbcBoundPropagation::Level level = CbcBoundPropagation::Fixpoint;
  const char *levelStr = "fixpoint";
  int maxRounds = 100;
  bool printHeader = true;
  bool headerOnly = false;
  const char *problemFile = nullptr;

  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--header-only") == 0) {
      headerOnly = true;
    } else if (strcmp(argv[i], "--no-header") == 0) {
      printHeader = false;
    } else if (strcmp(argv[i], "--level") == 0) {
      if (i + 1 >= argc) {
        fprintf(stderr, "Error: --level requires an argument\n");
        return 2;
      }
      ++i;
      if (strcmp(argv[i], "singletons") == 0) {
        level = CbcBoundPropagation::Singletons;
        levelStr = "singletons";
      } else if (strcmp(argv[i], "milpbt") == 0) {
        level = CbcBoundPropagation::MILPbt;
        levelStr = "milpbt";
      } else if (strcmp(argv[i], "fixpoint") == 0) {
        level = CbcBoundPropagation::Fixpoint;
        levelStr = "fixpoint";
      } else {
        fprintf(stderr, "Error: unknown level '%s' (expected: singletons, milpbt, fixpoint)\n", argv[i]);
        return 2;
      }
    } else if (strcmp(argv[i], "--max-rounds") == 0) {
      if (i + 1 >= argc) {
        fprintf(stderr, "Error: --max-rounds requires an argument\n");
        return 2;
      }
      maxRounds = atoi(argv[++i]);
      if (maxRounds <= 0) {
        fprintf(stderr, "Error: --max-rounds must be a positive integer\n");
        return 2;
      }
    } else if (argv[i][0] == '-' && argv[i][1] == '-') {
      fprintf(stderr, "Error: unknown option '%s'\n", argv[i]);
      printUsage(argv[0]);
      return 2;
    } else {
      if (problemFile) {
        fprintf(stderr, "Error: more than one problem file specified\n");
        printUsage(argv[0]);
        return 2;
      }
      problemFile = argv[i];
    }
  }

  if (headerOnly) {
    printf("%s\n", CSV_HEADER);
    return 0;
  }

  if (!problemFile) {
    fprintf(stderr, "Error: no problem file specified\n");
    printUsage(argv[0]);
    return 2;
  }

  /* ── Load problem ─────────────────────────────────────────────────── */

  silenceBegin();
  Cbc_Model *m = Cbc_newModel();
  int readErr = Cbc_readMps(m, problemFile);
  if (readErr) {
    silenceEnd();
    silenceBegin();
    readErr = Cbc_readLp(m, problemFile);
  }
  silenceEnd();

  if (readErr) {
    fprintf(stderr, "Error: could not read problem file '%s'\n", problemFile);
    Cbc_deleteModel(m);
    return 1;
  }

  OsiClpSolverInterface *solver =
    static_cast< OsiClpSolverInterface * >(Cbc_getSolverPtr(m));
  assert(solver);

  const int ncols = solver->getNumCols();
  const int nrows = solver->getNumRows();

  int nbin = 0;
  int nint = 0;
  for (int j = 0; j < ncols; ++j) {
    if (solver->isBinary(j))
      ++nbin;
    else if (solver->isInteger(j))
      ++nint;
  }

  /* ── Run bound propagation ────────────────────────────────────────── */

  CbcBoundPropagation bp;
  const double t0 = CoinGetTimeOfDay();
  bp.run(solver, /*handler=*/nullptr, /*logLevel=*/0,
    level, maxRounds, /*timeLimit=*/1e100, /*startTime=*/t0);
  const double elapsed = CoinGetTimeOfDay() - t0;

  const int singletonFixed = bp.nSingletonFixed();
  const int singletonTight = bp.nSingletonTightened();
  const int bpFixed = bp.nBoundPropFixed();
  const int totalFixed = bp.nFixed();
  const int rounds = bp.nRoundsRun();
  const int infeasible = bp.stopReason() == CbcBoundPropagation::InfeasibleDetected ? 1 : 0;
  const char *stopReason = stopReasonStr(bp.stopReason());

  /* ── Emit CSV row ─────────────────────────────────────────────────── */

  if (printHeader)
    printf("%s\n", CSV_HEADER);

  printf("%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s,%d,%.6f,%s\n",
    instanceBasename(problemFile).c_str(),
    ncols, nrows, nbin, nint,
    singletonFixed, singletonTight,
    bpFixed,
    totalFixed, singletonTight,
    rounds,
    stopReason, infeasible,
    elapsed,
    levelStr);

  Cbc_deleteModel(m);
  return 0;
}
