// Copyright (C) 2024 COIN-OR Foundation
// Authors: Cbc development team
// This code is licensed under the terms of the Eclipse Public License (EPL)

#include "CbcBoundPropagation.hpp"

#include "CoinBoundPropagation.hpp"
#include "CoinMessageHandler.hpp"
#include "CoinTime.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiSolverInterface.hpp"

// Reference solution loaded by -debugCuts; set before applyLpMethod() so that
// bound propagation can check fixings before the OsiRowCutDebugger is active.
extern double *debugSolution;
extern int debugNumberColumns;

#include <cassert>
#include <climits>
#include <cmath>
#include <string>
#include <vector>

CbcBoundPropagation::CbcBoundPropagation()
  : nSingletonTightened_(0)
  , nSingletonFixed_(0)
  , nBoundPropFixed_(0)
  , nFBBTTightened_(0)
  , nRoundsRun_(0)
  , stopReason_(NotRun)
  , timeUsed_(0.0)
  , infeasibleRow_(-1)
  , infeasibleCol_(-1)
  , nonBinaryFBBT_(true)
{
}

bool CbcBoundPropagation::run(OsiSolverInterface *solver,
  CoinMessageHandler * /*handler*/,
  int logLevel,
  Level level,
  int maxRounds,
  double timeLimit,
  double startTime)
{
  assert(level != Off);

  if (!solver || solver->getNumCols() == 0)
    return true;

  const double t0 = CoinGetTimeOfDay();

  // If a debugCuts reference solution is loaded, check every bound fixing we
  // apply against that solution.
  //
  // Before the LP solve the OsiRowCutDebugger is not yet active, so we fall
  // back to the debugSolution global which is populated from the -debugCuts
  // sol file at the very start of applyLpMethod().
  //
  // During B&B we use getRowCutDebugger (without "Always"): it returns NULL
  // at subtree nodes whose branching has already excluded the reference
  // solution — fixing variables differently from the global optimal is correct
  // there, so we must NOT flag those.  Only on the optimal path is a
  // contradictory fixing a bug.
  //
  // Fallback to debugSolution only when no OsiRowCutDebugger has been
  // activated at all (pre-B&B bound propagation). Once the debugger is
  // active, debugger==NULL means we are on a wrong subtree — never check.
  const OsiRowCutDebugger *debugger = solver->getRowCutDebugger();
  const OsiRowCutDebugger *debuggerAlways = solver->getRowCutDebuggerAlways();
  const double *optSol = debugger
    ? debugger->optimalSolution()
    : (!debuggerAlways && debugSolution && debugNumberColumns == solver->getNumCols()
         ? debugSolution
         : nullptr);

  // Declare these early so they are in scope for the checkFixing lambda.
  // colType and curLB/curUB are set to their real values before phase 2.
  const int nCols = solver->getNumCols();
  // Local owned copy — not a pointer into the solver's internal cache.
  // Refreshed before each propagation round from current bounds.
  std::vector< char > colTypeBuf;
  const char *colType = nullptr;       // points to colTypeBuf.data() after phase-2 init
  std::vector< double > curLB, curUB;  // initialised before phase-2 loop

  // diagRound == -1: singleton phase;  >= 0: CoinBoundPropagation round.
  int diagRound = -1;

  // Check one fixing against the reference solution (optimal path only).
  // Called with the new bounds; curLB/curUB must still hold OLD bounds when
  // called from the propagation loop (call before updating curLB/curUB).
  auto checkFixing = [&](int col, double newLB, double newUB,
                         const char *phase) {
    if (!optSol || !solver->isInteger(col))
      return;
    const double sv = optSol[col];
    if (newLB > sv + 0.5 || newUB < sv - 0.5) {
      const std::string name = solver->getColName(col);
      const char ct = colType ? colType[col] : char(-1);
      const char *ctName = (ct == 1) ? "binary"
                         : (ct == 2) ? "general-integer"
                         : (ct == 0) ? "continuous" : "?";
      // curLB/curUB contain old bounds when called from propagation loop.
      const double oldLB = (col < static_cast< int >(curLB.size()))
                             ? curLB[col] : newLB;
      const double oldUB = (col < static_cast< int >(curUB.size()))
                             ? curUB[col] : newUB;
      printf("nodeBoundProp BAD FIXING (%s, round %d): col %d (%s)"
             " type=%s old=[%g,%g] new=[%g,%g] but optimal has %g\n",
             phase, diagRound, col, name.c_str(), ctName,
             oldLB, oldUB, newLB, newUB, sv);
      fflush(stdout);
    }
  };

  // ---------------------------------------------------------------
  // Phase 1: singleton row tightening
  // ---------------------------------------------------------------
  {
    int nFixed = 0;
    const int nTightened = solver->tightenBoundsFromSingletonRows(nFixed);

    if (nTightened < 0) {
      // infeasibility detected — singleton API does not expose row/col source
      stopReason_ = InfeasibleDetected;
      timeUsed_ = (CoinGetTimeOfDay()) - t0;

      if (logLevel >= 1)
        printf("  Bound propagation: INFEASIBLE (singleton tightening), "
               "%.3f s.\n",
          timeUsed_);

      return false;
    }

    nSingletonFixed_ = nFixed;
    nSingletonTightened_ = nTightened - nFixed;

    // Check singleton fixings against the reference solution (optimal path only).
    // curLB/curUB are not yet populated; oldLB/oldUB will show newLB/newUB.
    if (optSol) {
      // Capture stable copies — the raw pointers from getColLower/Upper() may
      // be invalidated by intermediate solver calls (e.g. getColName(),
      // isInteger()) that trigger internal array reallocation.
      const double *rawLB = solver->getColLower();
      const double *rawUB = solver->getColUpper();
      const std::vector<double> lbVec(rawLB, rawLB + nCols);
      const std::vector<double> ubVec(rawUB, rawUB + nCols);
      for (int j = 0; j < nCols; j++)
        checkFixing(j, lbVec[j], ubVec[j], "singleton");
    }

    if (logLevel >= 2 && nTightened > 0)
      printf("  Bound propagation (singletons): %d tightened, %d fixed.\n",
        nSingletonTightened_, nSingletonFixed_);
  }

  if (level == Singletons) {
    stopReason_ = ReachedFixpoint;
    timeUsed_ = (CoinGetTimeOfDay()) - t0;

    if (logLevel >= 1)
      printf("  Bound propagation fixed %d vars in %.3f s.\n",
        nSingletonFixed_, timeUsed_);

    return true;
  }

  // ---------------------------------------------------------------
  // Phase 2: CoinBoundPropagation — iterate until fixpoint or limits
  // ---------------------------------------------------------------
  const int roundLimit = (level == Fixpoint) ? INT_MAX : maxRounds;

  // Initialise phase-2 data (after singleton tightening so bounds are current).
  // Copy colType into a local buffer so we own the data (not a pointer into
  // the solver's internal cache which may be stale from a prior invocation).
  // Use refresh=true to force recomputation from current bounds — a stale
  // cached value can misclassify GeneralInteger variables as Binary, causing
  // wrong fixings (root cause of the miclsp arm64 soundness bug).
  {
    const char *ct = solver->getColType(true);
    colTypeBuf.assign(ct, ct + nCols);
  }
  colType = colTypeBuf.data();
  const CoinPackedMatrix *matByRow = solver->getMatrixByRow();
  const char *rowSense = solver->getRowSense();
  const double *rhs = solver->getRightHandSide();
  const double *range = solver->getRowRange();
  double primalTol = 1e-7;
  solver->getDblParam(OsiPrimalTolerance, primalTol);
  const double infinity = solver->getInfinity();

  // Working copies of bounds (updated after each round).
  curLB.assign(solver->getColLower(), solver->getColLower() + nCols);
  curUB.assign(solver->getColUpper(), solver->getColUpper() + nCols);

  // Refresh colTypeBuf from current bounds before each round.
  // Cheaper than calling solver->getColType() and avoids stale-cache issues.
  // A GeneralInteger whose bounds have been tightened to [0,1] is reclassified
  // as Binary, enabling stronger propagation in later rounds.
  auto refreshColType = [&]() {
    for (int j = 0; j < nCols; ++j) {
      if (colTypeBuf[j] == 0)
        continue; // continuous: never changes
      colTypeBuf[j] = (curLB[j] >= 0.0 && curUB[j] <= 1.0) ? 1 : 2;
    }
  };

  // ── Dirty-row tracking for FBBT (lazy) ───────────────────────────────────
  // Round 1 always runs without a dirty hint (CoinBoundPropagation builds its
  // own rowHasNonBinary).  The O(nnz) colToRows adjacency list and the dirty
  // bitvector are built lazily after round 1 — only if there will be a round 2.
  // This avoids paying the build cost for problems that converge in one round.
  const int nRows = matByRow->getNumRows();
  const int *matIdxs = matByRow->getIndices();
  const CoinBigIndex *matStart = matByRow->getVectorStarts();
  const int *matLen = matByRow->getVectorLengths();

  bool hasNonBinaryVars = false;
  if (nonBinaryFBBT_) {
    for (int j = 0; j < nCols && !hasNonBinaryVars; ++j)
      if (colTypeBuf[j] != 1) // 1 = Binary
        hasNonBinaryVars = true;
  }

  // Built lazily after round 1 (if a second round is needed).
  bool dirtyInfraBuilt = false;
  std::vector< CoinBigIndex > colRowStart;
  std::vector< int > colRowList;
  std::vector< double > colRowCoef; // coefficient of col j in each (col,row) entry
  std::vector< bool > rowHasNonBinaryBP;
  std::vector< char > rowHasBinaryBP; // char to allow .data() (vector<bool> lacks it)
  std::vector< bool > dirtyRowsFBBT;

  // Cached row min activity (for ≤ constraint FBBT, incremental updates).
  // rowCachedMinAct[r] = sum(a>0: a*lb_j, a<0: a*ub_j) over all vars in row r.
  // rowCachedNUnbLB[r] = number of vars with a*lb_j = -inf (contributes -inf to minAct).
  // Valid when actCacheBuilt == true; updated incrementally from bound changes.
  bool actCacheBuilt = false;
  std::vector< double > rowCachedMinAct;
  std::vector< int > rowCachedNUnbLB;

  // Per-round buffer recording (col, oldLB, oldUB, newLB, newUB) for each committed
  // bound change. Used in the post-round single adjacency pass that simultaneously
  // updates the min-activity cache and marks dirty rows. Pre-allocated once, reused
  // every round to avoid repeated heap allocations.
  struct BoundChange {
    int col;
    double oldLB, oldUB, newLB, newUB;
  };
  std::vector< BoundChange > changedBounds;

  // Build colToRows (column → rows adjacency) and rowHasNonBinaryBP once.
  // Called at most once, the first time a round 2 is needed.
  auto buildDirtyInfra = [&]() {
    if (dirtyInfraBuilt)
      return;
    dirtyInfraBuilt = true;

    rowHasNonBinaryBP.assign(static_cast< size_t >(nRows), false);
    rowHasBinaryBP.assign(static_cast< size_t >(nRows), char(0));
    for (int r = 0; r < nRows; ++r) {
      const CoinBigIndex rs = matStart[r];
      const int len = matLen[r];
      for (int k = 0; k < len; ++k) {
        const int j = matIdxs[rs + k];
        if (colTypeBuf[j] != 1)
          rowHasNonBinaryBP[r] = true;
        else
          rowHasBinaryBP[r] = char(1);
      }
    }

    // colToRows: for each column j, list the rows and coefficients.
    colRowStart.assign(static_cast< size_t >(nCols + 1), CoinBigIndex(0));
    for (int r = 0; r < nRows; ++r) {
      const CoinBigIndex rs = matStart[r];
      const int len = matLen[r];
      for (int k = 0; k < len; ++k)
        ++colRowStart[matIdxs[rs + k] + 1];
    }
    for (int j = 0; j < nCols; ++j)
      colRowStart[j + 1] += colRowStart[j];
    colRowList.resize(static_cast< size_t >(colRowStart[nCols]));
    colRowCoef.resize(static_cast< size_t >(colRowStart[nCols]));
    {
      std::vector< CoinBigIndex > pos(colRowStart.begin(),
        colRowStart.begin() + nCols);
      const double *matCoefs = matByRow->getElements();
      for (int r = 0; r < nRows; ++r) {
        const CoinBigIndex rs = matStart[r];
        const int len = matLen[r];
        for (int k = 0; k < len; ++k) {
          const int j = matIdxs[rs + k];
          const CoinBigIndex p = pos[j]++;
          colRowList[p] = r;
          colRowCoef[p] = matCoefs[rs + k];
        }
      }
    }

    // Initialise the row min-activity cache from the current bounds (curLB/curUB).
    rowCachedMinAct.assign(static_cast< size_t >(nRows), 0.0);
    rowCachedNUnbLB.assign(static_cast< size_t >(nRows), 0);
    {
      const double *matCoefs = matByRow->getElements();
      for (int r = 0; r < nRows; ++r) {
        const CoinBigIndex rs = matStart[r];
        const int len = matLen[r];
        double minA = 0.0;
        int nUnbLB = 0;
        for (int k = 0; k < len; ++k) {
          const int j = matIdxs[rs + k];
          const double a = matCoefs[rs + k];
          const double lb = curLB[j], ub = curUB[j];
          if (a > 0.0) {
            if (lb <= -infinity)
              ++nUnbLB;
            else
              minA += a * lb;
          } else if (a < 0.0) {
            if (ub >= infinity)
              ++nUnbLB;
            else
              minA += a * ub;
          }
        }
        rowCachedMinAct[r] = minA;
        rowCachedNUnbLB[r] = nUnbLB;
      }
    }
    actCacheBuilt = true;

    dirtyRowsFBBT.assign(static_cast< size_t >(nRows), false);
  };

  for (int round = 0; round < roundLimit; ++round) {
    diagRound = round;
    refreshColType(); // update binary/general-integer from current curLB/curUB
    // Time-limit check before starting this round
    const double tNow = CoinGetTimeOfDay();
    if (tNow - startTime >= timeLimit) {
      stopReason_ = HitTimeLimit;
      timeUsed_ = tNow - t0;

      if (logLevel >= 1)
        printf("  Bound propagation: fixed %d (%d singleton + %d propagation, "
               "%d round(s), TIME LIMIT), %.3f s.\n",
          nSingletonFixed_ + nBoundPropFixed_, nSingletonFixed_, nBoundPropFixed_,
          nRoundsRun_, timeUsed_);

      return true;
    }

    // Round 0 and 1 always use nullptr (no dirty hint; CoinBoundPropagation
    // builds rowHasNonBinary internally).  Round 2+ use the dirty set built
    // after round 1 completes.
    const std::vector< bool > *dirtyHint =
      (hasNonBinaryVars && dirtyInfraBuilt) ? &dirtyRowsFBBT : nullptr;
    const double *cachedMinAct = actCacheBuilt ? rowCachedMinAct.data() : nullptr;
    const int *cachedNUnbLB = actCacheBuilt ? rowCachedNUnbLB.data() : nullptr;
    const bool *hasBinaryRow =
      dirtyInfraBuilt ? reinterpret_cast< const bool * >(rowHasBinaryBP.data()) : nullptr;
    CoinBoundPropagation bt(nCols, colType,
      curLB.data(), curUB.data(),
      matByRow, rowSense, rhs, range,
      primalTol, infinity,
      /*maxRowNz=*/-1, /*collectCases=*/false,
      nonBinaryFBBT_, dirtyHint,
      cachedMinAct, nullptr, cachedNUnbLB, nullptr,
      hasBinaryRow);

    ++nRoundsRun_;

    if (bt.isInfeasible()) {
      infeasibleRow_ = bt.infeasibleRow();
      infeasibleCol_ = bt.infeasibleCol();
      stopReason_ = InfeasibleDetected;
      timeUsed_ = (CoinGetTimeOfDay()) - t0;

      if (logLevel >= 1) {
        if (infeasibleRow_ >= 0 && infeasibleCol_ >= 0) {
          const std::string rowName = (infeasibleRow_ < solver->getNumRows())
            ? solver->getRowName(infeasibleRow_)
            : "(unknown)";
          const std::string colName = (infeasibleCol_ < solver->getNumCols())
            ? solver->getColName(infeasibleCol_)
            : "(unknown)";
          printf("  Bound propagation: INFEASIBLE in round %d — "
                 "row %d (%s), col %d (%s), %.3f s.\n",
            nRoundsRun_, infeasibleRow_, rowName.c_str(),
            infeasibleCol_, colName.c_str(), timeUsed_);
        } else if (infeasibleRow_ >= 0) {
          const std::string rowName = (infeasibleRow_ < solver->getNumRows())
            ? solver->getRowName(infeasibleRow_)
            : "(unknown)";
          printf("  Bound propagation: INFEASIBLE in round %d — "
                 "row %d (%s), %.3f s.\n",
            nRoundsRun_, infeasibleRow_, rowName.c_str(), timeUsed_);
        } else {
          printf("  Bound propagation: INFEASIBLE in round %d, %.3f s.\n",
            nRoundsRun_, timeUsed_);
        }
      }

      return false;
    }

    // Apply the fixings from this round to curLB/curUB and to the solver.
    // newBounds_ contains binary fixings (first) then FBBT tightenings (appended
    // at end of constructor), but we apply all in one loop.
    const auto &bounds = bt.updatedBounds();
    const int nNew = static_cast< int >(bounds.size());
    const int nFBBT = bt.nContinuousTightened();
    const int nFixed = nNew - nFBBT;

    if (nNew == 0) {
      stopReason_ = ReachedFixpoint;
      timeUsed_ = (CoinGetTimeOfDay()) - t0;

      if (logLevel >= 2)
        printf("  Bound propagation: fixpoint reached after %d "
               "round(s).\n",
          nRoundsRun_);

      break;
    }

    // Commit bound changes and record old/new for the post-round adjacency pass.
    changedBounds.clear();
    for (const auto &p : bounds) {
      const int col = static_cast< int >(p.first);
      const double newLB = p.second.first;
      const double newUB = p.second.second;
      const double oldLB = curLB[col];
      const double oldUB = curUB[col];
      // Check BEFORE updating curLB/curUB so the lambda sees the old bounds.
      checkFixing(col, newLB, newUB, "propagation");
      curLB[col] = newLB;
      curUB[col] = newUB;
      solver->setColLower(col, newLB);
      solver->setColUpper(col, newUB);
      if (hasNonBinaryVars)
        changedBounds.push_back({ col, oldLB, oldUB, newLB, newUB });
    }

    nBoundPropFixed_ += nFixed;
    nFBBTTightened_ += nFBBT;

    // Post-round single adjacency pass: update min-activity cache AND mark dirty rows.
    // Deferred until round >= 1 so fast 1-2 round problems never pay the build cost.
    if (hasNonBinaryVars && round >= 1) {
      const bool cacheAlreadyBuilt = actCacheBuilt; // false on first call (round 1)
      buildDirtyInfra(); // builds cache from current curLB/curUB (first time only)
      std::fill(dirtyRowsFBBT.begin(), dirtyRowsFBBT.end(), false);
      for (const BoundChange &bc : changedBounds) {
        const double dLB = bc.newLB - bc.oldLB;
        const double dUB = bc.newUB - bc.oldUB;
        for (CoinBigIndex ri = colRowStart[bc.col];
             ri < colRowStart[bc.col + 1]; ++ri) {
          const int r = colRowList[ri];
          // Update min-activity cache — only when the cache was already built
          // BEFORE this round's changes (i.e., not on the initial build round).
          // On the initial build, curLB/curUB already include round 1's changes,
          // so applying the deltas again would double-count.
          if (cacheAlreadyBuilt) {
            const double a = colRowCoef[ri];
            if (a > 0.0 && dLB != 0.0) {
              if (bc.oldLB <= -infinity) {
                --rowCachedNUnbLB[r];
                rowCachedMinAct[r] += a * bc.newLB;
              } else {
                rowCachedMinAct[r] += a * dLB;
              }
            } else if (a < 0.0 && dUB != 0.0) {
              if (bc.oldUB >= infinity) {
                --rowCachedNUnbLB[r];
                rowCachedMinAct[r] += a * bc.newUB;
              } else {
                rowCachedMinAct[r] += a * dUB;
              }
            }
          }
          // Mark dirty rows for next round.
          if (rowHasNonBinaryBP[r])
            dirtyRowsFBBT[r] = true;
        }
      }
    }

    if (logLevel >= 2)
      printf("  Bound propagation: round %d fixed %d vars, FBBT tightened %d"
             " (total fixed %d).\n",
        nRoundsRun_, nFixed, nFBBT, nBoundPropFixed_);
  }

  if (stopReason_ == NotRun) {
    // Exited by roundLimit without fixpoint
    stopReason_ = HitMaxRounds;
  }

  timeUsed_ = (CoinGetTimeOfDay()) - t0;

  if (logLevel >= 1) {
    const int totalFixed = nSingletonFixed_ + nBoundPropFixed_;
    printf("  Bound propagation fixed %d vars, FBBT tightened %d in %.3f s.\n",
      totalFixed, nFBBTTightened_, timeUsed_);
  }

  return true;
}
