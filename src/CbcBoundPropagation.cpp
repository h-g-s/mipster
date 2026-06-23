// Copyright (C) 2024 COIN-OR Foundation
// Authors: Cbc development team
// This code is licensed under the terms of the Eclipse Public License (EPL)

#include "CbcBoundPropagation.hpp"

#include "CoinBoundPropagation.hpp"
#include "CoinFBBT.hpp"
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
{
}

bool CbcBoundPropagation::run(OsiSolverInterface *solver,
  CoinMessageHandler * /*handler*/,
  int logLevel,
  Level level,
  int maxRounds,
  double timeLimit,
  double startTime,
  bool enableFBBT)
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

    // Construction runs the algorithm; results are immediately available.
    CoinBoundPropagation bt(nCols, colType,
      curLB.data(), curUB.data(),
      matByRow, rowSense, rhs, range,
      primalTol, infinity);

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

    // Apply the fixings from this round to curLB/curUB and to the solver
    const auto &bounds = bt.updatedBounds();
    const int nNew = static_cast< int >(bounds.size());

    if (nNew == 0) {
      stopReason_ = ReachedFixpoint;
      timeUsed_ = (CoinGetTimeOfDay()) - t0;

      if (logLevel >= 2)
        printf("  Bound propagation: fixpoint reached after %d "
               "round(s).\n",
          nRoundsRun_);

      break;
    }

    for (const auto &p : bounds) {
      const int col = static_cast< int >(p.first);
      // Check BEFORE updating curLB/curUB so the lambda sees the old bounds.
      checkFixing(col, p.second.first, p.second.second, "propagation");
      curLB[col] = p.second.first;
      curUB[col] = p.second.second;
      solver->setColLower(col, p.second.first);
      solver->setColUpper(col, p.second.second);
    }

    nBoundPropFixed_ += nNew;

    if (logLevel >= 2)
      printf("  Bound propagation: round %d fixed %d variables "
             "(total %d).\n",
        nRoundsRun_, nNew, nBoundPropFixed_);
  }

  if (stopReason_ == NotRun) {
    // Exited by roundLimit without fixpoint
    stopReason_ = HitMaxRounds;
  }

  // ---------------------------------------------------------------
  // Phase 3 (optional): FBBT — tighten general integer and continuous bounds
  // ---------------------------------------------------------------
  if (enableFBBT && stopReason_ != InfeasibleDetected) {
    // Refresh colType from current bounds (binaries may have been fixed).
    refreshColType();

    // Pre-compute per-column flags.
    std::vector< uint8_t > colFlagsVec(static_cast< size_t >(nCols));
    CoinFBBT::buildColFlags(
      nCols, colType,
      curLB.data(), curUB.data(),
      infinity, colFlagsVec.data());

    CoinFBBT fbbt(
      nCols, colFlagsVec.data(),
      curLB.data(), curUB.data(),
      matByRow, rowSense, rhs, range,
      primalTol, infinity);

    if (fbbt.isInfeasible()) {
      infeasibleRow_ = fbbt.infeasibleRow();
      infeasibleCol_ = -1;
      stopReason_ = InfeasibleDetected;
      timeUsed_ = (CoinGetTimeOfDay()) - t0;

      if (logLevel >= 1)
        printf("  Bound propagation (FBBT): INFEASIBLE in row %d, %.3f s.\n",
          infeasibleRow_, timeUsed_);

      return false;
    }

    const auto &fbbtBounds = fbbt.updatedBounds();
    nFBBTTightened_ = static_cast< int >(fbbtBounds.size());

    for (const auto &p : fbbtBounds) {
      const int col = p.first;
      checkFixing(col, p.second.first, p.second.second, "fbbt");
      curLB[col] = p.second.first;
      curUB[col] = p.second.second;
      solver->setColLower(col, p.second.first);
      solver->setColUpper(col, p.second.second);
    }

    if (logLevel >= 2 && nFBBTTightened_ > 0)
      printf("  Bound propagation (FBBT): tightened %d variable bounds.\n",
        nFBBTTightened_);
  }

  timeUsed_ = (CoinGetTimeOfDay()) - t0;

  if (logLevel >= 1) {
    const int totalFixed = nSingletonFixed_ + nBoundPropFixed_;
    printf("  Bound propagation fixed %d vars in %.3f s.\n",
      totalFixed, timeUsed_);
  }

  return true;
}
