// Copyright (C) 2026, MIPster contributors.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcHeuristicNoRelRepair_H
#define CbcHeuristicNoRelRepair_H

#include "CbcHeuristic.hpp"

/** No-relaxation repair heuristic.

    This heuristic targets flat-objective feasibility models where LP-guided
    primal heuristics have little signal. It searches directly on row
    violations, using compound moves for exact-one binary rows and coordinate
    repair moves for remaining variables.
*/
class CBCLIB_EXPORT CbcHeuristicNoRelRepair : public CbcHeuristic {
public:
  CbcHeuristicNoRelRepair();
  CbcHeuristicNoRelRepair(CbcModel &model);
  CbcHeuristicNoRelRepair(const CbcHeuristicNoRelRepair &rhs);
  ~CbcHeuristicNoRelRepair();

  CbcHeuristicNoRelRepair &operator=(const CbcHeuristicNoRelRepair &rhs);

  virtual CbcHeuristic *clone() const override;
  virtual void resetModel(CbcModel *model) override;
  virtual void setModel(CbcModel *model) override;
  virtual bool shouldHeurRun(int whereFrom) override;

  using CbcHeuristic::solution;
  virtual int solution(double &objectiveValue, double *newSolution) override;

  inline void setMaxRestarts(int value) { maxRestarts_ = value; }
  inline int maxRestarts() const { return maxRestarts_; }

  inline void setMaxIterations(int value) { maxIterations_ = value; }
  inline int maxIterations() const { return maxIterations_; }

  inline void setMaxSeconds(double value) { maxSeconds_ = value; }
  inline double maxSeconds() const { return maxSeconds_; }

protected:
  int maxRestarts_;
  int maxIterations_;
  double maxSeconds_;
};

#endif
