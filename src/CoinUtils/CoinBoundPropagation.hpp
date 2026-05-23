/* -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; -*-
 *
 * This file is part of the COIN-OR CoinUtils package
 *
 * @file   CoinBoundPropagation.hpp
 * @brief  Bound propagation for binary and general integer variables in a MILP.
 *
 * Copyright (C) 2025
 * All rights reserved.
 *
 * This code is licensed under the terms of the Eclipse Public License (EPL).
 */

#ifndef COIN_BOUND_PROPAGATION_HPP
#define COIN_BOUND_PROPAGATION_HPP

#include <cstddef>
#include <utility>
#include <vector>

#include "CoinUtilsConfig.h"

class CoinPackedMatrix;

#ifdef COIN_BT_STATS
/**
 * @brief Per-row statistics collected when COIN_BT_STATS is defined.
 *
 * Available via CoinBoundPropagation::rowStats().
 * Rows that were skipped (unbounded, sense 'N', zero binary variables)
 * have timeSeconds == 0 and nFixings == 0.
 */
struct CoinBTRowStats {
  size_t rowIdx;         ///< Index of the constraint row.
  double timeSeconds;    ///< Wall-clock time spent processing this row (both knapsack iterations).
  size_t nFixings;       ///< Number of new fixings produced by this row.
  bool   skipped;        ///< True if the row was skipped (unbounded / no binary vars).
};
#endif // COIN_BT_STATS

/**
 * @brief Fast bound tightening for binary and general integer variables.
 *
 * Scans every constraint row of a MILP once and tightens variable bounds:
 *
 * **General integer bound tightening** — for each row in L-form
 * (sum c_i x_i <= b, after applying row sense multiplier):
 *   - minActivity = sum_{c_i>0} c_i*lb_i + sum_{c_i<0} c_i*ub_i
 *   - rowSlack = b - minActivity
 *   - If rowSlack < -tolerance: row is infeasible.
 *   - For general integer j with c_j > 0:  new_ub = lb + floor(rowSlack / c_j)
 *   - For general integer j with c_j < 0:  new_lb = ub - floor(rowSlack / |c_j|)
 *   All sign combinations (positive/negative coefficient, positive/negative
 *   bounds) are handled correctly.  An overflow guard (division-based, avoids
 *   c*range multiplication overflow) skips entries where no tightening is
 *   possible.
 *
 * **Binary fixing** — identifies binary variables that can be fixed
 * (lower bound raised to 1 or upper bound lowered to 0) based on the
 * knapsack structure of each row:
 *
 *   - Each row is converted to the form  sum a_i x_i <= b  (all a_i >= 0)
 *     by discounting fixed and continuous/general-integer variables from the
 *     RHS and complementing variables with negative coefficients.
 *   - A binary variable x_j with coefficient a_j > b + tolerance must be 0
 *     (setting it to 1 alone already violates the constraint).  In complemented
 *     form, x_j appears with a positive coefficient; if that exceeds the RHS
 *     the *original* variable must be 1.
 *
 * General integer tightening runs first; updated bounds immediately feed into
 * the binary fixing pass for the same row, enabling propagation chains
 * (e.g., a tightened general integer tightens the effective RHS, which then
 * forces a binary indicator variable).
 *
 * Fixings discovered early are propagated within the single forward pass:
 * the internal mutable bounds are updated immediately so that subsequent
 * rows benefit from tighter coefficient sums.
 *
 * **Infeasibility detection** — three conditions are checked:
 *   1. Row infeasibility: minActivity > effRHS + tolerance.
 *   2. Direct binary row infeasibility: effective RHS is finite but negative.
 *   3. Contradictory binary fixings: same variable implied to be both 0 and 1.
 *
 * On any condition the scan is aborted immediately.  Call isInfeasible()
 * to check whether this occurred.
 *
 * Note: only the "coefficient exceeds RHS" fixing mechanism is implemented
 * for binaries.  More elaborate probing (fixing a variable and propagating)
 * is out of scope; this class is intentionally lightweight.
 */
class COINUTILSLIB_EXPORT CoinBoundPropagation {
public:
  /**
   * @brief Construct and run the bound-tightening pass.
   *
   * @param numCols   Number of structural columns.
   * @param colType   Column types (see CoinColumnType.hpp).
   * @param colLB     Column lower bounds (length @p numCols).
   * @param colUB     Column upper bounds (length @p numCols).
   * @param matrixByRow Row-wise constraint matrix.
   * @param sense     Row senses ('L', 'G', 'E', 'R', 'N').
   * @param rowRHS    Row right-hand-side values.
   * @param rowRange  Row range values (may be nullptr for non-ranged models).
   * @param primalTolerance Numerical tolerance for tightness checks.
   * @param infinity  Value used as practical infinity when discounting
   *                  continuous variables from the RHS.
   * @param maxRowNz  Hard nonzero limit (>= 0): rows with more nonzeros than
   *                  this are skipped entirely — no infeasibility check, no
   *                  fixings.  Use -1 (the default) for no limit, i.e. the
   *                  complete algorithm.  When any rows are skipped,
   *                  isComplete() returns false.
   *
   * @note The pre-check (which inspects each row without calling processRow
   *       and skips it when no fixing is possible) is always on and always
   *       sound.  For rows where all free binary coefficients are negative in
   *       the multiplier direction (e.g. G-direction rows with only positive
   *       original coefficients), it exits early as soon as rhs ≥ maxCoef ≥ 0,
   *       which is typically after just two elements for dense covering rows.
   */
  CoinBoundPropagation(
    int numCols,
    const char *colType,
    const double *colLB,
    const double *colUB,
    const CoinPackedMatrix *matrixByRow,
    const char *sense,
    const double *rowRHS,
    const double *rowRange,
    double primalTolerance = 1e-7,
    double infinity = 1e50,
    int maxRowNz = -1);

  ~CoinBoundPropagation() = default;

  CoinBoundPropagation(const CoinBoundPropagation &) = delete;
  CoinBoundPropagation &operator=(const CoinBoundPropagation &) = delete;
  CoinBoundPropagation(CoinBoundPropagation &&) = delete;
  CoinBoundPropagation &operator=(CoinBoundPropagation &&) = delete;

  /**
   * @brief Tightened bounds discovered during construction.
   *
   * Each entry has the form (column_index, (new_lb, new_ub)).
   * For binary variables the new bound is either (0,0) or (1,1).
   *
   * If isInfeasible() returns true the contents of this vector reflect
   * only the fixings recorded *before* the infeasibility was detected.
   *
   * If isComplete() returns false some rows were skipped (due to maxRowNz),
   * so this list may be missing fixings or infeasibilities from those rows.
   */
  const std::vector< std::pair< size_t, std::pair< double, double > > > &updatedBounds() const
  {
    return newBounds_;
  }

  /** @brief Number of variable fixings discovered. */
  size_t nFixings() const { return newBounds_.size(); }

  /**
   * @brief True if the bound-tightening pass detected infeasibility.
   *
   * Infeasibility is reported when either:
   *  - a processed knapsack row has a finite but negative effective RHS
   *    (the constraint is impossible to satisfy), or
   *  - the same variable is implied to be both 0 and 1 by two different rows.
   *
   * Always false when isComplete() is false and the actual infeasibility
   * resides in a skipped row.
   */
  bool isInfeasible() const { return infeasible_; }

  /**
   * @brief Row index that triggered infeasibility, or -1 if unknown.
   *
   * Valid only when isInfeasible() is true.  Set when the effective RHS of
   * a knapsack row becomes strictly negative (the constraint cannot be
   * satisfied by any binary assignment), OR when a contradictory fixing is
   * detected (row that first observed the conflict).
   */
  int infeasibleRow() const { return infeasibleRow_; }

  /**
   * @brief Column index involved in a contradictory fixing, or -1 if unknown.
   *
   * Valid only when isInfeasible() is true.  Set when the same binary
   * variable is implied to be both 0 and 1 by two different rows.
   * Returns -1 for row-RHS infeasibility (see infeasibleRow()).
   */
  int infeasibleCol() const { return infeasibleCol_; }

  /**
   * @brief True if every row was processed (or skipped by the sound
   *        pre-check only).
   *
   * Returns false when the @p maxRowNz constructor parameter caused at
   * least one row to be hard-skipped.  In that case nFixings() and
   * isInfeasible() may be incomplete.
   */
  bool isComplete() const { return complete_; }

#ifdef COIN_BT_STATS
  /**
   * @brief Per-row statistics (only available when compiled with COIN_BT_STATS).
   *
   * One entry per row of the constraint matrix, in row order.
   * Use this to identify rows that consume significant time but contribute
   * no fixings (skipped == false && nFixings == 0 && timeSeconds is large).
   */
  const std::vector< CoinBTRowStats > &rowStats() const { return rowStats_; }
#endif // COIN_BT_STATS

private:
  std::vector< std::pair< size_t, std::pair< double, double > > > newBounds_;
  bool infeasible_;
  int infeasibleRow_; ///< row that triggered infeasibility (-1 if none / unknown)
  int infeasibleCol_; ///< column in contradictory fixing (-1 if none / unknown)
  bool complete_;     ///< false if any row was hard-skipped via maxRowNz
#ifdef COIN_BT_STATS
  std::vector< CoinBTRowStats > rowStats_;
#endif
};

#endif // COIN_BOUND_PROPAGATION_HPP
