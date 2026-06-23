/* -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; -*-
 *
 * This file is part of the COIN-OR CoinUtils package
 *
 * @file   CoinFBBT.hpp
 * @brief  Feasibility-Based Bounds Tightening (FBBT) for all variable types.
 *
 * Copyright (C) 2025
 * All rights reserved.
 *
 * This code is licensed under the terms of the Eclipse Public License (EPL).
 */

#ifndef COIN_FBBT_HPP
#define COIN_FBBT_HPP

#include <cstdint>
#include <utility>
#include <vector>

#include "CoinUtilsConfig.h"

class CoinPackedMatrix;

/**
 * @brief Feasibility-Based Bounds Tightening (FBBT) for a MILP.
 *
 * For each constraint row `Σ aᵢxᵢ ≤ b` (and its complement for G/E/R rows),
 * computes the minimum activity of the row and uses it to derive tighter
 * bounds for every variable in the row.  All variable types are handled:
 * binary, general integer, and continuous.
 *
 * ### Free and unbounded variables
 *
 * A variable is an *unbounded contributor* for a given row direction when its
 * bound in the constraining direction is infinite:
 *   - coefficient > 0  and  lb = -∞, or
 *   - coefficient < 0  and  ub = +∞.
 *
 * Three cases:
 *   - **0 unbounded contributors** — minimum activity is finite; standard
 *     FBBT applies; infeasibility is detectable when min_act > rhs + tol.
 *   - **1 unbounded contributor** (variable k) — minimum activity is -∞, so
 *     no tightening of other variables is possible.  However, variable k
 *     itself can be given a finite bound from the finite contributions of all
 *     other variables.  This is particularly valuable for free variables.
 *   - **≥ 2 unbounded contributors** — nothing can be tightened; row skipped.
 *
 * ### Within-pass propagation
 *
 * Tightenings are applied to internal working copies of the bounds and flags
 * as they are discovered, so later rows in the same pass benefit from earlier
 * tightenings.  If a previously unbounded variable receives a finite bound,
 * its COL_LB_INF / COL_UB_INF flag is cleared immediately, enabling it to
 * contribute to min-activity in subsequent rows.
 *
 * ### Integer rounding
 *
 * For integer variables, the raw tightened bound is rounded conservatively:
 *   - upper bound: `floor(raw_ub + primalTol)`
 *   - lower bound: `ceil(raw_lb − primalTol)`
 *
 * ### Usage
 * @code
 *   // Build flags once (refresh if bounds change significantly)
 *   std::vector<uint8_t> flags(ncols);
 *   CoinFBBT::buildColFlags(ncols, colType, lb, ub, infinity, flags.data());
 *
 *   // Run one FBBT pass
 *   CoinFBBT fbbt(ncols, flags.data(), lb, ub,
 *                 matByRow, sense, rhs, range, primalTol, infinity);
 *
 *   if (fbbt.isInfeasible()) { ... }
 *   for (auto& [col, bounds] : fbbt.updatedBounds()) {
 *     solver->setColLower(col, bounds.first);
 *     solver->setColUpper(col, bounds.second);
 *   }
 * @endcode
 */
class COINUTILSLIB_EXPORT CoinFBBT {
public:
  // ── Column flags pre-computed by buildColFlags() ──────────────────────────

  static constexpr uint8_t COL_FIXED   = 0x01; ///< lb == ub: discount from RHS, never tighten
  static constexpr uint8_t COL_SEMI    = 0x02; ///< semi-continuous or semi-integer: skip row
  static constexpr uint8_t COL_LB_INF  = 0x04; ///< lb ≤ −infinity: unbounded contributor when c > 0
  static constexpr uint8_t COL_UB_INF  = 0x08; ///< ub ≥ +infinity: unbounded contributor when c < 0
  static constexpr uint8_t COL_INTEGER = 0x10; ///< integer (binary or general): apply floor/ceil

  /**
   * @brief Pre-compute per-column flags from column type and bound arrays.
   *
   * Call this once before constructing CoinFBBT.  Refresh whenever bounds
   * change (e.g. after each binary-propagation round).
   *
   * @param numCols  Number of structural columns.
   * @param colType  Array of CoinColumnType::Code values (length numCols).
   * @param colLB    Column lower bounds (length numCols).
   * @param colUB    Column upper bounds (length numCols).
   * @param infinity Value used as practical infinity (e.g. 1e50).
   * @param colFlags Output array (length numCols); caller allocates.
   */
  static void buildColFlags(
    int numCols,
    const char *colType,
    const double *colLB,
    const double *colUB,
    double infinity,
    uint8_t *colFlags);

  /**
   * @brief Run one FBBT pass over all rows.
   *
   * Tightenings are collected internally; the input arrays colLB/colUB are
   * NOT modified.  Use updatedBounds() to retrieve results and apply them to
   * the solver.
   *
   * @param numCols        Number of structural columns.
   * @param colFlags       Pre-computed column flags (from buildColFlags()).
   * @param colLB          Column lower bounds.
   * @param colUB          Column upper bounds.
   * @param matByRow       Row-wise constraint matrix.
   * @param rowSense       Row senses ('L', 'G', 'E', 'R', 'N').
   * @param rowRHS         Row right-hand-side values.
   * @param rowRange       Row range values (may be nullptr for non-ranged models).
   * @param primalTolerance Numerical tolerance.
   * @param infinity       Practical infinity value.
   */
  CoinFBBT(
    int numCols,
    const uint8_t *colFlags,
    const double *colLB,
    const double *colUB,
    const CoinPackedMatrix *matByRow,
    const char *rowSense,
    const double *rowRHS,
    const double *rowRange,
    double primalTolerance = 1e-7,
    double infinity = 1e50);

  ~CoinFBBT() = default;
  CoinFBBT(const CoinFBBT &) = delete;
  CoinFBBT &operator=(const CoinFBBT &) = delete;
  CoinFBBT(CoinFBBT &&) = delete;
  CoinFBBT &operator=(CoinFBBT &&) = delete;

  /**
   * @brief Tightened bounds found during the pass.
   *
   * Each entry is (column_index, (new_lb, new_ub)).  Both bounds are always
   * present; a variable that only got its UB tightened still carries the
   * original LB in the pair, and vice versa.
   */
  const std::vector<std::pair<int, std::pair<double, double>>> &updatedBounds() const
  {
    return newBounds_;
  }

  /// Number of variables with at least one bound tightened.
  size_t nTightenings() const { return newBounds_.size(); }

  /// True if infeasibility was proved during the pass.
  bool isInfeasible() const { return infeasible_; }

  /// Row that triggered infeasibility, or -1 if not detected.
  int infeasibleRow() const { return infeasibleRow_; }

private:
  std::vector<std::pair<int, std::pair<double, double>>> newBounds_;
  bool infeasible_;
  int infeasibleRow_;
};

#endif // COIN_FBBT_HPP
