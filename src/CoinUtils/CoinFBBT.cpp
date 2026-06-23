/* -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; -*-
 *
 * This file is part of the COIN-OR CoinUtils package
 *
 * @file   CoinFBBT.cpp
 * @brief  Feasibility-Based Bounds Tightening (FBBT) for all variable types.
 *
 * Copyright (C) 2025
 * All rights reserved.
 *
 * This code is licensed under the terms of the Eclipse Public License (EPL).
 */

#include "CoinFBBT.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "CoinColumnType.hpp"
#include "CoinKnapsackRow.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPragma.hpp"

void CoinFBBT::buildColFlags(
  int numCols,
  const char *colType,
  const double *colLB,
  const double *colUB,
  double infinity,
  uint8_t *colFlags)
{
  for (int j = 0; j < numCols; ++j) {
    uint8_t f = 0;
    const double lb = colLB[j];
    const double ub = colUB[j];
    const char ct = colType[j];

    if (lb == ub)
      f |= COL_FIXED;

    if (ct == CoinColumnType::SemiContinuous || ct == CoinColumnType::SemiInteger)
      f |= COL_SEMI;

    if (lb <= -infinity)
      f |= COL_LB_INF;

    if (ub >= infinity)
      f |= COL_UB_INF;

    if (ct == CoinColumnType::Binary
      || ct == CoinColumnType::GeneralInteger
      || ct == CoinColumnType::SemiInteger)
      f |= COL_INTEGER;

    colFlags[j] = f;
  }
}

CoinFBBT::CoinFBBT(
  int numCols,
  const uint8_t *colFlags,
  const double *colLB,
  const double *colUB,
  const CoinPackedMatrix *matByRow,
  const char *rowSense,
  const double *rowRHS,
  const double *rowRange,
  double primalTolerance,
  double infinity)
  : newBounds_()
  , infeasible_(false)
  , infeasibleRow_(-1)
{
  // Working copies of bounds and flags updated during the pass.
  std::vector< double > curLB(colLB, colLB + numCols);
  std::vector< double > curUB(colUB, colUB + numCols);
  std::vector< uint8_t > flags(colFlags, colFlags + numCols);
  double *lb = curLB.data();
  double *ub = curUB.data();
  uint8_t *fl = flags.data();

  // Track which columns had their bounds improved (to build newBounds_ at end).
  // -1 = not touched, >= 0 = touched (index into the working arrays).
  std::vector< int > touched(static_cast< size_t >(numCols), -1);

  const int *idxs = matByRow->getIndices();
  const double *coefs = matByRow->getElements();
  const CoinBigIndex *start = matByRow->getVectorStarts();
  const int *length = matByRow->getVectorLengths();
  const int nRows = matByRow->getNumRows();

  double multipliers[2];
  double rhsAdjustments[2];

  for (int idxRow = 0; idxRow < nRows; ++idxRow) {
    const char sense = rowSense[idxRow];
    if (sense == 'N')
      continue;

    const CoinBigIndex rowStart = start[idxRow];
    const int rowLen = length[idxRow];
    const int *rowIdx = idxs + rowStart;
    const double *rowCoef = coefs + rowStart;
    const double range = rowRange ? rowRange[idxRow] : 0.0;

    const int numIter = CoinKnapsackRow::rowIterations(
      sense, rowRHS[idxRow], range, multipliers, rhsAdjustments);

    for (int it = 0; it < numIter; ++it) {
      const double mult = multipliers[it];
      // Effective constraint: Σ_j (mult * coef_j) * x_j ≤ b_eff
      double b_eff = mult * rhsAdjustments[it];

      // ── Pass 1: compute finite minimum activity, count unbounded contributors ──
      double minAct = 0.0;
      int nUnbounded = 0;
      int unboundedRowPos = -1; // position within the row arrays
      bool hasSemi = false;

      for (int k = 0; k < rowLen; ++k) {
        const int col = rowIdx[k];
        const double c = rowCoef[k] * mult;
        if (c == 0.0)
          continue;

        const uint8_t f = fl[col];

        if (f & COL_FIXED) {
          // Fixed variable: discount from effective RHS.
          b_eff -= c * lb[col];
          continue;
        }

        if (f & COL_SEMI) {
          hasSemi = true;
          break;
        }

        if (c > 0.0) {
          if (f & COL_LB_INF) {
            // Variable contributes -∞ to min activity.
            ++nUnbounded;
            unboundedRowPos = k;
            if (nUnbounded >= 2)
              break;
          } else {
            minAct += c * lb[col];
          }
        } else {
          // c < 0
          if (f & COL_UB_INF) {
            ++nUnbounded;
            unboundedRowPos = k;
            if (nUnbounded >= 2)
              break;
          } else {
            minAct += c * ub[col];
          }
        }
      }

      if (hasSemi || nUnbounded >= 2)
        continue;

      const double slack = b_eff - minAct;

      // ── Infeasibility check (only valid when min activity is finite) ──
      if (nUnbounded == 0 && slack < -primalTolerance) {
        infeasible_ = true;
        infeasibleRow_ = idxRow;
        return;
      }

      // ── Pass 2: tighten bounds ──
      if (nUnbounded == 0) {
        // Standard FBBT: tighten every non-fixed variable.
        for (int k = 0; k < rowLen; ++k) {
          const int col = rowIdx[k];
          const double c = rowCoef[k] * mult;
          if (c == 0.0)
            continue;

          const uint8_t f = fl[col];
          if (f & (COL_FIXED | COL_SEMI))
            continue;

          if (c > 0.0) {
            // Tighten upper bound: new_ub = lb[col] + slack / c
            double newUB = lb[col] + slack / c;
            if (f & COL_INTEGER)
              newUB = std::floor(newUB + primalTolerance);
            if (newUB < ub[col] - primalTolerance) {
              ub[col] = newUB;
              fl[col] = fl[col] & ~COL_UB_INF; // now finite
              touched[col] = col;
            }
          } else {
            // Tighten lower bound: new_lb = ub[col] + slack / c
            // (slack / c is ≤ 0 since c < 0 and slack ≥ 0 for feasible rows)
            double newLB = ub[col] + slack / c;
            if (f & COL_INTEGER)
              newLB = std::ceil(newLB - primalTolerance);
            if (newLB > lb[col] + primalTolerance) {
              lb[col] = newLB;
              fl[col] = fl[col] & ~COL_LB_INF; // now finite
              touched[col] = col;
            }
          }
        }
      } else {
        // nUnbounded == 1: tighten only the one unbounded variable.
        const int col = rowIdx[unboundedRowPos];
        const double c = rowCoef[unboundedRowPos] * mult;

        if (c > 0.0) {
          // lb[col] = -inf; derive new upper bound.
          // Constraint:  c * x_col + finite_others ≤ b_eff
          //   c * x_col ≤ b_eff - finite_others = slack  (since minAct = finite_others here)
          //   x_col ≤ slack / c
          double newUB = slack / c;
          if (fl[col] & COL_INTEGER)
            newUB = std::floor(newUB + primalTolerance);
          if (newUB < ub[col] - primalTolerance) {
            ub[col] = newUB;
            fl[col] = fl[col] & ~COL_UB_INF;
            touched[col] = col;
          }
        } else {
          // c < 0 and ub[col] = +inf; derive new lower bound.
          // c * x_col ≤ slack  =>  x_col ≥ slack / c  (flip sign: c < 0)
          double newLB = slack / c;
          if (fl[col] & COL_INTEGER)
            newLB = std::ceil(newLB - primalTolerance);
          if (newLB > lb[col] + primalTolerance) {
            lb[col] = newLB;
            fl[col] = fl[col] & ~COL_LB_INF;
            touched[col] = col;
          }
        }
      }
    } // for each row direction
  } // for each row

  // Collect tightened variables by comparing working bounds against originals.
  for (int j = 0; j < numCols; ++j) {
    if (touched[j] < 0)
      continue;
    if (lb[j] > colLB[j] + primalTolerance || ub[j] < colUB[j] - primalTolerance)
      newBounds_.push_back({ j, { lb[j], ub[j] } });
  }
}
