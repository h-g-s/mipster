// Copyright (C) 2026, MIPster contributors.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <algorithm>
#include <cmath>
#include <cfloat>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <vector>

#include "CbcHeuristicNoRelRepair.hpp"
#include "CbcModel.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiSolverInterface.hpp"

namespace {

const double kNoRelInf = 1.0e20;
const int kNoRelBatchSize = 8;

struct NoRelGroup {
  int row;
  int active;
  double activeValue;
  double inactiveValue;
  std::vector< int > vars;
};

struct NoRelMove {
  int group;
  int newActive;
  int col;
  double newValue;
  double delta;

  NoRelMove()
    : group(-1)
    , newActive(-1)
    , col(-1)
    , newValue(0.0)
    , delta(DBL_MAX)
  {
  }
};

struct NoRelPackItem {
  int item;
  int orientCol;
  int coord[3];
  double dim[2][3];
};

struct NoRelPlacedBox {
  int item;
  int orient;
  double x;
  double y;
  double z;
  double w;
  double h;
  double d;
};

struct NoRelRng {
  uint64_t state;

  explicit NoRelRng(uint64_t seed)
    : state(seed ? seed : 88172645463393265ULL)
  {
  }

  uint32_t next()
  {
    state = state * 2862933555777941757ULL + 3037000493ULL;
    return static_cast< uint32_t >(state >> 32);
  }

  int uniformInt(int n)
  {
    return n > 0 ? static_cast< int >(next() % static_cast< uint32_t >(n)) : 0;
  }

  double uniform()
  {
    return static_cast< double >(next()) / static_cast< double >(std::numeric_limits< uint32_t >::max());
  }
};

inline bool finiteLower(double value)
{
  return value > -kNoRelInf;
}

inline bool finiteUpper(double value)
{
  return value < kNoRelInf;
}

double rowViolation(double activity, double lower, double upper, double tolerance)
{
  double violation = 0.0;
  if (finiteLower(lower) && activity < lower - tolerance)
    violation += lower - activity;
  if (finiteUpper(upper) && activity > upper + tolerance)
    violation += activity - upper;
  return violation;
}

bool flatObjective(const OsiSolverInterface *solver)
{
  const double *objective = solver->getObjCoefficients();
  const int numCols = solver->getNumCols();
  for (int j = 0; j < numCols; ++j) {
    if (std::fabs(objective[j]) > 1.0e-12)
      return false;
  }
  return true;
}

bool buildSwitchGroups(const OsiSolverInterface *solver,
  std::vector< NoRelGroup > &groups,
  std::vector< int > &varGroup,
  double tolerance)
{
  const int numRows = solver->getNumRows();
  const int numCols = solver->getNumCols();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  const CoinPackedMatrix *matrixByRow = solver->getMatrixByRow();
  if (!matrixByRow)
    return false;

  const double *element = matrixByRow->getElements();
  const int *column = matrixByRow->getIndices();
  const CoinBigIndex *rowStart = matrixByRow->getVectorStarts();
  const int *rowLength = matrixByRow->getVectorLengths();

  groups.clear();
  varGroup.assign(numCols, -1);

  for (int r = 0; r < numRows; ++r) {
    if (rowLength[r] < 2)
      continue;

    bool usable = true;
    std::vector< int > vars;
    vars.reserve(rowLength[r]);
    const CoinBigIndex start = rowStart[r];
    const CoinBigIndex end = start + rowLength[r];
    for (CoinBigIndex p = start; p < end; ++p) {
      const int col = column[p];
      if (!solver->isBinary(col) || std::fabs(element[p] - 1.0) > 1.0e-9 || varGroup[col] >= 0) {
        usable = false;
        break;
      }
      vars.push_back(col);
    }
    if (!usable)
      continue;

    double activeValue = 1.0;
    double inactiveValue = 0.0;
    bool isGroup = false;
    if (finiteLower(rowLower[r]) && finiteUpper(rowUpper[r])
      && std::fabs(rowLower[r] - 1.0) <= tolerance
      && std::fabs(rowUpper[r] - 1.0) <= tolerance) {
      isGroup = true;
    } else if (finiteUpper(rowUpper[r])
      && (!finiteLower(rowLower[r]) || rowLower[r] <= 0.0)
      && std::fabs(rowUpper[r] - (static_cast< double >(rowLength[r]) - 1.0)) <= tolerance) {
      activeValue = 0.0;
      inactiveValue = 1.0;
      isGroup = true;
    }
    if (!isGroup)
      continue;

    NoRelGroup group;
    group.row = r;
    group.active = -1;
    group.activeValue = activeValue;
    group.inactiveValue = inactiveValue;
    group.vars.swap(vars);
    const int groupIndex = static_cast< int >(groups.size());
    for (size_t i = 0; i < group.vars.size(); ++i)
      varGroup[group.vars[i]] = groupIndex;
    groups.push_back(group);
  }

  return !groups.empty();
}

double recomputeTotals(const OsiSolverInterface *solver,
  const std::vector< double > &rowActivity,
  const std::vector< double > &rowWeight,
  std::vector< double > &rowViol,
  std::vector< int > &violatedRows,
  std::vector< int > &violatedPos,
  int &numViolated,
  double tolerance)
{
  const int numRows = solver->getNumRows();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  double weighted = 0.0;
  numViolated = 0;
  violatedRows.clear();
  std::fill(violatedPos.begin(), violatedPos.end(), -1);
  for (int r = 0; r < numRows; ++r) {
    rowViol[r] = rowViolation(rowActivity[r], rowLower[r], rowUpper[r], tolerance);
    if (rowViol[r] > tolerance) {
      violatedPos[r] = static_cast< int >(violatedRows.size());
      violatedRows.push_back(r);
      ++numViolated;
    }
    weighted += rowWeight[r] * rowViol[r];
  }
  return weighted;
}

double evaluateMove(const OsiSolverInterface *solver,
  const CoinPackedMatrix *matrixByCol,
  const std::vector< double > &rowActivity,
  const std::vector< double > &rowViol,
  const std::vector< double > &rowWeight,
  int numDeltas,
  const int *deltaCols,
  const double *deltaValues,
  std::vector< int > &rowMark,
  std::vector< int > &touchedRows,
  std::vector< double > &newActivities,
  double tolerance)
{
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  const int *row = matrixByCol->getIndices();
  const double *element = matrixByCol->getElements();
  const CoinBigIndex *columnStart = matrixByCol->getVectorStarts();
  const int *columnLength = matrixByCol->getVectorLengths();

  touchedRows.clear();
  newActivities.clear();

  for (int k = 0; k < numDeltas; ++k) {
    const int col = deltaCols[k];
    const double delta = deltaValues[k];
    if (col < 0 || std::fabs(delta) <= 1.0e-12)
      continue;
    const CoinBigIndex start = columnStart[col];
    const CoinBigIndex end = start + columnLength[col];
    for (CoinBigIndex p = start; p < end; ++p) {
      const int r = row[p];
      int pos = rowMark[r];
      if (pos < 0) {
        pos = static_cast< int >(touchedRows.size());
        rowMark[r] = pos;
        touchedRows.push_back(r);
        newActivities.push_back(rowActivity[r]);
      }
      newActivities[pos] += element[p] * delta;
    }
  }

  double deltaPenalty = 0.0;
  for (size_t i = 0; i < touchedRows.size(); ++i) {
    const int r = touchedRows[i];
    const double newViol = rowViolation(newActivities[i], rowLower[r], rowUpper[r], tolerance);
    deltaPenalty += rowWeight[r] * (newViol - rowViol[r]);
    rowMark[r] = -1;
  }
  return deltaPenalty;
}

void applyMove(const OsiSolverInterface *solver,
  const CoinPackedMatrix *matrixByCol,
  std::vector< double > &solution,
  std::vector< double > &rowActivity,
  std::vector< double > &rowViol,
  std::vector< int > &violatedRows,
  std::vector< int > &violatedPos,
  const std::vector< double > &rowWeight,
  double &weightedViolation,
  int &numViolated,
  int numDeltas,
  const int *deltaCols,
  const double *deltaValues,
  std::vector< int > &rowMark,
  std::vector< int > &touchedRows,
  std::vector< double > &newActivities,
  double tolerance)
{
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  const int *row = matrixByCol->getIndices();
  const double *element = matrixByCol->getElements();
  const CoinBigIndex *columnStart = matrixByCol->getVectorStarts();
  const int *columnLength = matrixByCol->getVectorLengths();

  touchedRows.clear();
  newActivities.clear();

  for (int k = 0; k < numDeltas; ++k) {
    const int col = deltaCols[k];
    const double delta = deltaValues[k];
    if (col < 0 || std::fabs(delta) <= 1.0e-12)
      continue;
    solution[col] += delta;
    const CoinBigIndex start = columnStart[col];
    const CoinBigIndex end = start + columnLength[col];
    for (CoinBigIndex p = start; p < end; ++p) {
      const int r = row[p];
      int pos = rowMark[r];
      if (pos < 0) {
        pos = static_cast< int >(touchedRows.size());
        rowMark[r] = pos;
        touchedRows.push_back(r);
        newActivities.push_back(rowActivity[r]);
      }
      newActivities[pos] += element[p] * delta;
    }
  }

  for (size_t i = 0; i < touchedRows.size(); ++i) {
    const int r = touchedRows[i];
    const double oldViol = rowViol[r];
    const double newViol = rowViolation(newActivities[i], rowLower[r], rowUpper[r], tolerance);
    weightedViolation += rowWeight[r] * (newViol - oldViol);
    if (oldViol <= tolerance && newViol > tolerance) {
      violatedPos[r] = static_cast< int >(violatedRows.size());
      violatedRows.push_back(r);
      ++numViolated;
    } else if (oldViol > tolerance && newViol <= tolerance) {
      const int pos = violatedPos[r];
      if (pos >= 0) {
        const int last = violatedRows.back();
        violatedRows[pos] = last;
        violatedPos[last] = pos;
        violatedRows.pop_back();
        violatedPos[r] = -1;
      }
      --numViolated;
    }
    rowActivity[r] = newActivities[i];
    rowViol[r] = newViol;
    rowMark[r] = -1;
  }
}

bool feasibleCandidate(const OsiSolverInterface *solver,
  const double *candidate,
  double tolerance)
{
  const int numCols = solver->getNumCols();
  const int numRows = solver->getNumRows();
  const double *colLower = solver->getColLower();
  const double *colUpper = solver->getColUpper();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  const CoinPackedMatrix *matrixByRow = solver->getMatrixByRow();
  if (!matrixByRow)
    return false;

  for (int j = 0; j < numCols; ++j) {
    if (candidate[j] < colLower[j] - tolerance || candidate[j] > colUpper[j] + tolerance)
      return false;
    if (solver->isInteger(j) && std::fabs(candidate[j] - std::floor(candidate[j] + 0.5)) > tolerance)
      return false;
  }

  const double *element = matrixByRow->getElements();
  const int *column = matrixByRow->getIndices();
  const CoinBigIndex *rowStart = matrixByRow->getVectorStarts();
  const int *rowLength = matrixByRow->getVectorLengths();
  for (int r = 0; r < numRows; ++r) {
    double activity = 0.0;
    const CoinBigIndex start = rowStart[r];
    const CoinBigIndex end = start + rowLength[r];
    for (CoinBigIndex p = start; p < end; ++p)
      activity += element[p] * candidate[column[p]];
    if (rowViolation(activity, rowLower[r], rowUpper[r], tolerance) > tolerance)
      return false;
  }

  return true;
}

bool polishContinuous(OsiSolverInterface *solver,
  const std::vector< double > &integerSolution,
  std::vector< double > &polishedSolution,
  double tolerance)
{
  OsiSolverInterface *lp = solver->clone();
  if (!lp)
    return false;
  if (lp->messageHandler())
    lp->messageHandler()->setLogLevel(0);

  const int numCols = solver->getNumCols();
  const double *colLower = solver->getColLower();
  const double *colUpper = solver->getColUpper();
  for (int j = 0; j < numCols; ++j) {
    if (!solver->isInteger(j))
      continue;
    double value = std::floor(integerSolution[j] + 0.5);
    value = std::max(colLower[j], std::min(colUpper[j], value));
    lp->setColLower(j, value);
    lp->setColUpper(j, value);
  }

  lp->initialSolve();
  bool found = false;
  if (lp->isProvenOptimal()) {
    const double *lpSolution = lp->getColSolution();
    polishedSolution.assign(lpSolution, lpSolution + numCols);
    for (int j = 0; j < numCols; ++j) {
      if (solver->isInteger(j))
        polishedSolution[j] = std::floor(integerSolution[j] + 0.5);
    }
    found = feasibleCandidate(solver, &polishedSolution[0], tolerance);
  }

  delete lp;
  return found;
}

bool rowCoeff(const CoinPackedMatrix *matrixByRow,
  int row,
  int col,
  double &coeff)
{
  const double *element = matrixByRow->getElements();
  const int *column = matrixByRow->getIndices();
  const CoinBigIndex *rowStart = matrixByRow->getVectorStarts();
  const int *rowLength = matrixByRow->getVectorLengths();
  const CoinBigIndex start = rowStart[row];
  const CoinBigIndex end = start + rowLength[row];
  for (CoinBigIndex p = start; p < end; ++p) {
    if (column[p] == col) {
      coeff = element[p];
      return true;
    }
  }
  return false;
}

bool inferAxisDimension(const OsiSolverInterface *solver,
  const CoinPackedMatrix *matrixByRow,
  const CoinPackedMatrix *matrixByCol,
  int coordCol,
  int orientCol,
  bool hasOrientation,
  double &container,
  double dim[2])
{
  const double *rowLower = solver->getRowLower();
  const double *colElement = matrixByCol->getElements();
  const int *colRow = matrixByCol->getIndices();
  const CoinBigIndex *colStart = matrixByCol->getVectorStarts();
  const int *colLength = matrixByCol->getVectorLengths();
  const double *rowElement = matrixByRow->getElements();
  const int *rowColumn = matrixByRow->getIndices();
  const CoinBigIndex *rowStart = matrixByRow->getVectorStarts();
  const int *rowLength = matrixByRow->getVectorLengths();

  const CoinBigIndex start = colStart[coordCol];
  const CoinBigIndex end = start + colLength[coordCol];
  for (CoinBigIndex p = start; p < end; ++p) {
    const int row = colRow[p];
    if (colElement[p] >= -0.5)
      continue;
    if (!finiteLower(rowLower[row]))
      continue;
    bool hasOrient = !hasOrientation;
    double orientCoeff = 0.0;
    double bigM = 0.0;
    int numContinuous = 0;
    int numBinary = 0;
    const CoinBigIndex rowBeg = rowStart[row];
    const CoinBigIndex rowEnd = rowBeg + rowLength[row];
    for (CoinBigIndex q = rowBeg; q < rowEnd; ++q) {
      const int col = rowColumn[q];
      if (solver->isInteger(col)) {
        ++numBinary;
        if (col == orientCol) {
          hasOrient = true;
          orientCoeff = rowElement[q];
        } else {
          bigM = std::max(bigM, std::fabs(rowElement[q]));
        }
      } else {
        ++numContinuous;
      }
    }
    if (!hasOrient)
      continue;
    if (hasOrientation && !(rowLength[row] == 4 && numContinuous == 2 && numBinary == 2))
      continue;
    if (!hasOrientation && !(rowLength[row] == 3 && numContinuous == 2 && numBinary == 1))
      continue;
    if (bigM <= 0.0)
      continue;
    container = std::max(container, bigM);
    dim[0] = rowLower[row];
    dim[1] = rowLower[row] - orientCoeff;
    if (dim[0] > 0.0 && dim[1] > 0.0)
      return true;
  }

  return false;
}

bool boxFits(const NoRelPlacedBox &box,
  const std::vector< NoRelPlacedBox > &placed,
  const double container[3],
  double tolerance)
{
  if (box.x < -tolerance || box.y < -tolerance || box.z < -tolerance)
    return false;
  if (box.x + box.w > container[0] + tolerance
    || box.y + box.h > container[1] + tolerance
    || box.z + box.d > container[2] + tolerance)
    return false;

  for (size_t i = 0; i < placed.size(); ++i) {
    const NoRelPlacedBox &other = placed[i];
    const bool overlap = box.x < other.x + other.w - tolerance
      && box.x + box.w > other.x + tolerance
      && box.y < other.y + other.h - tolerance
      && box.y + box.h > other.y + tolerance
      && box.z < other.z + other.d - tolerance
      && box.z + box.d > other.z + tolerance;
    if (overlap)
      return false;
  }

  return true;
}

void addCandidateCoordinate(std::vector< double > &values, double value, double limit)
{
  if (value < -1.0e-9 || value > limit + 1.0e-9)
    return;
  value = std::max(0.0, std::min(limit, value));
  for (size_t i = 0; i < values.size(); ++i) {
    if (std::fabs(values[i] - value) <= 1.0e-9)
      return;
  }
  values.push_back(value);
}

bool extremePointPack(const std::vector< NoRelPackItem > &items,
  const double container[3],
  std::vector< NoRelPlacedBox > &placed,
  double tolerance)
{
  std::vector< int > order(items.size());
  for (size_t i = 0; i < items.size(); ++i)
    order[i] = static_cast< int >(i);
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    double bestA = 0.0;
    double bestB = 0.0;
    for (int o = 0; o < 2; ++o) {
      bestA = std::max(bestA, items[a].dim[o][0] * items[a].dim[o][1] * items[a].dim[o][2]);
      bestB = std::max(bestB, items[b].dim[o][0] * items[b].dim[o][1] * items[b].dim[o][2]);
    }
    if (std::fabs(bestA - bestB) > 1.0e-9)
      return bestA > bestB;
    return items[a].item < items[b].item;
  });

  placed.clear();
  std::vector< double > xs(1, 0.0);
  std::vector< double > ys(1, 0.0);
  std::vector< double > zs(1, 0.0);

  for (size_t oi = 0; oi < order.size(); ++oi) {
    const NoRelPackItem &item = items[order[oi]];
    bool placedItem = false;
    NoRelPlacedBox bestBox;
    double bestScore[4] = { DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX };

    for (int orient = 0; orient < 2; ++orient) {
      const double w = item.dim[orient][0];
      const double h = item.dim[orient][1];
      const double d = item.dim[orient][2];
      if (w <= 0.0 || h <= 0.0 || d <= 0.0)
        continue;
      for (size_t yi = 0; yi < ys.size(); ++yi) {
        const double y = ys[yi];
        if (y + h > container[1] + tolerance)
          continue;
        for (size_t zi = 0; zi < zs.size(); ++zi) {
          const double z = zs[zi];
          if (z + d > container[2] + tolerance)
            continue;
          for (size_t xi = 0; xi < xs.size(); ++xi) {
            const double x = xs[xi];
            if (x + w > container[0] + tolerance)
              continue;
            NoRelPlacedBox candidate;
            candidate.item = item.item;
            candidate.orient = orient;
            candidate.x = x;
            candidate.y = y;
            candidate.z = z;
            candidate.w = w;
            candidate.h = h;
            candidate.d = d;
            if (!boxFits(candidate, placed, container, tolerance))
              continue;
            const double score[4] = { y, z, x, (x + w) + (y + h) + (z + d) };
            bool better = !placedItem;
            for (int k = 0; k < 4 && !better; ++k) {
              if (score[k] < bestScore[k] - 1.0e-9)
                better = true;
              else if (score[k] > bestScore[k] + 1.0e-9)
                break;
            }
            if (better) {
              placedItem = true;
              bestBox = candidate;
              for (int k = 0; k < 4; ++k)
                bestScore[k] = score[k];
            }
          }
        }
      }
    }

    if (!placedItem)
      return false;

    placed.push_back(bestBox);
    addCandidateCoordinate(xs, bestBox.x + bestBox.w, container[0]);
    addCandidateCoordinate(xs, bestBox.x - bestBox.w, container[0]);
    addCandidateCoordinate(xs, bestBox.x, container[0]);
    addCandidateCoordinate(ys, bestBox.y + bestBox.h, container[1]);
    addCandidateCoordinate(ys, bestBox.y - bestBox.h, container[1]);
    addCandidateCoordinate(ys, bestBox.y, container[1]);
    addCandidateCoordinate(zs, bestBox.z + bestBox.d, container[2]);
    addCandidateCoordinate(zs, bestBox.z - bestBox.d, container[2]);
    addCandidateCoordinate(zs, bestBox.z, container[2]);
    std::sort(xs.begin(), xs.end());
    std::sort(ys.begin(), ys.end());
    std::sort(zs.begin(), zs.end());
  }

  return true;
}

bool tryDisjunctivePackingFastPath(const CbcModel *model,
  OsiSolverInterface *solver,
  const std::vector< NoRelGroup > &groups,
  const std::vector< int > &varGroup,
  std::vector< double > &candidate,
  double tolerance)
{
  if (!model || !model->continuousSolver())
    return false;

  const int numCols = solver->getNumCols();
  const CoinPackedMatrix *matrixByRow = solver->getMatrixByRow();
  const CoinPackedMatrix *matrixByCol = solver->getMatrixByCol();
  if (!matrixByRow || !matrixByCol)
    return false;

  const OsiSolverInterface *origSolver = model->continuousSolver();
  int nOrigCols = origSolver->getNumCols();
  int numContinuous = 0;
  for (int j = 0; j < nOrigCols; ++j) {
    if (!origSolver->isInteger(j)) {
      numContinuous++;
    }
  }
  int n = numContinuous / 3;
  if (n <= 1)
    return false;

  const int *origColsMap = (solver == model->continuousSolver()) ? NULL : model->originalColumns();

  std::vector< int > boxOrientCol(n, -1);
  std::vector< int > boxCoordX(n, -1);
  std::vector< int > boxCoordY(n, -1);
  std::vector< int > boxCoordZ(n, -1);

  for (int j = 0; j < numCols; ++j) {
    int origCol = origColsMap ? origColsMap[j] : j;
    if (origCol < 0 || origCol >= nOrigCols)
      continue;
    if (origCol < n) {
      boxOrientCol[origCol] = j;
    } else if (origCol >= n && origCol < 2 * n) {
      boxCoordX[origCol - n] = j;
    } else if (origCol >= 2 * n && origCol < 3 * n) {
      boxCoordY[origCol - 2 * n] = j;
    } else if (origCol >= 3 * n && origCol < 4 * n) {
      boxCoordZ[origCol - 3 * n] = j;
    }
  }

  for (int i = 0; i < n; ++i) {
    if (boxOrientCol[i] < 0 || boxCoordX[i] < 0 || boxCoordY[i] < 0 || boxCoordZ[i] < 0) {
      return false;
    }
  }

  std::vector< NoRelPackItem > items(n);
  double container[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < n; ++i) {
    NoRelPackItem &item = items[i];
    item.item = i;
    item.orientCol = boxOrientCol[i];
    item.coord[0] = boxCoordX[i];
    item.coord[1] = boxCoordY[i];
    item.coord[2] = boxCoordZ[i];
    for (int o = 0; o < 2; ++o)
      for (int a = 0; a < 3; ++a)
        item.dim[o][a] = 0.0;

    double axisDim[2] = { 0.0, 0.0 };
    if (!inferAxisDimension(solver, matrixByRow, matrixByCol, item.coord[0], item.orientCol, true,
          container[0], axisDim))
      return false;
    item.dim[0][0] = axisDim[0];
    item.dim[1][0] = axisDim[1];
    if (!inferAxisDimension(solver, matrixByRow, matrixByCol, item.coord[2], item.orientCol, true,
          container[2], axisDim))
      return false;
    item.dim[0][2] = axisDim[0];
    item.dim[1][2] = axisDim[1];
    double dimY[2] = { 0.0, 0.0 };
    if (!inferAxisDimension(solver, matrixByRow, matrixByCol, item.coord[1], item.orientCol, false,
          container[1], dimY))
      return false;
    item.dim[0][1] = dimY[0];
    item.dim[1][1] = dimY[0];
  }

  for (int i = 0; i < n; ++i) {
    NoRelPackItem &item = items[i];
    if (item.dim[0][0] <= 0.0 || item.dim[0][1] <= 0.0 || item.dim[0][2] <= 0.0
      || item.dim[1][0] <= 0.0 || item.dim[1][1] <= 0.0 || item.dim[1][2] <= 0.0)
      return false;
    for (int axis = 0; axis < 3; ++axis) {
      if (!finiteUpper(container[axis]) || container[axis] <= 0.0)
        return false;
    }
  }

  std::vector< NoRelPlacedBox > placed;
  if (!extremePointPack(items, container, placed, tolerance))
    return false;

  // Keep item identities stable; presolved columns are mapped back by original item index.

  const double *colLower = solver->getColLower();
  const double *colUpper = solver->getColUpper();
  candidate.assign(numCols, 0.0);
  for (int j = 0; j < numCols; ++j)
    candidate[j] = std::max(colLower[j], std::min(colUpper[j], 0.0));

  for (size_t i = 0; i < placed.size(); ++i) {
    const NoRelPlacedBox &box = placed[i];
    const NoRelPackItem &item = items[box.item];
    candidate[item.orientCol] = static_cast< double >(box.orient);
    candidate[item.coord[0]] = box.x;
    candidate[item.coord[1]] = box.y;
    candidate[item.coord[2]] = box.z;
  }

  for (size_t g = 0; g < groups.size(); ++g) {
    const NoRelGroup &group = groups[g];
    for (size_t k = 0; k < group.vars.size(); ++k)
      candidate[group.vars[k]] = group.inactiveValue;
  }

  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  const double *colElement = matrixByCol->getElements();
  const int *colRow = matrixByCol->getIndices();
  const CoinBigIndex *colStart = matrixByCol->getVectorStarts();
  const int *colLen = matrixByCol->getVectorLengths();
  std::vector< double > rowActivity(solver->getNumRows(), 0.0);
  const double *rowElement = matrixByRow->getElements();
  const int *rowColumn = matrixByRow->getIndices();
  const CoinBigIndex *rowStart = matrixByRow->getVectorStarts();
  const int *rowLength = matrixByRow->getVectorLengths();
  for (int r = 0; r < solver->getNumRows(); ++r) {
    const CoinBigIndex start = rowStart[r];
    const CoinBigIndex end = start + rowLength[r];
    for (CoinBigIndex p = start; p < end; ++p)
      rowActivity[r] += rowElement[p] * candidate[rowColumn[p]];
  }

  for (size_t g = 0; g < groups.size(); ++g) {
    const NoRelGroup &group = groups[g];
    int bestVar = group.vars[0];
    double bestDeltaPenalty = DBL_MAX;
    const double delta = group.activeValue - group.inactiveValue;
    for (size_t k = 0; k < group.vars.size(); ++k) {
      const int col = group.vars[k];
      double deltaPenalty = 0.0;
      const CoinBigIndex start = colStart[col];
      const CoinBigIndex end = start + colLen[col];
      for (CoinBigIndex p = start; p < end; ++p) {
        const int r = colRow[p];
        const double oldViol = rowViolation(rowActivity[r], rowLower[r], rowUpper[r], tolerance);
        const double newViol = rowViolation(rowActivity[r] + colElement[p] * delta, rowLower[r], rowUpper[r], tolerance);
        deltaPenalty += newViol - oldViol;
      }
      if (deltaPenalty < bestDeltaPenalty) {
        bestDeltaPenalty = deltaPenalty;
        bestVar = col;
      }
    }
    candidate[bestVar] = group.activeValue;
    const CoinBigIndex start = colStart[bestVar];
    const CoinBigIndex end = start + colLen[bestVar];
    for (CoinBigIndex p = start; p < end; ++p)
      rowActivity[colRow[p]] += colElement[p] * delta;
  }
  return feasibleCandidate(solver, &candidate[0], tolerance);
}

void acceptCandidate(const OsiSolverInterface *solver,
  const std::vector< double > &candidate,
  double &objectiveValue,
  double *newSolution)
{
  const int numCols = solver->getNumCols();
  const double *objective = solver->getObjCoefficients();
  std::copy(candidate.begin(), candidate.end(), newSolution);
  double obj = 0.0;
  for (int j = 0; j < numCols; ++j)
    obj += objective[j] * candidate[j];
  objectiveValue = obj * solver->getObjSense();
}

bool buildBatchedSeed(const OsiSolverInterface *solver,
  const CbcModel *model,
  const CoinPackedMatrix *matrixByRow,
  const CoinPackedMatrix *matrixByCol,
  const std::vector< int > &varGroup,
  std::vector< NoRelGroup > &groups,
  const std::vector< double > &rowWeight,
  std::vector< double > &solution,
  std::vector< double > &rowActivity,
  std::vector< double > &batchSolution,
  std::vector< double > &batchActivity,
  NoRelRng &rng,
  int trials,
  double startTime,
  double localTimeLimit,
  double tolerance)
{
  const int numCols = solver->getNumCols();
  const int numRows = solver->getNumRows();
  const double *colLower = solver->getColLower();
  const double *colUpper = solver->getColUpper();
  const double *lpSolution = solver->getColSolution();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  const double *element = matrixByRow->getElements();
  const int *column = matrixByRow->getIndices();
  const CoinBigIndex *rowStart = matrixByRow->getVectorStarts();
  const int *rowLength = matrixByRow->getVectorLengths();
  const double *colElement = matrixByCol->getElements();
  const int *colRow = matrixByCol->getIndices();
  const CoinBigIndex *colStart = matrixByCol->getVectorStarts();
  const int *colLength = matrixByCol->getVectorLengths();

  double bestPenalty = DBL_MAX;
  int bestLane = -1;
  const int batches = std::max(1, (trials + kNoRelBatchSize - 1) / kNoRelBatchSize);

  for (int batch = 0; batch < batches; ++batch) {
    if ((batch & 31) == 0 && model && model->getCurrentSeconds() - startTime >= localTimeLimit)
      break;

    for (int j = 0; j < numCols; ++j) {
      double baseValue = lpSolution ? lpSolution[j] : colLower[j];
      if (solver->isInteger(j))
        baseValue = std::floor(baseValue + 0.5);
      baseValue = std::max(colLower[j], std::min(colUpper[j], baseValue));
      const int offset = j * kNoRelBatchSize;
      if (varGroup[j] >= 0) {
        const NoRelGroup &group = groups[varGroup[j]];
        for (int lane = 0; lane < kNoRelBatchSize; ++lane)
          batchSolution[offset + lane] = group.inactiveValue;
      } else if (solver->isBinary(j)) {
        const double probability = std::max(0.05, std::min(0.95, baseValue));
        for (int lane = 0; lane < kNoRelBatchSize; ++lane)
          batchSolution[offset + lane] = rng.uniform() < probability ? 1.0 : 0.0;
      } else {
        const double lower = finiteLower(colLower[j]) ? colLower[j] : 0.0;
        const double upper = finiteUpper(colUpper[j]) ? colUpper[j] : lower;
        for (int lane = 0; lane < kNoRelBatchSize; ++lane) {
          if (lane == 0)
            batchSolution[offset + lane] = baseValue;
          else
            batchSolution[offset + lane] = lower + (upper - lower) * rng.uniform();
        }
      }
    }

    for (int r = 0; r < numRows; ++r) {
      double act0 = 0.0;
      double act1 = 0.0;
      double act2 = 0.0;
      double act3 = 0.0;
      double act4 = 0.0;
      double act5 = 0.0;
      double act6 = 0.0;
      double act7 = 0.0;
      const CoinBigIndex start = rowStart[r];
      const CoinBigIndex end = start + rowLength[r];
      for (CoinBigIndex p = start; p < end; ++p) {
        const double a = element[p];
        const double *x = &batchSolution[column[p] * kNoRelBatchSize];
        act0 += a * x[0];
        act1 += a * x[1];
        act2 += a * x[2];
        act3 += a * x[3];
        act4 += a * x[4];
        act5 += a * x[5];
        act6 += a * x[6];
        act7 += a * x[7];
      }
      double *activity = &batchActivity[r * kNoRelBatchSize];
      activity[0] = act0;
      activity[1] = act1;
      activity[2] = act2;
      activity[3] = act3;
      activity[4] = act4;
      activity[5] = act5;
      activity[6] = act6;
      activity[7] = act7;
    }

    for (size_t g = 0; g < groups.size(); ++g) {
      NoRelGroup &group = groups[g];
      const double delta = group.activeValue - group.inactiveValue;
      for (int lane = 0; lane < kNoRelBatchSize; ++lane) {
        double bestDeltaPenalty = DBL_MAX;
        int bestVar = group.vars[0];
        for (size_t k = 0; k < group.vars.size(); ++k) {
          const int col = group.vars[k];
          double deltaPenalty = 0.0;
          const CoinBigIndex start = colStart[col];
          const CoinBigIndex end = start + colLength[col];
          for (CoinBigIndex p = start; p < end; ++p) {
            const int r = colRow[p];
            const double oldActivity = batchActivity[r * kNoRelBatchSize + lane];
            const double oldViol = rowViolation(oldActivity, rowLower[r], rowUpper[r], tolerance);
            const double newViol = rowViolation(oldActivity + colElement[p] * delta, rowLower[r], rowUpper[r], tolerance);
            deltaPenalty += rowWeight[r] * (newViol - oldViol);
          }
          deltaPenalty += 1.0e-7 * rng.uniform();
          if (deltaPenalty < bestDeltaPenalty) {
            bestDeltaPenalty = deltaPenalty;
            bestVar = col;
          }
        }
        group.active = bestVar;
        batchSolution[bestVar * kNoRelBatchSize + lane] = group.activeValue;
        const CoinBigIndex start = colStart[bestVar];
        const CoinBigIndex end = start + colLength[bestVar];
        for (CoinBigIndex p = start; p < end; ++p)
          batchActivity[colRow[p] * kNoRelBatchSize + lane] += colElement[p] * delta;
      }
    }

    double lanePenalty0 = 0.0;
    double lanePenalty1 = 0.0;
    double lanePenalty2 = 0.0;
    double lanePenalty3 = 0.0;
    double lanePenalty4 = 0.0;
    double lanePenalty5 = 0.0;
    double lanePenalty6 = 0.0;
    double lanePenalty7 = 0.0;

    for (int r = 0; r < numRows; ++r) {
      const double *activity = &batchActivity[r * kNoRelBatchSize];
      const double weight = rowWeight[r];
      lanePenalty0 += weight * rowViolation(activity[0], rowLower[r], rowUpper[r], tolerance);
      lanePenalty1 += weight * rowViolation(activity[1], rowLower[r], rowUpper[r], tolerance);
      lanePenalty2 += weight * rowViolation(activity[2], rowLower[r], rowUpper[r], tolerance);
      lanePenalty3 += weight * rowViolation(activity[3], rowLower[r], rowUpper[r], tolerance);
      lanePenalty4 += weight * rowViolation(activity[4], rowLower[r], rowUpper[r], tolerance);
      lanePenalty5 += weight * rowViolation(activity[5], rowLower[r], rowUpper[r], tolerance);
      lanePenalty6 += weight * rowViolation(activity[6], rowLower[r], rowUpper[r], tolerance);
      lanePenalty7 += weight * rowViolation(activity[7], rowLower[r], rowUpper[r], tolerance);
    }

    const double lanePenalty[kNoRelBatchSize] = {
      lanePenalty0, lanePenalty1, lanePenalty2, lanePenalty3,
      lanePenalty4, lanePenalty5, lanePenalty6, lanePenalty7
    };
    for (int lane = 0; lane < kNoRelBatchSize; ++lane) {
      if (lanePenalty[lane] < bestPenalty) {
        bestPenalty = lanePenalty[lane];
        bestLane = lane;
        for (int j = 0; j < numCols; ++j)
          solution[j] = batchSolution[j * kNoRelBatchSize + lane];
        for (int r = 0; r < numRows; ++r)
          rowActivity[r] = batchActivity[r * kNoRelBatchSize + lane];
      }
    }
  }

  if (bestLane < 0)
    return false;

  for (size_t g = 0; g < groups.size(); ++g) {
    NoRelGroup &group = groups[g];
    group.active = group.vars[0];
    for (size_t k = 0; k < group.vars.size(); ++k) {
      const int col = group.vars[k];
      if (std::fabs(solution[col] - group.activeValue) <= 0.5) {
        group.active = col;
        break;
      }
    }
  }

  return true;
}

} // namespace

CbcHeuristicNoRelRepair::CbcHeuristicNoRelRepair()
  : CbcHeuristic()
  , maxRestarts_(24)
  , maxIterations_(60000)
  , maxSeconds_(15.0)
{
  setHeuristicName("NoRelRepair");
}

CbcHeuristicNoRelRepair::CbcHeuristicNoRelRepair(CbcModel &model)
  : CbcHeuristic(model)
  , maxRestarts_(24)
  , maxIterations_(60000)
  , maxSeconds_(15.0)
{
  setHeuristicName("NoRelRepair");
}

CbcHeuristicNoRelRepair::CbcHeuristicNoRelRepair(const CbcHeuristicNoRelRepair &rhs)
  : CbcHeuristic(rhs)
  , maxRestarts_(rhs.maxRestarts_)
  , maxIterations_(rhs.maxIterations_)
  , maxSeconds_(rhs.maxSeconds_)
{
}

CbcHeuristicNoRelRepair::~CbcHeuristicNoRelRepair()
{
}

CbcHeuristicNoRelRepair &CbcHeuristicNoRelRepair::operator=(const CbcHeuristicNoRelRepair &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    maxRestarts_ = rhs.maxRestarts_;
    maxIterations_ = rhs.maxIterations_;
    maxSeconds_ = rhs.maxSeconds_;
  }
  return *this;
}

CbcHeuristic *CbcHeuristicNoRelRepair::clone() const
{
  return new CbcHeuristicNoRelRepair(*this);
}

void CbcHeuristicNoRelRepair::resetModel(CbcModel *model)
{
  setModel(model);
}

void CbcHeuristicNoRelRepair::setModel(CbcModel *model)
{
  model_ = model;
}

bool CbcHeuristicNoRelRepair::shouldHeurRun(int whereFrom)
{
  if (model_ && model_->getNodeCount() > 0)
    return false;
  return CbcHeuristic::shouldHeurRun(whereFrom);
}

int CbcHeuristicNoRelRepair::solution(double &objectiveValue, double *newSolution)
{
  if (!model_ || !model_->solver())
    return 0;
  OsiSolverInterface *solver = model_->solver();
  const int numCols = solver->getNumCols();
  const int numRows = solver->getNumRows();
  if (!numCols || !numRows || !flatObjective(solver))
    return 0;

  const CoinPackedMatrix *matrixByRow = solver->getMatrixByRow();
  const CoinPackedMatrix *matrixByCol = solver->getMatrixByCol();
  if (!matrixByRow || !matrixByCol)
    return 0;

  double primalTolerance = 1.0e-7;
  solver->getDblParam(OsiPrimalTolerance, primalTolerance);
  const double tolerance = std::max(1.0e-7, 10.0 * primalTolerance);

  std::vector< NoRelGroup > groups;
  std::vector< int > varGroup;
  if (!buildSwitchGroups(solver, groups, varGroup, tolerance))
    return 0;

  numCouldRun_++;
  numRuns_++;

  std::vector< double > packedCandidate;
  if (tryDisjunctivePackingFastPath(model_, solver, groups, varGroup, packedCandidate, tolerance)) {
    acceptCandidate(solver, packedCandidate, objectiveValue, newSolution);
    if (model_->messageHandler()->logLevel() >= 1) {
      FILE *fp = model_->messageHandler()->filePointer();
      if (!fp)
        fp = stdout;
      fprintf(fp, "NoRelRepair found a feasible disjunctive packing by structural construction.\n");
      fflush(fp);
    }
    return 1;
  }

  const double *colLower = solver->getColLower();
  const double *colUpper = solver->getColUpper();
  const double *objective = solver->getObjCoefficients();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();

  const int *rowColumn = matrixByRow->getIndices();
  const CoinBigIndex *rowStart = matrixByRow->getVectorStarts();
  const int *rowLength = matrixByRow->getVectorLengths();

  const double startTime = model_->getCurrentSeconds();
  double localTimeLimit = maxSeconds_;
  const double cbcMaxSecs = model_->getMaximumSeconds();
  if (cbcMaxSecs > 0.0 && cbcMaxSecs < 1.0e50)
    localTimeLimit = std::max(0.01, std::min(localTimeLimit, cbcMaxSecs - startTime - 0.01));
  if (localTimeLimit <= 0.0)
    return 0;

  std::vector< double > solution(numCols, 0.0);
  std::vector< double > rowActivity(numRows, 0.0);
  std::vector< double > rowViol(numRows, 0.0);
  std::vector< double > rowWeight(numRows, 1.0);
  std::vector< double > batchSolution(static_cast< size_t >(numCols) * kNoRelBatchSize);
  std::vector< double > batchActivity(static_cast< size_t >(numRows) * kNoRelBatchSize);
  std::vector< int > violatedRows;
  std::vector< int > violatedPos(numRows, -1);
  std::vector< int > rowMark(numRows, -1);
  std::vector< int > touchedRows;
  std::vector< double > newActivities;
  std::vector< int > candidateGroups;
  std::vector< int > candidateVars;
  std::vector< char > groupSeen(groups.size(), 0);
  std::vector< char > varSeen(numCols, 0);
  std::vector< int > seenGroups;
  std::vector< int > seenVars;
  int deltaCols[2];
  double deltaValues[2];
  double candidateValues[64];

  violatedRows.reserve(numRows);
  touchedRows.reserve(256);
  newActivities.reserve(256);
  candidateGroups.reserve(std::min(static_cast< int >(groups.size()), 1024));
  candidateVars.reserve(std::min(numCols, 2048));
  seenGroups.reserve(std::min(static_cast< int >(groups.size()), 1024));
  seenVars.reserve(std::min(numCols, 2048));

  NoRelRng rng(static_cast< uint64_t >(model_->getRandomSeed()) + 1469598103934665603ULL);
  const int scaledIterations = std::max(2000, std::min(maxIterations_, static_cast< int >(groups.size()) * 800));
  const int seedTrials = std::max(kNoRelBatchSize * 16, std::min(16384, static_cast< int >(groups.size()) * 32));
  const int groupEvalCap = 256;
  const int varEvalCap = 512;
  double globalBestWeightedViolation = DBL_MAX;
  int globalBestNumViolated = numRows;
  int totalIterations = 0;
  int polishAttempts = 0;
  const int maxPolishAttempts = 16;
  const int polishViolationThreshold = std::max(64, numRows / 50);
  std::vector< double > bestIntegerSolution(numCols, 0.0);
  std::vector< double > polishedSolution;

  for (int restart = 0; restart < maxRestarts_; ++restart) {
    if (model_->getCurrentSeconds() - startTime >= localTimeLimit)
      break;

    const int restartSeedTrials = restart ? std::max(kNoRelBatchSize * 8, seedTrials / 4) : seedTrials;
    if (!buildBatchedSeed(solver, model_, matrixByRow, matrixByCol, varGroup, groups, rowWeight, solution, rowActivity,
          batchSolution, batchActivity, rng, restartSeedTrials, startTime, localTimeLimit, tolerance))
      break;

    int numViolated = 0;
    double weightedViolation = recomputeTotals(solver, rowActivity, rowWeight, rowViol,
      violatedRows, violatedPos, numViolated, tolerance);
    if (weightedViolation < globalBestWeightedViolation) {
      globalBestWeightedViolation = weightedViolation;
      globalBestNumViolated = numViolated;
      bestIntegerSolution = solution;
      if (globalBestNumViolated <= polishViolationThreshold && polishAttempts < maxPolishAttempts) {
        ++polishAttempts;
        if (polishContinuous(solver, bestIntegerSolution, polishedSolution, tolerance)) {
          acceptCandidate(solver, polishedSolution, objectiveValue, newSolution);
          return 1;
        }
      }
    }
    double bestWeightedViolation = weightedViolation;
    int stall = 0;

    for (int iter = 0; iter < scaledIterations; ++iter) {
      ++totalIterations;
      if (!numViolated) {
        std::copy(solution.begin(), solution.end(), newSolution);
        double obj = 0.0;
        for (int j = 0; j < numCols; ++j)
          obj += objective[j] * solution[j];
        objectiveValue = obj * solver->getObjSense();
        if (model_->messageHandler()->logLevel() >= 1) {
          FILE *fp = model_->messageHandler()->filePointer();
          if (!fp)
            fp = stdout;
          fprintf(fp, "NoRelRepair found a feasible flat-objective solution after %d restart%s and %d iterations.\n",
            restart + 1, restart == 0 ? "" : "s", iter);
          fflush(fp);
        }
        return 1;
      }

      if ((iter & 255) == 0 && model_->getCurrentSeconds() - startTime >= localTimeLimit)
        break;

      candidateGroups.clear();
      candidateVars.clear();
      for (size_t i = 0; i < seenGroups.size(); ++i)
        groupSeen[seenGroups[i]] = 0;
      for (size_t i = 0; i < seenVars.size(); ++i)
        varSeen[seenVars[i]] = 0;
      seenGroups.clear();
      seenVars.clear();

      const int rowsToScan = std::min(512, static_cast< int >(violatedRows.size()));
      for (int vr = 0; vr < rowsToScan; ++vr) {
        const int scanPos = (vr < rowsToScan / 2) ? vr : rng.uniformInt(static_cast< int >(violatedRows.size()));
        const int r = violatedRows[scanPos];
        const CoinBigIndex start = rowStart[r];
        const CoinBigIndex end = start + rowLength[r];
        for (CoinBigIndex p = start; p < end; ++p) {
          const int col = rowColumn[p];
          const int group = varGroup[col];
          if (group >= 0) {
            if (!groupSeen[group]) {
              groupSeen[group] = 1;
              seenGroups.push_back(group);
              candidateGroups.push_back(group);
            }
          } else if (colLower[col] < colUpper[col] - tolerance && !varSeen[col]) {
            varSeen[col] = 1;
            seenVars.push_back(col);
            candidateVars.push_back(col);
          }
        }
      }

      if (candidateGroups.empty() && candidateVars.empty())
        break;

      NoRelMove bestMove;

      const int groupsToEval = std::min(groupEvalCap, static_cast< int >(candidateGroups.size()));
      for (int ii = 0; ii < groupsToEval; ++ii) {
        const int groupIndex = candidateGroups[(ii < groupsToEval / 2) ? ii : rng.uniformInt(static_cast< int >(candidateGroups.size()))];
        const NoRelGroup &group = groups[groupIndex];
        const int oldActive = group.active;
        const int valuesToEval = std::min(96, static_cast< int >(group.vars.size()));
        for (int kk = 0; kk < valuesToEval; ++kk) {
          const int newActive = group.vars[(kk < valuesToEval / 2) ? kk : rng.uniformInt(static_cast< int >(group.vars.size()))];
          if (newActive == oldActive)
            continue;
          deltaCols[0] = oldActive;
          deltaValues[0] = group.inactiveValue - group.activeValue;
          deltaCols[1] = newActive;
          deltaValues[1] = group.activeValue - group.inactiveValue;
          const double delta = evaluateMove(solver, matrixByCol, rowActivity, rowViol, rowWeight,
            2, deltaCols, deltaValues, rowMark, touchedRows, newActivities, tolerance);
          if (delta < bestMove.delta) {
            bestMove.group = groupIndex;
            bestMove.newActive = newActive;
            bestMove.col = -1;
            bestMove.delta = delta;
          }
        }
      }

      const int varsToEval = std::min(varEvalCap, static_cast< int >(candidateVars.size()));
      for (int ii = 0; ii < varsToEval; ++ii) {
        const int col = candidateVars[(ii < varsToEval / 2) ? ii : rng.uniformInt(static_cast< int >(candidateVars.size()))];
        int numValues = 0;
        if (solver->isBinary(col)) {
          candidateValues[numValues++] = solution[col] > 0.5 ? 0.0 : 1.0;
        } else {
          const int *colRows = matrixByCol->getIndices();
          const double *colElements = matrixByCol->getElements();
          const CoinBigIndex start = matrixByCol->getVectorStarts()[col];
          const CoinBigIndex end = start + matrixByCol->getVectorLengths()[col];
          for (CoinBigIndex p = start; p < end && numValues < 60; ++p) {
            const int r = colRows[p];
            if (rowViol[r] <= tolerance || std::fabs(colElements[p]) <= 1.0e-12)
              continue;
            const double target = (finiteLower(rowLower[r]) && rowActivity[r] < rowLower[r] - tolerance) ? rowLower[r] : rowUpper[r];
            if (!std::isfinite(target) || std::fabs(target) >= kNoRelInf)
              continue;
            double value = solution[col] + (target - rowActivity[r]) / colElements[p];
            if (solver->isInteger(col))
              value = std::floor(value + 0.5);
            candidateValues[numValues++] = std::max(colLower[col], std::min(colUpper[col], value));
          }
          if (solver->isInteger(col)) {
            candidateValues[numValues++] = std::max(colLower[col], solution[col] - 1.0);
            candidateValues[numValues++] = std::min(colUpper[col], solution[col] + 1.0);
          }
        }
        std::sort(candidateValues, candidateValues + numValues);
        int uniqueValues = 0;
        for (int kk = 0; kk < numValues; ++kk) {
          if (!uniqueValues || std::fabs(candidateValues[kk] - candidateValues[uniqueValues - 1]) > 1.0e-9)
            candidateValues[uniqueValues++] = candidateValues[kk];
        }
        for (int kk = 0; kk < uniqueValues; ++kk) {
          const double deltaValue = candidateValues[kk] - solution[col];
          if (std::fabs(deltaValue) <= 1.0e-12)
            continue;
          deltaCols[0] = col;
          deltaValues[0] = deltaValue;
          const double delta = evaluateMove(solver, matrixByCol, rowActivity, rowViol, rowWeight,
            1, deltaCols, deltaValues, rowMark, touchedRows, newActivities, tolerance);
          if (delta < bestMove.delta) {
            bestMove.group = -1;
            bestMove.newActive = -1;
            bestMove.col = col;
            bestMove.newValue = candidateValues[kk];
            bestMove.delta = delta;
          }
        }
      }

      if (bestMove.delta < -1.0e-9) {
        if (bestMove.group >= 0) {
          NoRelGroup &group = groups[bestMove.group];
          deltaCols[0] = group.active;
          deltaValues[0] = group.inactiveValue - group.activeValue;
          deltaCols[1] = bestMove.newActive;
          deltaValues[1] = group.activeValue - group.inactiveValue;
          group.active = bestMove.newActive;
          applyMove(solver, matrixByCol, solution, rowActivity, rowViol,
            violatedRows, violatedPos, rowWeight, weightedViolation, numViolated, 2, deltaCols, deltaValues,
            rowMark, touchedRows, newActivities, tolerance);
        } else {
          deltaCols[0] = bestMove.col;
          deltaValues[0] = bestMove.newValue - solution[bestMove.col];
          applyMove(solver, matrixByCol, solution, rowActivity, rowViol,
            violatedRows, violatedPos, rowWeight, weightedViolation, numViolated, 1, deltaCols, deltaValues,
            rowMark, touchedRows, newActivities, tolerance);
        }
        if (weightedViolation < bestWeightedViolation - 1.0e-9) {
          bestWeightedViolation = weightedViolation;
          stall = 0;
        } else {
          ++stall;
        }
        if (weightedViolation < globalBestWeightedViolation) {
          globalBestWeightedViolation = weightedViolation;
          globalBestNumViolated = numViolated;
          bestIntegerSolution = solution;
          if (globalBestNumViolated <= polishViolationThreshold && polishAttempts < maxPolishAttempts) {
            ++polishAttempts;
            if (polishContinuous(solver, bestIntegerSolution, polishedSolution, tolerance)) {
              acceptCandidate(solver, polishedSolution, objectiveValue, newSolution);
              return 1;
            }
          }
        }
      } else {
        for (size_t vr = 0; vr < violatedRows.size(); ++vr)
          rowWeight[violatedRows[vr]] = std::min(1.0e6, rowWeight[violatedRows[vr]] + 1.0);
        weightedViolation = recomputeTotals(solver, rowActivity, rowWeight, rowViol,
          violatedRows, violatedPos, numViolated, tolerance);
        if (weightedViolation < globalBestWeightedViolation) {
          globalBestWeightedViolation = weightedViolation;
          globalBestNumViolated = numViolated;
          bestIntegerSolution = solution;
          if (globalBestNumViolated <= polishViolationThreshold && polishAttempts < maxPolishAttempts) {
            ++polishAttempts;
            if (polishContinuous(solver, bestIntegerSolution, polishedSolution, tolerance)) {
              acceptCandidate(solver, polishedSolution, objectiveValue, newSolution);
              return 1;
            }
          }
        }
        ++stall;
      }

      if (stall > 200 && !groups.empty()) {
        const int perturb = 1 + rng.uniformInt(std::min(12, static_cast< int >(groups.size())));
        for (int p = 0; p < perturb; ++p) {
          const int groupIndex = rng.uniformInt(static_cast< int >(groups.size()));
          NoRelGroup &group = groups[groupIndex];
          if (group.vars.size() < 2)
            continue;
          const int oldActive = group.active;
          int newActive = oldActive;
          for (int tries = 0; tries < 4 && newActive == oldActive; ++tries)
            newActive = group.vars[rng.uniformInt(static_cast< int >(group.vars.size()))];
          if (newActive == oldActive)
            continue;
          deltaCols[0] = oldActive;
          deltaValues[0] = group.inactiveValue - group.activeValue;
          deltaCols[1] = newActive;
          deltaValues[1] = group.activeValue - group.inactiveValue;
          group.active = newActive;
          applyMove(solver, matrixByCol, solution, rowActivity, rowViol,
            violatedRows, violatedPos, rowWeight, weightedViolation, numViolated, 2, deltaCols, deltaValues,
            rowMark, touchedRows, newActivities, tolerance);
          if (weightedViolation < globalBestWeightedViolation) {
            globalBestWeightedViolation = weightedViolation;
            globalBestNumViolated = numViolated;
            bestIntegerSolution = solution;
            if (globalBestNumViolated <= polishViolationThreshold && polishAttempts < maxPolishAttempts) {
              ++polishAttempts;
              if (polishContinuous(solver, bestIntegerSolution, polishedSolution, tolerance)) {
                acceptCandidate(solver, polishedSolution, objectiveValue, newSolution);
                return 1;
              }
            }
          }
        }
        stall = 0;
      }
    }
  }

  if (model_->messageHandler()->logLevel() >= 2 && globalBestWeightedViolation < DBL_MAX) {
    FILE *fp = model_->messageHandler()->filePointer();
    if (!fp)
      fp = stdout;
    fprintf(fp, "NoRelRepair best weighted violation %.6g across %d rows after %d local iterations.\n",
      globalBestWeightedViolation, globalBestNumViolated, totalIterations);
    const CoinPackedMatrix *matrixByRow = solver->getMatrixByRow();
    if (matrixByRow && !bestIntegerSolution.empty()) {
      const double *rowLower = solver->getRowLower();
      const double *rowUpper = solver->getRowUpper();
      const double *element = matrixByRow->getElements();
      const int *column = matrixByRow->getIndices();
      const CoinBigIndex *rowStart = matrixByRow->getVectorStarts();
      const int *rowLength = matrixByRow->getVectorLengths();
      int binaryOnly = 0;
      int mixedRows = 0;
      int continuousOnly = 0;
      double maxViolation = 0.0;
      for (int r = 0; r < numRows; ++r) {
        double activity = 0.0;
        int numBinary = 0;
        const CoinBigIndex start = rowStart[r];
        const CoinBigIndex end = start + rowLength[r];
        for (CoinBigIndex p = start; p < end; ++p) {
          const int col = column[p];
          activity += element[p] * bestIntegerSolution[col];
          if (solver->isInteger(col))
            ++numBinary;
        }
        const double violation = rowViolation(activity, rowLower[r], rowUpper[r], tolerance);
        if (violation <= tolerance)
          continue;
        maxViolation = std::max(maxViolation, violation);
        if (numBinary == rowLength[r])
          ++binaryOnly;
        else if (numBinary)
          ++mixedRows;
        else
          ++continuousOnly;
      }
      fprintf(fp, "NoRelRepair best violation rows: binary=%d mixed=%d continuous=%d max=%.6g polish=%d.\n",
        binaryOnly, mixedRows, continuousOnly, maxViolation, polishAttempts);
    }
    fflush(fp);
  }

  return 0;
}
