#include <math.h>
#include <stdio.h>

#include "Cbc_C_Interface.h"
#include "test_utils.h"

static Cbc_Model *build_pair_blocks(int blocks, int flat_objective)
{
  Cbc_Model *m = Cbc_newModel();
  char name[64];

  for (int b = 0; b < blocks; ++b) {
    for (int j = 0; j < 4; ++j) {
      snprintf(name, sizeof(name), "x_%d_%d", b, j);
      Cbc_addCol(m, name, 0.0, 1.0, (!flat_objective && b == 0 && j == 0) ? 1.0 : 0.0, 1, 0, NULL, NULL);
    }
  }

  for (int b = 0; b < blocks; ++b) {
    const int base = 4 * b;
    int idx[2];
    double coef[2] = { 1.0, 1.0 };

    idx[0] = base;
    idx[1] = base + 1;
    snprintf(name, sizeof(name), "choose_a_%d", b);
    Cbc_addRow(m, name, 2, idx, coef, 'E', 1.0);

    idx[0] = base + 2;
    idx[1] = base + 3;
    snprintf(name, sizeof(name), "choose_b_%d", b);
    Cbc_addRow(m, name, 2, idx, coef, 'E', 1.0);

    idx[0] = base;
    idx[1] = base + 2;
    snprintf(name, sizeof(name), "cap_0_%d", b);
    Cbc_addRow(m, name, 2, idx, coef, 'L', 1.0);

    idx[0] = base + 1;
    idx[1] = base + 3;
    snprintf(name, sizeof(name), "cap_1_%d", b);
    Cbc_addRow(m, name, 2, idx, coef, 'L', 1.0);
  }

  return m;
}

static Cbc_Model *build_feasibility_probe(void)
{
  Cbc_Model *m = Cbc_newModel();
  int idx[3];
  double coef[3];

  Cbc_addCol(m, "x", 0.0, 10.0, 0.0, 0, 0, NULL, NULL);
  Cbc_addCol(m, "y", 0.0, 2.0, 0.0, 1, 0, NULL, NULL);
  Cbc_addCol(m, "z", 0.0, 1.0, 0.0, 1, 0, NULL, NULL);

  idx[0] = 0;
  idx[1] = 1;
  idx[2] = 2;
  coef[0] = 1.0;
  coef[1] = 1.0;
  coef[2] = 1.0;
  Cbc_addRow(m, "demand", 3, idx, coef, 'G', 2.0);

  coef[0] = 1.0;
  coef[1] = 2.0;
  coef[2] = 1.0;
  Cbc_addRow(m, "capacity", 3, idx, coef, 'L', 5.0);

  idx[0] = 1;
  idx[1] = 2;
  coef[0] = 1.0;
  coef[1] = 1.0;
  Cbc_addRow(m, "balance", 2, idx, coef, 'E', 1.0);

  Cbc_setDblParam(m, DBL_PARAM_PRIMAL_TOL, 1.0e-5);
  Cbc_setDblParam(m, DBL_PARAM_INT_TOL, 1.0e-5);
  return m;
}

static void configure_root_only(Cbc_Model *m, const char *no_rel_repair)
{
  Cbc_setLogLevel(m, 0);
  Cbc_setMaximumNodes(m, 0);
  Cbc_setMaximumSeconds(m, 5.0);
  Cbc_setParameter(m, "preprocess", "off");
  Cbc_setParameter(m, "heuristics", "off");
  Cbc_setParameter(m, "feasibilityJump", "off");
  Cbc_setParameter(m, "noRelRepair", no_rel_repair);
}

static void configure_root_only_with_preprocess(Cbc_Model *m, const char *no_rel_repair)
{
  Cbc_setLogLevel(m, 0);
  Cbc_setMaximumNodes(m, 0);
  Cbc_setMaximumSeconds(m, 5.0);
  Cbc_setParameter(m, "preprocess", "on");
  Cbc_setParameter(m, "heuristics", "off");
  Cbc_setParameter(m, "feasibilityJump", "off");
  Cbc_setParameter(m, "noRelRepair", no_rel_repair);
}

static int expect_feasibility(Cbc_Model *m, const double *solution, int expected, const char *tag)
{
  double maxViolRow = 0.0;
  double maxViolCol = 0.0;
  int rowIdx = -1;
  int colIdx = -1;
  const int feasible = Cbc_checkFeasibility(m, solution, &maxViolRow, &rowIdx, &maxViolCol, &colIdx);
  if (feasible != expected) {
    printf("  FAIL feasibility_detector:%s expected=%d got=%d row=%g(%d) col=%g(%d)\n",
      tag, expected, feasible, maxViolRow, rowIdx, maxViolCol, colIdx);
    return 0;
  }
  return 1;
}

static int require_feasible_pool(Cbc_Model *m, const char *tag)
{
  const int nsol = Cbc_numberSavedSolutions(m);
  if (nsol <= 0) {
    printf("  FAIL %s: expected a saved solution\n", tag);
    return 0;
  }

  const double *best = Cbc_getColSolution(m);
  double maxViolRow = 0.0;
  double maxViolCol = 0.0;
  int rowIdx = -1;
  int colIdx = -1;
  if (!Cbc_checkFeasibility(m, best, &maxViolRow, &rowIdx, &maxViolCol, &colIdx)) {
    printf("  FAIL %s: incumbent infeasible row=%g(%d) col=%g(%d)\n",
      tag, maxViolRow, rowIdx, maxViolCol, colIdx);
    return 0;
  }

  if (validate_all_saved_solutions(m, 0.0, 1.0e-7, tag))
    return 0;

  if (Cbc_getNodeCount(m) != 0) {
    printf("  FAIL %s: expected root-only solve, got %d nodes\n", tag, Cbc_getNodeCount(m));
    return 0;
  }

  return 1;
}

static int test_control_without_norel(void)
{
  Cbc_Model *m = build_pair_blocks(1, 1);
  configure_root_only(m, "off");
  Cbc_solve(m);

  const int nsol = Cbc_numberSavedSolutions(m);
  const int ok = nsol == 0 && Cbc_getNodeCount(m) == 0;
  if (!ok)
    printf("  FAIL control_without_norel: nsol=%d nodes=%d\n", nsol, Cbc_getNodeCount(m));
  else
    printf("  PASS control_without_norel\n");

  Cbc_deleteModel(m);
  return ok;
}

static int test_root_heuristic_finds_solution(void)
{
  Cbc_Model *m = build_pair_blocks(1, 1);
  configure_root_only(m, "on");
  Cbc_solve(m);

  const int ok = require_feasible_pool(m, "root_heuristic");
  if (ok)
    printf("  PASS root_heuristic\n");

  Cbc_deleteModel(m);
  return ok;
}

static int test_nonflat_self_gate(void)
{
  Cbc_Model *m = build_pair_blocks(1, 0);
  configure_root_only(m, "on");
  Cbc_solve(m);

  const int nsol = Cbc_numberSavedSolutions(m);
  const int ok = nsol == 0 && Cbc_getNodeCount(m) == 0;
  if (!ok)
    printf("  FAIL nonflat_self_gate: nsol=%d nodes=%d\n", nsol, Cbc_getNodeCount(m));
  else
    printf("  PASS nonflat_self_gate\n");

  Cbc_deleteModel(m);
  return ok;
}

static int test_batched_speed(void)
{
  Cbc_Model *m = build_pair_blocks(96, 1);
  configure_root_only(m, "on");

  const double t0 = perf_wall_time();
  Cbc_solve(m);
  const double elapsed = perf_wall_time() - t0;

  int ok = require_feasible_pool(m, "batched_speed");
  if (elapsed > 2.0) {
    printf("  FAIL batched_speed: elapsed %.3fs exceeds 2.000s\n", elapsed);
    ok = 0;
  }

  if (ok)
    printf("  PASS batched_speed %.3fs\n", elapsed);

  Cbc_deleteModel(m);
  return ok;
}

static int test_feasibility_detector_cases(void)
{
  Cbc_Model *m = build_feasibility_probe();
  int ok = 1;

  const double valid[] = { 1.0, 1.0, 0.0 };
  const double rowLowerViolation[] = { 0.5, 1.0, 0.0 };
  const double rowUpperViolation[] = { 4.0, 1.0, 0.0 };
  const double colLowerViolation[] = { -0.001, 1.0, 0.0 };
  const double colUpperViolation[] = { 10.001, 1.0, 0.0 };
  const double integerViolation[] = { 1.0, 0.5, 0.5 };
  const double rowWithinTolerance[] = { 1.0 - 5.0e-6, 1.0, 0.0 };
  const double rowOutsideTolerance[] = { 1.0 - 5.0e-4, 1.0, 0.0 };
  const double integerWithinTolerance[] = { 1.0, 1.0 + 5.0e-6, 0.0 };
  const double integerOutsideTolerance[] = { 1.0, 1.0 + 5.0e-4, 0.0 };

  ok &= expect_feasibility(m, valid, 1, "valid");
  ok &= expect_feasibility(m, rowLowerViolation, 0, "row_lower");
  ok &= expect_feasibility(m, rowUpperViolation, 0, "row_upper");
  ok &= expect_feasibility(m, colLowerViolation, 0, "col_lower");
  ok &= expect_feasibility(m, colUpperViolation, 0, "col_upper");
  ok &= expect_feasibility(m, integerViolation, 0, "integer");
  ok &= expect_feasibility(m, rowWithinTolerance, 1, "row_within_tolerance");
  ok &= expect_feasibility(m, rowOutsideTolerance, 0, "row_outside_tolerance");
  ok &= expect_feasibility(m, integerWithinTolerance, 1, "integer_within_tolerance");
  ok &= expect_feasibility(m, integerOutsideTolerance, 0, "integer_outside_tolerance");

  if (ok)
    printf("  PASS feasibility_detector_cases\n");

  Cbc_deleteModel(m);
  return ok;
}

static int test_preprocess_postprocess_feasible(void)
{
  Cbc_Model *m = build_pair_blocks(8, 1);
  configure_root_only_with_preprocess(m, "on");
  Cbc_solve(m);

  const int ok = require_feasible_pool(m, "preprocess_postprocess");
  if (ok)
    printf("  PASS preprocess_postprocess\n");

  Cbc_deleteModel(m);
  return ok;
}

int main(void)
{
  int pass = 0;
  int fail = 0;

  if (test_control_without_norel())
    ++pass;
  else
    ++fail;

  if (test_root_heuristic_finds_solution())
    ++pass;
  else
    ++fail;

  if (test_nonflat_self_gate())
    ++pass;
  else
    ++fail;

  if (test_batched_speed())
    ++pass;
  else
    ++fail;

  if (test_feasibility_detector_cases())
    ++pass;
  else
    ++fail;

  if (test_preprocess_postprocess_feasible())
    ++pass;
  else
    ++fail;

  printf(fail == 0 ? "=== All %d NoRelRepair tests PASSED ===\n" : "=== %d passed, %d FAILED ===\n",
    pass, fail);
  return fail == 0 ? 0 : 1;
}
