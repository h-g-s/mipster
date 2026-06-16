/* CInterfaceTest_trdtaunimep.c
 *
 * Regression test for the trdtaunimep instance (MIPLIB 2017).
 *
 * This instance was used to detect a bug in slack_singleton_action where
 * integer singleton columns with coefficient -1 in equality rows were
 * incorrectly substituted out, dropping the col_j = col_i equality from
 * the presolved model and causing the solver to claim a wrong optimal.
 *
 * Known values:
 *   LP relaxation:  2,132,019.xx  (fractional)
 *   MIP optimal:    2,780,434.839 (certified)
 *
 * Instance structure: 23268 rows, 14741 cols
 *   13448 binary, 1293 general integer, 0 continuous
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Cbc_C_Interface.h"

#include "test_utils.h"
#define TRDTAUNIMEP_MPS fixture_path("trdtaunimep.mps.gz")

#define MIP_OPT     2780434.839
#define MIP_TOL     1.0

static int tests_run    = 0;
static int tests_passed = 0;

#define CHECK(cond, msg)                                                \
  do {                                                                  \
    ++tests_run;                                                        \
    if (cond) {                                                         \
      ++tests_passed;                                                   \
    } else {                                                            \
      fprintf(stderr, "FAIL [%s:%d]: %s\n", __FILE__, __LINE__, msg);  \
    }                                                                   \
  } while (0)

static Cbc_Model *load_model(void)
{
  Cbc_Model *m = Cbc_newModel();
  Cbc_setLogLevel(m, 0);
  int rc = Cbc_readMps(m, TRDTAUNIMEP_MPS);
  assert(rc == 0);
  return m;
}

/* ------------------------------------------------------------------ */
static void test_model_dimensions(void)
{
  printf("test_model_dimensions\n");
  Cbc_Model *m = load_model();

  CHECK(Cbc_getNumRows(m) == 23268, "expected 23268 rows");
  CHECK(Cbc_getNumCols(m) == 14741, "expected 14741 cols");

  int nbin = 0, nint = 0, ncont = 0;
  const double *lb = Cbc_getColLower(m);
  const double *ub = Cbc_getColUpper(m);
  for (int j = 0; j < Cbc_getNumCols(m); ++j) {
    if (Cbc_isInteger(m, j)) {
      if (lb[j] == 0.0 && ub[j] == 1.0) nbin++;
      else nint++;
    } else {
      ncont++;
    }
  }
  CHECK(nbin  == 13448, "expected 13448 binary vars");
  CHECK(nint  ==  1293, "expected 1293 general integer vars");
  CHECK(ncont ==     0, "expected 0 continuous vars");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
static void test_mip(void)
{
  printf("test_mip\n");
  Cbc_Model *m = load_model();

  Cbc_setIntParam(m, INT_PARAM_MAX_NODES, 300);
  Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, 300.0);
  Cbc_setDblParam(m, DBL_PARAM_GAP_RATIO, 0.0);

  Cbc_solve(m);

  int is_optimal = Cbc_isProvenOptimal(m);

  /* Feasibility of all solutions in the pool is always required */
  int fails = validate_all_saved_solutions(m, MIP_OPT, MIP_TOL, "trdtaunimep");
  CHECK(fails == 0, "all saved solutions feasible and objective on correct side");

  /* If solver claims optimality, verify the objective matches the certified value */
  if (is_optimal) {
    double obj = Cbc_getObjValue(m);
    CHECK(fabs(obj - MIP_OPT) < MIP_TOL,
          "claimed optimal obj must match certified 2,780,434.839");
  } else {
    /* Not proven optimal within 300 nodes — acceptable as long as solutions are feasible */
    printf("  note: not proven optimal within 300 nodes (feasibility checked above)\n");
    ++tests_run;
    ++tests_passed; /* not a failure */
  }

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
int main(void)
{
  printf("=== trdtaunimep C-interface tests ===\n");

  test_model_dimensions();
  test_mip();

  printf("=== %d / %d tests passed ===\n", tests_passed, tests_run);
  return (tests_passed == tests_run) ? 0 : 1;
}
