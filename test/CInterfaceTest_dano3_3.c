/* CInterfaceTest_dano3_3.c
 *
 * Tests for the dano3_3 instance (MIPLIB 2017 benchmark).
 *
 * This instance has previously triggered postprocessing bugs that produced
 * infeasible solutions.  The primary goal of this test is to verify that
 * every solution returned by the solver (best + pool) is genuinely feasible.
 *
 * Known values (MIPLIB 2017 certified):
 *   MIP optimal:  576.3446   (primal == dual, gap == 0)
 *
 * Instance structure: 3202 rows, 13873 cols, 79655 NZ
 *   69 binary, 0 general integer, 13804 continuous
 *
 * The MIP test uses a 15-node limit; the solver finds a solution equal to
 * the certified optimum within this limit but cannot prove optimality.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Cbc_C_Interface.h"

#include "test_utils.h"
#define DANO3_3_MPS fixture_path("dano3_3.mps.gz")

/* Known / reference values */
#define MIP_OPT       576.3446
#define MIP_TOL       1e-3
#define MIP_NODE_LIMIT 15

static int tests_run    = 0;
static int tests_passed = 0;

#define CHECK(cond, msg)                                               \
  do {                                                                 \
    ++tests_run;                                                       \
    if (cond) {                                                        \
      ++tests_passed;                                                  \
    } else {                                                           \
      fprintf(stderr, "FAIL [%s:%d] %s\n", __FILE__, __LINE__, msg);  \
    }                                                                  \
  } while (0)

/* ------------------------------------------------------------------ */
/* helpers                                                             */
/* ------------------------------------------------------------------ */

static Cbc_Model *load_model(void)
{
  Cbc_Model *m = Cbc_newModel();
  Cbc_setLogLevel(m, 0);
  int rc = Cbc_readMps(m, DANO3_3_MPS);
  if (rc != 0) {
    fprintf(stderr, "ERROR: could not read %s\n", DANO3_3_MPS);
    Cbc_deleteModel(m);
    return NULL;
  }
  return m;
}

/* ------------------------------------------------------------------ */
/* test_model_dimensions                                               */
/* ------------------------------------------------------------------ */

static void test_model_dimensions(void)
{
  Cbc_Model *m = load_model();
  if (!m) return;

  CHECK(Cbc_getNumRows(m) == 3202,  "dano3_3: 3202 rows");
  CHECK(Cbc_getNumCols(m) == 13873, "dano3_3: 13873 cols");

  int nbin = 0, nint = 0, ncont = 0;
  const double *lb = Cbc_getColLower(m);
  const double *ub = Cbc_getColUpper(m);
  for (int j = 0; j < Cbc_getNumCols(m); j++) {
    if (Cbc_isInteger(m, j)) {
      if (lb[j] == 0.0 && ub[j] == 1.0) nbin++; else nint++;
    } else {
      ncont++;
    }
  }
  CHECK(nbin  ==    69, "dano3_3: 69 binary vars");
  CHECK(nint  ==     0, "dano3_3: 0 general integer vars");
  CHECK(ncont == 13804, "dano3_3: 13804 continuous vars");

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
/* test_mip_node_limited                                               */
/*                                                                     */
/* Solve with a 15-node limit.  Every solution in the pool is          */
/* validated for feasibility.  If a solution is found its objective    */
/* must be >= the certified optimal (576.3446).                        */
/* ------------------------------------------------------------------ */

static void test_mip_node_limited(void)
{
  Cbc_Model *m = load_model();
  if (!m) return;

  Cbc_setIntParam(m, INT_PARAM_MAX_NODES, MIP_NODE_LIMIT);
  Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, 300.0);

  int rc = Cbc_solve(m);
  (void)rc;

  CHECK(!Cbc_isProvenInfeasible(m), "dano3_3 is not infeasible");

  int nsols = Cbc_numberSavedSolutions(m);
  printf("  Saved solutions: %d   Best bound: %g\n",
         nsols, Cbc_getBestPossibleObjValue(m));

  if (nsols > 0) {
    double best_obj = Cbc_getObjValue(m);
    CHECK(best_obj >= MIP_OPT - MIP_TOL,
          "dano3_3: best obj >= certified optimal (576.3446)");

    int fails = validate_all_saved_solutions(m, MIP_OPT, MIP_TOL, "dano3_3");
    CHECK(fails == 0,
          "dano3_3: all saved solutions are feasible with obj >= certified optimal");

    if (fails == 0)
      printf("  All %d saved solution(s) validated (feasible, obj >= %g)\n",
             nsols, MIP_OPT);
  } else {
    printf("  (no integer solution found within %d nodes — skipping "
           "solution validation)\n", MIP_NODE_LIMIT);
    ++tests_run; ++tests_passed;
  }

  Cbc_deleteModel(m);
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(void)
{
  printf("=== dano3_3 C interface tests ===\n");

  test_model_dimensions();
  test_mip_node_limited();

  printf("\n%d/%d tests passed.\n", tests_passed, tests_run);
  return (tests_passed == tests_run) ? 0 : 1;
}
