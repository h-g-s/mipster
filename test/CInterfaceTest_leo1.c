/* CInterfaceTest_leo1.c
 *
 * Node-limited smoke test for the leo1 MIPLIB 2017 instance.
 *
 * Known values:
 *   MIP optimal (MIPLIB certified): 404227536.16
 *
 * Instance structure: minimisation, mixed-integer.
 *
 * Test strategy (per testing guide):
 *   - Node limit (200) is the primary termination criterion.
 *   - Wall-clock time limit (600 s) is a loose fallback only.
 *   - At 200 nodes the instance is not yet proven optimal; we validate
 *     every solution in the pool for feasibility.
 *   - If the solver happens to claim optimality (e.g. presolve collapses
 *     the problem on some platforms), we also verify the objective matches
 *     the certified value within tolerance.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Cbc_C_Interface.h"
#include "test_utils.h"

#define LEO1_MPS fixture_path("leo1.mps.gz")

/* MIPLIB certified optimal */
#define MIP_OPT  404227536.16
/* Absolute tolerance for objective comparison */
#define MIP_TOL  1.0

#define NODE_LIMIT 200

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
  int rc = Cbc_readMps(m, LEO1_MPS);
  assert(rc == 0);
  return m;
}

/* ── Test: node-limited MIP ─────────────────────────────────────────── */
static void test_mip_nodelimited(void)
{
  printf("test_mip_nodelimited\n");
  Cbc_Model *m = load_model();

  /* Node limit is the deterministic primary stop. */
  Cbc_setMaximumNodes(m, NODE_LIMIT);
  /* Loose wall-clock fallback only — never the primary criterion. */
  Cbc_setDblParam(m, DBL_PARAM_TIME_LIMIT, 600.0);

  Cbc_solve(m);

  int is_optimal = Cbc_isProvenOptimal(m);

  /* Every solution in the pool must be feasible. */
  int fails = validate_all_saved_solutions(m,
    is_optimal ? MIP_OPT : NAN,
    MIP_TOL, "leo1");
  CHECK(fails == 0,
        "all saved solutions are feasible and consistent with certified opt");

  /* If the solver claims optimality, verify the objective. */
  if (is_optimal) {
    double obj = Cbc_getObjValue(m);
    CHECK(fabs(obj - MIP_OPT) < MIP_TOL,
          "if proven optimal, obj must match MIPLIB certified 404227536.16");
  }

  Cbc_deleteModel(m);
}

/* ── Main ────────────────────────────────────────────────────────────── */
int main(void)
{
  printf("=== leo1 C-interface tests ===\n");
  test_mip_nodelimited();
  printf("=== %d / %d tests passed ===\n", tests_passed, tests_run);
  return (tests_passed == tests_run) ? 0 : 1;
}
