// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).


/*! \file CbcLinkedUtils.cpp
    \brief Stub for the (now removed) AMPL/ASL hooks used by CbcLinked.

    MIPster no longer links against the AMPL Solver Library (ASL), so the
    nonlinear-objective loading path is unavailable. ClpSimplex_loadNonLinear
    is kept as a stub so callers that still reference it continue to link;
    it aborts if ever invoked.
*/

#include "ClpConfig.h"
#include "ClpSimplex.hpp"
#include "ClpConstraint.hpp"

int ClpSimplex_loadNonLinear(ClpSimplex &cs, void *, int &,
  ClpConstraint **&)
{
  abort();
  return 0;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
