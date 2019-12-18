#include <R.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>

#include "dem_discr.h"
#include "macros.h"

SEXP starDiscrepancyC(SEXP r_points);
SEXP starDiscrepancyTAC(SEXP r_points, SEXP r_iter, SEXP r_max_trials);
