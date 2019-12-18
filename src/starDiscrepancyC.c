#include <R.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>

#include "dem_discr.h"
#include "TA_shrink_delta.h"
#include "TA_shrink_bardelta.h"
#include "starDiscrepancyC.h"
#include "macros.h"

SEXP starDiscrepancyC(SEXP r_points) {
  // first unpack R structures
  EXTRACT_NUMERIC_MATRIX(r_points, c_points, dim, n_points);

  double **pointset, upper, lower;

  // ugly way to convert SEXP to double array
  pointset = malloc(n_points * sizeof(double*));
  for (int i = 0; i < n_points; i++) {
    pointset[i] = malloc(dim * sizeof(double));

    // get start
    double *start = c_points + i * dim;

    for (int j = 0; j < dim; j++) {
      // newline counts as whitespace
      pointset[i][j] = start[j];
    }
  }

  upper = oydiscr(pointset, dim, n_points, &lower);

  return ScalarReal(upper);
}

SEXP starDiscrepancyTAC(SEXP r_points, SEXP r_iter, SEXP r_max_trials) {
  // first unpack R structures
  EXTRACT_NUMERIC_MATRIX(r_points, c_points, dim, n_points);
  EXTRACT_INTEGER(r_iter, iter);
  EXTRACT_INTEGER(r_max_trials, max_trials);

  double **pointset, delta_result, bardelta_result;

  // ugly way to convert SEXP to double array
  pointset = malloc(n_points * sizeof(double*));
  for (int i = 0; i < n_points; i++) {
    pointset[i] = malloc(dim * sizeof(double));

    // get start
    double *start = c_points + i * dim;

    for (int j = 0; j < dim; j++) {
      // newline counts as whitespace
      pointset[i][j] = start[j];
    }
  }

  //delta_result = delta_calc(pointset, n_points, dim);
  bardelta_result = bardelta_calc(pointset, n_points, dim);
  return ScalarReal(bardelta_result);

  //return ScalarReal(MIN(delta_result, bardelta_result));
}
