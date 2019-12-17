#include <R.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>

#include "dem_discr.h"
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
