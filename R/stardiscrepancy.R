#' @title
#' Calculate star discrepancy of a set of points.
#'
#' @description
#' To be written ...
#'
#' @param x [\code{matrix(n, d)}]\cr
#'   An \eqn{n \times d} matrix where \eqn{n} is the number of points and \eqn{d}
#'   is the dimension of the search space.
#' @return [\numeric(1)]\cr Star discrepancy of \code{x}.
#' @export
stardiscrepancy = function(x) {
  if (checkmate::testDataFrame(x))
    x = unname(as.matrix(x))

  checkmate::assertMatrix(x, min.rows = 2L, min.cols = 2L, any.missing = FALSE, all.missing = FALSE, mode = "numeric")
  return(.Call("starDiscrepancyC", t(x)))
}
