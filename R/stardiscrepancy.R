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
stardiscrepancy = function(x, force.exact = FALSE, iter = 1e4, trials = 10L) {
  if (checkmate::testDataFrame(x))
    x = unname(as.matrix(x))

  checkmate::assertMatrix(x, min.rows = 2L, min.cols = 2L, any.missing = FALSE, all.missing = FALSE, mode = "numeric")
  n = nrow(x)
  k = ncol(x)
  # if (n^(1+ceiling(k/2)) > 1e2) {

  #   BBmisc::messagef("[sampling::stardiscrepancy] Exact star-discrepancy calculation requires
  #     time O(n^{1+k/2}) time.\n We thus apply heuristic TA-algorithm.")

  #   if (!force.exact)
  #     return(.Call("starDiscrepancyTAC", t(x)))
  # }

  #FIXME: seeding
  return(.Call("starDiscrepancyTAC", t(x), as.integer(iter), as.integer(trials)))

  # return(.Call("starDiscrepancyC", t(x)))
}
