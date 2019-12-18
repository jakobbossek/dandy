#' @title
#' Star discrepancy calculation.
#'
#' @description
#' Method for calculating the (star) discrepancy of \eqn{n} points in \eqn{d}
#' dimensions. The function offers an exact approach with running time
#' \eqn{O(n^{1+d/2})} and a sophisticated approximation method based on
#' threshold accepting introduced by Gnewuch, Wahlström, and Winzen [1].
#'
#' @param x [\code{matrix(n, d)}]\cr
#'   An \eqn{n \times d} matrix where \eqn{n} is the number of points and \eqn{d}
#'   is the dimension of the search space.
#' @param method [\code{character(1)}]\cr
#'   Use option \dQuote{exact} for exact star discrepancy calculation and
#'   \dQuote{ta} for the threshold accepting based discrepancy algorithm by
#'   Gnewuch, Wahlstroem and Winzen [1].
#' @param iters [\code{integer(1)}]\cr
#'   Number of iterations for the treshold accepting discrepancy approximation.
#'   Default is 100000.
#' @param trials [\code{integer(1)}]\cr
#'   Number of independent trials of threshold accepting discrepancy approximation.
#'   Default is 10.
#' @return [\code{numeric}(1)]\cr Star discrepancy of \code{x}.
#'
#' @references [1] Gnewuch, Michael, Magnus Wahlström, and Carola Winzen. "A NEW RANDOMIZED
#' ALGORITHM TO APPROXIMATE THE STAR DISCREPANCY BASED ON THRESHOLD ACCEPTING."
#' SIAM Journal on Numerical Analysis 50, no. 2 (2012): 781-807.
#' www.jstor.org/stable/41582760.
#'
#' @examples
#' d = design(n = 20, k = 3, method = "uniform")
#' stardiscrepancy(d)
#' \dontrun{
#' stardiscrepancy(d, method = "ta", iter = 100, trials = 3)
#' }
#' @export
stardiscrepancy = function(x, method = "exact", iters = 1e5, trials = 10) {
  if (checkmate::testDataFrame(x))
    x = unname(as.matrix(x))

  checkmate::assertMatrix(x, min.rows = 2L, min.cols = 2L, any.missing = FALSE, all.missing = FALSE, mode = "numeric")
  n = nrow(x)
  d = ncol(x)

  checkmate::assertChoice(method, choices = c("exact", "ta"))
  if (method == "exact") {
    if (n^(1 + ceiling(d/2)) > 1e6) {
      BBmisc::messagef("[sampling::stardiscrepancy] Exact star-discrepancy calculation requires
        time O(n^{1+d/2}) time.\n Go grab yourself a coffee. This may take some time.")
    }
    return(.Call("starDiscrepancyC", t(x[])))
  }

  #SEEDING
  return(.Call("starDiscrepancyTAC", t(x[]), as.integer(iters), as.integer(trials)))
}
