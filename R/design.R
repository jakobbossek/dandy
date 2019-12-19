#' @title Design generator.
#'
#' @description Generate pseudo-random (uniform, latin-hypercube samle) or quasi-random
#' (Halton and Sobol sequences) designs. Basically a wrapper
#' around different functions from packages \pkg{lhs} and \pkg{randtoolbox}.
#'
#' @param n [\code{integer(1)}]\cr
#'   Number of design points (rows).
#' @param k [\code{integer(1)}]\cr
#'   Number of variables (columns).
#' @param l [\code{numeric}]\cr
#'   Lower bound for design points. Either a single numeric value
#'   or a vector of length \code{k}.
#'   Default is 0.
#' @param u [\code{numeric}]\cr
#'   Upper bound for design points. Either a single numeric value
#'   or a vector of length \code{k}.
#'   Default is 1.
#' @param method [\code{character(1)}]\cr
#'   Name of method used to generate the design. Possible values are
#'   \describe{
#'     \item{uniform}{Uniform random sampling.}
#'     \item{improvedlhs}{Delegate to \code{\link[lhs]{improvedLHS}}.}
#'     \item{maximinlhs}{Delegate to \code{\link[lhs]{maximinLHS}}.}
#'     \item{geneticlhs}{Delegate to \code{\link[lhs]{geneticLHS}}.}
#'     \item{halton}{Delegate to \code{\link[randtoolbox]{halton}}.}
#'     \item{sobol}{Delegate to \code{\link[randtoolbox]{sobol}} with option scrambling=3.}
#'   }
#' @param as.df [\code{logical(1)}]\cr
#'   Return points as data frame?
#'   Default is \code{TRUE}.
#' @param ... [any]\cr
#'   Further parameters passed down to generator if applicable.
#' @return Design points as a matrix or data.frame (see parameter \code{as.df}).
#' @examples
#' methods = getSupportedMethods()
#' designs = lapply(methods, function(method) {
#'   design(n = 100, k = 2, method = method, l = -5, u = c(5, 10), as.df = FALSE)
#' })
#'
#' # pass down options to generator
#' d = design(n = 50, k = 4, method = "sobol", scrambling = 2, seed = 123)
#' print(d)
#' @export
design = function(n, k, method, l = 0, u = 1, as.df = TRUE, ...) {
  n = checkmate::asInt(n, lower = 2L, na.ok = FALSE)
  k = checkmate::asInt(k, lower = 1L, na.ok = FALSE)
  checkmate::assertChoice(method, choices = getSupportedMethods())
  checkmate::assertFlag(as.df)

  if (length(l) == 1L)
    l = rep(l, k)
  if (length(u) == 1L)
    u = rep(u, k)

  checkmate::assertNumeric(l, len = k, any.missing = FALSE, all.missing = FALSE)
  checkmate::assertNumeric(u, len = k, any.missing = FALSE, all.missing = FALSE)

  if (any(l >= u))
    BBmisc::stopf("[dandy::design] Lower bounds must be strictly lower than upper bounds.")

  des = if (method == "uniform") {
    matrix(runif(n * k), nrow = n)
  } else if (method == "improvedlhs") {
    lhs::improvedLHS(n = n, k = k, ...)
  } else if (method == "maximinlhs") {
    lhs::maximinLHS(n = n, k = k, ...)
  } else if (method == "geneticlhs") {
    lhs::geneticLHS(n = n, k = k, ...)
  } else if (method == "halton") {
    randtoolbox::halton(n = n, dim = k, ...)
  } else if (method == "sobol") {
    args = BBmisc::insert(list(scrambling = 3, seed = ceiling(runif(1L, min = 1, max = 10000000))), list(...))
    args$n = n; args$dim = k
    do.call(randtoolbox::sobol, args)
  } else {
    BBmisc::stopf("[dandy::design] Unsupported method '%s'.", method)
  }

  # linear transformation
  des = t(t(des) * (u - l) + l)

  if (as.df) {
    des = as.data.frame(des)
    colnames(des) = paste0("x", seq_len(k))
  }

  return(des)
}

#' Get character vector of supported design generators,
#' i.e., possible value for parameter \code{method} of
#' function \code{\link{design}}.
#' @export
getSupportedMethods = function() {
  c("uniform", "improvedlhs", "maximinlhs", "geneticlhs", "halton", "sobol")
}
