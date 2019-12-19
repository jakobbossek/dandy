context("discrepancy")

test_that("discrepancy calculation works", {
  dims = 2:4
  n = 20L

  for (d in dims) {
    des = design(n = n, k = d, method = "uniform", as.df = TRUE)
    checkmate::expect_number(discrepancy(des, method = "exact"), lower = 0, upper = 1)
    checkmate::expect_number(discrepancy(des, method = "ta", iters = 100L, trials = 1L), lower = 0, upper = 1)
  }
})
