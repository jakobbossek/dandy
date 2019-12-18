context("design generation")

test_that("designs are generated correctly", {
  methods = getSupportedMethods()
  dims = 2:4
  n = 20L

  # check basic functionality
  for (method in methods) {
    for (d in dims) {
      des = design(n = n, k = d, method = method, as.df = TRUE)
      checkmate::expect_data_frame(des, types = "numeric", nrows = n, ncols = d,
        any.missing = FALSE, all.missing = FALSE)
    }
  }

  # check if passing of further arguments works
  # Here we deactivate scrambling and check whether it works
  d1 = design(n = n, k = 2L, method = "sobol", as.df = FALSE, scrambling = 0)
  d2 = design(n = n, k = 2L, method = "sobol", as.df = FALSE, scrambling = 0)
  testthat::expect_true(all(d1 == d2))
})
