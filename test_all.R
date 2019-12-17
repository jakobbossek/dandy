library(methods)
library(devtools)
library(testthat)

if (interactive()) {
  load_all(".")
} else {
  library(sampling)
}

test_dir("tests/testthat")
