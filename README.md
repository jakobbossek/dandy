
# dandy: Desings ANd DiscrepancY

[![Status](https://img.shields.io/badge/Status-experimental-red.svg)](https://GitHub.com/jakobbossek/dandy)
[![CRAN Status
Badge](http://www.r-pkg.org/badges/version/dandy)](http://cran.r-project.org/web/packages/dandy)
[![CRAN
checks](https://cranchecks.info/badges/worst/dandy)](https://cran.r-project.org/web/checks/check_results_dandy.html)
[![CRAN
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/dandy?color=orange)](http://cran.rstudio.com/web/packages/dandy/index.html)

[![Build
Status](https://travis-ci.org/jakobbossek/dandy.svg?branch=master)](https://travis-ci.org/jakobbossek/dandy)
[![Build
status](https://ci.appveyor.com/api/projects/status/eu0nns2dsgocwntw/branch/master?svg=true)](https://ci.appveyor.com/project/jakobbossek/dandy/branch/master)
[![Coverage
Status](https://coveralls.io/repos/github/jakobbossek/dandy/badge.svg?branch=master)](https://coveralls.io/github/jakobbossek/dandy?branch=master)

## Introduction

This package offers a simple method for the generation of random
uniform, LHS, Sobol and Halton sequences of n points in d dimension.
Basically the package delegates to other packages and offers a unified
interface. Additionally, an exact method and an approximate method based
on threshold accepting optimization for the calculation of the star
discrepancy measure is included.

## Example

Here, we generate each one uniform, latin hypercube sample, Halton and
Sobol design with n=20 points in k=2 dimensions and calculate its exact
disrepancy.

``` r
library(dandy)
library(ggplot2)

n = 50
methods = c("uniform", "improvedlhs", "halton", "sobol")
designs = lapply(methods, function(method) {
  dandy::design(n = n, k = 5, method = method)
})

# caluclate star discrepancy
discr.exact  = sapply(designs, dandy::stardiscrepancy, method = "exact")
discr.approx = sapply(designs, dandy::stardiscrepancy, method = "ta")

# visualize
designs = do.call(rbind, designs)
groups = sprintf("%s (%.4f)", methods, discr.exact)
designs$method = factor(rep(groups, each = n))

pl = ggplot(designs, aes(x = x1, y = x2))
pl = pl + geom_point()
pl = pl + facet_grid(. ~ method)
print(pl)
```

## Development Team

The package is a one-man project by [Jakob
Bossek](https://researchers.adelaide.edu.au/profile/jakob.bossek) at the
moment of writing. However, the package interfaces a neat implementation
of approximate star-discrepancy caluclation by [Magnus
Wahlström](https://pure.royalholloway.ac.uk/portal/en/persons/magnus-wahlstroem\(a0940b3f-c15b-404f-b9f6-d26cb5664829\).html).

Gnewuch, Michael, Magnus Wahlström, and Carola Winzen. “A NEW RANDOMIZED
ALGORITHM TO APPROXIMATE THE STAR DISCREPANCY BASED ON THRESHOLD
ACCEPTING.” SIAM Journal on Numerical Analysis 50, no. 2 (2012):
781-807. www.jstor.org/stable/41582760.

## How to contribute?

You can contribute by identifing annoying bugs in the [issue
tracker](http://github.com/jakobbossek/dandy). This is also the
preferred place to ask questions and raise feature requests. Moreover,
users can contribute even more by
[forking](https://help.github.com/en/github/getting-started-with-github/fork-a-repo)
the dandy repository, implementing features or bugfixes and raising a
[pull
request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests).

## Installation Instructions

The package will be available at [CRAN](http://cran.r-project.org) *when
it is done*. If you are interested in trying out and playing around with
the current github developer version use the
[devtools](https://github.com/hadley/devtools) package and type the
following command in R:

``` r
remotes::install_github("jakobbossek/dandy")
```

## Getting help

Please address questions and missing features about the *dandy* as weell
as annoying bug reports in the [issue
tracker](https://github.com/jakobbossek/dandy/issues). Pay attention to
explain your problem as good as possible. At its best you provide an
example, so I can reproduce your problem quickly. Please avoid sending
e-mails.
