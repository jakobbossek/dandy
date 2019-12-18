
# Sampling and (Approximate) Star Discrepancy Calculation

[![Status](https://img.shields.io/badge/Status-experimental-red.svg)](https://GitHub.com/jakobbossek/findAName)
[![CRAN Status
Badge](http://www.r-pkg.org/badges/version/findAName)](http://cran.r-project.org/web/packages/findAName)
[![CRAN
checks](https://cranchecks.info/badges/worst/findAName)](https://cran.r-project.org/web/checks/check_results_findAName.html)
[![CRAN
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/findAName?color=orange)](http://cran.rstudio.com/web/packages/findAName/index.html)

[![Build
Status](https://travis-ci.org/jakobbossek/findAName.svg?branch=master)](https://travis-ci.org/jakobbossek/findAName)
[![Build
status](https://ci.appveyor.com/api/projects/status/eu0nns2dsgocwntw/branch/master?svg=true)](https://ci.appveyor.com/project/jakobbossek/findAName/branch/master)
[![Coverage
Status](https://coveralls.io/repos/github/jakobbossek/findAName/badge.svg?branch=master)](https://coveralls.io/github/jakobbossek/findAName?branch=master)

## Introduction

## Example

Here, we generate each one uniform, latin hypercube sample, Halton and
Sobol design with n=20 points in k=2 dimensions and calculate its exact
disrepancy.

``` r
library(findAName)
library(ggplot2)

n = 50
methods = c("uniform", "improvedlhs", "halton", "sobol")
designs = lapply(methods, function(method) {
  findAName::design(n = n, k = 5, method = method)
})

# caluclate star discrepancy
discrepancies = sapply(designs, findAName::stardiscrepancy, method = "exact")
discrepanciesta = sapply(designs, findAName::stardiscrepancy, method = "ta")

# visualize
designs = do.call(rbind, designs)
groups = sprintf("%s (%.4f)", methods, discrepancies)
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

## How to contribute?

You can contribute by identifing annoying bugs in the [issue
tracker](http://github.com/jakobbossek/findAName). This is also the
preferred place to ask questions and raise feature requests. Moreover,
users can contribute even more by
[forking](https://help.github.com/en/github/getting-started-with-github/fork-a-repo)
the findAName repository, implementing features or bugfixes and raising
a [pull
request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests).

## Installation Instructions

The package will be available at [CRAN](http://cran.r-project.org) *when
it is done*. If you are interested in trying out and playing around with
the current github developer version use the
[devtools](https://github.com/hadley/devtools) package and type the
following command in R:

``` r
remotes::install_github("jakobbossek/findAName")
```

## Getting help

Please address questions and missing features about the *findAName* as
weell as annoying bug reports in the [issue
tracker](https://github.com/jakobbossek/findAName/issues). Pay attention
to explain your problem as good as possible. At its best you provide an
example, so I can reproduce your problem quickly. Please avoid sending
e-mails.
