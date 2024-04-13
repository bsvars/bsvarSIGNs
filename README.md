
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bsvarSIGNs

Developing an R package for Bayesian Structural VARs identified by sign
and narrative restrictions

# Installation

``` r
devtools::install_github("bsvars/bsvarSIGNs")
```

# Example

A repliaction of Antolín-Díaz & Rubio-Ramírez (2018).

``` r
data(oil)

# try ?specify_bsvarSIGN
sign_narrative = matrix(c(2, -1, 3, 2, 236, 0), ncol = 6)
sign_irf       = array(matrix(c(-1, -1, 1, 1, 1, 1, 1, -1, 1), nrow = 3), dim = c(3, 3, 1))

specification  = specify_bsvarSIGN$new(oil,
                                        p = 12,
                                        sign_irf = sign_irf,
                                        sign_narrative = sign_narrative)

posterior      = estimate(specification, S = 100)
irf            = compute_impulse_responses(posterior, horizon = 24)
plot(irf)
```

<!-- badges: start -->

[![R-CMD-check](https://github.com/bsvars/bsvarSIGNs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bsvars/bsvarSIGNs/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
