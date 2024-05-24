
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bsvarSIGNs

Developing an R package for Bayesian Structural VARs identified by zero, sign,
and narrative restrictions

# Installation

``` r
devtools::install_github("bsvars/bsvarSIGNs")
```

# Example

A repliaction of Antolín-Díaz and Rubio-Ramírez (2018).

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

A repliaction of Arias, Rubio-Ramírez and Waggoner (2018).

``` r
data(optimism)

# try ?specify_bsvarSIGN
zero_irf          = matrix(0, nrow = 5, ncol = 5)
zero_irf[1, 1]    = 1
sign_irf          = array(0, dim = c(5, 5, 1))
sign_irf[2, 1, 1] = 1

specification  = specify_bsvarSIGN$new(optimism*100,
                                       p = 4,
                                       sign_irf = sign_irf,
                                       zero_irf = zero_irf
                                       )

posterior      = estimate(specification, S = 10000)
irf            = compute_impulse_responses(posterior, horizon = 40)
plot(irf, probability = 0.68)
```

<!-- badges: start -->

[![R-CMD-check](https://github.com/bsvars/bsvarSIGNs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bsvars/bsvarSIGNs/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
