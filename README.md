
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bsvarSIGNs

Developing an R package for Bayesian Structural VARs identified by zero,
sign and narrative restrictions.

# Installation

``` r
devtools::install_github("bsvars/bsvarSIGNs")
```

# Example

A repliaction of Antolín-Díaz and Rubio-Ramírez (2018).

``` r
data(monetary)

# contractionary monetary policy shock
sign_irf       = matrix(0, 6, 6)
sign_irf[, 1]  = c(0, -1, -1, 0, -1, 1)
sign_irf       = array(sign_irf, dim = c(6, 6, 5))

# in October 1979 the shock
sign_narrative = rbind(c(1, 1, NA, 1, 166, 0),   # is positive
                       c(3, 1, 6,  1, 166, 0))   # greatest historical decomposition

specification  = specify_bsvarSIGN$new(monetary       * 100,
                                       p              = 12,
                                       sign_irf       = sign_irf,
                                       sign_narrative = sign_narrative)

posterior      = estimate(specification, S = 100)
irf            = compute_impulse_responses(posterior, horizon = 60)
plot(irf, probability = 0.68)
```

A repliaction of Arias, Rubio-Ramírez and Waggoner (2018).

``` r
data(optimism)

# optimism shock
# no contemporaneous effect on productivity
zero_irf          = matrix(0, nrow = 5, ncol = 5)
zero_irf[1, 1]    = 1
# positive contemporaneous effect on stock prices
sign_irf          = array(0, dim = c(5, 5, 1))
sign_irf[2, 1, 1] = 1

specification  = specify_bsvarSIGN$new(optimism * 100,
                                       p        = 4,
                                       sign_irf = sign_irf,
                                       zero_irf = zero_irf)

posterior      = estimate(specification, S = 100)
irf            = compute_impulse_responses(posterior, horizon = 40)
plot(irf, probability = 0.68)
```

<!-- badges: start -->

[![R-CMD-check](https://github.com/bsvars/bsvarSIGNs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bsvars/bsvarSIGNs/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
