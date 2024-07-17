
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bsvarSIGNs

An **R** package for Bayesian Estimation of Structural Vector
Autoregressions Identified by Sign, Zero, and Narrative Restrictions

<!-- badges: start -->

[![R-CMD-check](https://github.com/bsvars/bsvarSIGNs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bsvars/bsvarSIGNs/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Implements state-of-the-art algorithms for the Bayesian analysis of
Structural Vector Autoregressions identified by sign, zero, and
narrative restrictions. The core model is based on the flexible Vector
Autoregression with the estimated hyper-parameters of the Minnesota
prior as in [Giannone, Lenza, Primiceri
(2015)](http://doi.org/10.1162/REST_a_00483). The sign restrictions are
implemented employing the methods outlined by [Rubio-Ramírez, Waggoner &
Zha (2010)](http://doi.org/10.1111/j.1467-937X.2009.00578.x), while
identification through sign and zero restrictions follows the approach
developed by [Arias, Rubio-Ramírez, & Waggoner
(2018)](http://doi.org/10.3982/ECTA14468). Furthermore, our tool
provides algorithms for identification via sign and narrative
restrictions, in line with the methods introduced by [Antolín-Díaz and
Rubio-Ramírez (2018)](http://doi.org/10.1257/aer.20161852). Users can
also estimate a model with sign, zero, and narrative restrictions
imposed at once. The package facilitates predictive and structural
analyses using impulse responses, forecast error variance and historical
decompositions, forecasting and conditional forecasting, as well as
analyses of structural shocks and fitted values. All this is
complemented by colourful plots, user-friendly summary functions, and
comprehensive documentation. The **bsvarSIGNs** package is aligned
regarding code structure, objects, and workflows with the **R** package
**bsvars** by [Woźniak
(2024)](http://doi.org/10.32614/CRAN.package.bsvars), and they
constitute an integrated toolset.

## Installation

#### The first time you install the package

You must have a **cpp** compiler. Follow the instructions from [Section
1.3. by Eddelbuettel & François
(2023)](https://cran.r-project.org/package=Rcpp/vignettes/Rcpp-FAQ.pdf).
In short, for **Windows:** install
[RTools](https://CRAN.R-project.org/bin/windows/Rtools/), for **macOS:**
install [Xcode Command Line
Tools](https://www.freecodecamp.org/news/install-xcode-command-line-tools/),
and for **Linux:** install the standard developement packages.

#### Once that’s done:

Just open your **R** and type:

    install.packages("bsvarSIGNs")

The developer’s version of the package with the newest features can be
installed by typing:

    devtools::install_github("bsvars/bsvarSIGNs")

## Development

The package is under intensive development. Your help is most welcome!
Please, have a look at the
[roadmap](https://github.com/bsvars/bsvarSIGNs/milestones) and [report a
bug](https://github.com/bsvars/bsvarSIGNs/issues). Thank you!

## Example

A replication of [Arias, Rubio-Ramírez, & Waggoner
(2018)](http://doi.org/10.3982/ECTA14468).

``` r
data(optimism)

# optimism shock
# no effect on productivity + positive effect on stock prices
sign_irf       = matrix(c(0, 1, rep(NA, 23)), 5, 5)

specification  = specify_bsvarSIGN$new(optimism * 100,
                                       p        = 4,
                                       sign_irf = sign_irf)

posterior      = estimate(specification, S = 100)
irf            = compute_impulse_responses(posterior, horizon = 40)
plot(irf, probability = 0.68)
```

A replication of [Antolín-Díaz and Rubio-Ramírez
(2018)](http://doi.org/10.1257/aer.20161852).

``` r
data(monetary)

# contractionary monetary policy shock
sign_irf       = matrix(NA, 6, 6)
sign_irf[, 1]  = c(NA, -1, -1, NA, -1, 1)
sign_irf       = array(sign_irf, dim = c(6, 6, 5))

# the shock is positive in October 1979
sign_narrative = t(c(1, 1, NA, 1, 166, 0))

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

specification     = specify_bsvarSIGN$new(optimism * 100,
                                          p        = 4,
                                          sign_irf = sign_irf,
                                          zero_irf = zero_irf)

posterior         = estimate(specification, S = 100)
irf               = compute_impulse_responses(posterior, horizon = 40)
plot(irf, probability = 0.68)
```
