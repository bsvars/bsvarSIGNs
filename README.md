
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bsvarSIGNs

Developing an R package for Bayesian Structural VARs identified by sign
and narrative restrictions

# Installation

``` r
devtools::install_github("bsvars/bsvarSIGNs")
#> Downloading GitHub repo bsvars/bsvarSIGNs@HEAD
#> RcppArmad... (0.12.6.6.1   -> 0.12.8.1.0  ) [CRAN]
#> coda         (0.19-4       -> 0.19-4.1    ) [CRAN]
#> stochvol     (3.2.3        -> 3.2.4       ) [CRAN]
#> bsvars       (23f015a04... -> d4dcc1e97...) [Git]
#> Installing 3 packages: RcppArmadillo, coda, stochvol
#> Installing packages into '/private/var/folders/_d/j7rqc71x2y54_69_6v6q33fr0000gn/T/RtmpuNc4T7/temp_libpathc879703baff7'
#> (as 'lib' is unspecified)
#> 
#> The downloaded binary packages are in
#>  /var/folders/_d/j7rqc71x2y54_69_6v6q33fr0000gn/T//RtmpKrQWw0/downloaded_packages
#> Downloading git repo https://github.com/bsvars/bsvars.git
#> '/usr/bin/git' clone --depth 1 --no-hardlinks https://github.com/bsvars/bsvars.git /var/folders/_d/j7rqc71x2y54_69_6v6q33fr0000gn/T//RtmpKrQWw0/filee4ff266f7188
#> RcppArmad... (0.12.8.0.0 -> 0.12.8.1.0) [CRAN]
#> stochvol     (3.2.3      -> 3.2.4     ) [CRAN]
#> Installing 2 packages: RcppArmadillo, stochvol
#> Installing packages into '/private/var/folders/_d/j7rqc71x2y54_69_6v6q33fr0000gn/T/RtmpuNc4T7/temp_libpathc879703baff7'
#> (as 'lib' is unspecified)
#> 
#> The downloaded binary packages are in
#>  /var/folders/_d/j7rqc71x2y54_69_6v6q33fr0000gn/T//RtmpKrQWw0/downloaded_packages
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/private/var/folders/_d/j7rqc71x2y54_69_6v6q33fr0000gn/T/RtmpKrQWw0/filee4ff266f7188/DESCRIPTION’ ... OK
#> * preparing ‘bsvars’:
#> * checking DESCRIPTION meta-information ... OK
#> * cleaning src
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * looking to see if a ‘data/datalist’ file should be added
#>   NB: this package now depends on R (>= 3.5.0)
#>   WARNING: Added dependency on R >= 3.5.0 because serialized objects in
#>   serialize/load version 3 cannot be read in older versions of R.
#>   File(s) containing such objects:
#>     ‘bsvars/data/us_fiscal_ex.rda’ ‘bsvars/data/us_fiscal_lsuw.rda’
#> * building ‘bsvars_2.1.0.9000.tar.gz’
#> Installing package into '/private/var/folders/_d/j7rqc71x2y54_69_6v6q33fr0000gn/T/RtmpuNc4T7/temp_libpathc879703baff7'
#> (as 'lib' is unspecified)
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/private/var/folders/_d/j7rqc71x2y54_69_6v6q33fr0000gn/T/RtmpKrQWw0/remotese4ff79fc617f/bsvars-bsvarSIGNs-a1cbe98/DESCRIPTION’ ... OK
#> * preparing ‘bsvarSIGNs’:
#> * checking DESCRIPTION meta-information ... OK
#> * cleaning src
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘bsvarSIGNs_0.1.9000.tar.gz’
#> Installing package into '/private/var/folders/_d/j7rqc71x2y54_69_6v6q33fr0000gn/T/RtmpuNc4T7/temp_libpathc879703baff7'
#> (as 'lib' is unspecified)
```

# Example

A repliaction of Antolín-Díaz & Rubio-Ramírez (2018).

    library(bsvarSIGNs)

    data(oil)

    # try ?specify_bsvarSIGN
    sign_narrative <- matrix(c(2, 0, 3, 2, 236, 0), ncol = 6)
    sign_irf       <- array(matrix(c(-1, -1, 1, 1, 1, 1, 1, -1, 1), nrow = 3), dim = c(3, 3, 1))

    specification  <- specify_bsvarSIGN$new(oil,
                                            p = 12,
                                            sign_irf = sign_irf,
                                            sign_narrative = sign_narrative)

    posterior      <- estimate(specification, S = 10000, thin = 10)
    irf            <- compute_irf(posterior, horizon = 24)

<!-- badges: start -->

[![R-CMD-check](https://github.com/bsvars/bsvarSIGNs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bsvars/bsvarSIGNs/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
