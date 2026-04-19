# Computes posterior draws of historical decompositions

Each of the draws from the posterior estimation of models from packages
bsvars or bsvarSIGNs is transformed into a draw from the posterior
distribution of the historical decompositions. IMPORTANT! The historical
decompositions are interpreted correctly for covariance stationary data.
Application to unit-root non-stationary data might result in
non-interpretable outcomes.

## Usage

``` r
# S3 method for class 'PosteriorBSVARSIGN'
compute_historical_decompositions(posterior, show_progress = TRUE)
```

## Arguments

- posterior:

  posterior estimation outcome - an object of class `PosteriorBSVARSIGN`
  obtained by running the `estimate` function.

- show_progress:

  a logical value, if `TRUE` the estimation progress bar is visible

## Value

An object of class `PosteriorHD`, that is, an `NxNxTxS` array with
attribute `PosteriorHD` containing `S` draws of the historical
decompositions.

## References

Kilian, L., & Lütkepohl, H. (2017). Structural VAR Tools, Chapter 4, In:
Structural vector autoregressive analysis. Cambridge University Press.

## See also

[`estimate.BSVARSIGN`](https://bsvars.org/bsvarSIGNs/dev/reference/estimate.BSVARSIGN.md),
[`summary`](https://rdrr.io/r/base/summary.html),
[`plot`](https://rdrr.io/r/graphics/plot.default.html)

## Author

Xiaolei Wang <adamwang15@gmail.com> and Tomasz Woźniak
<wozniak.tom@pm.me>

## Examples

``` r
# upload data
data(optimism)

# specify the model and set seed
set.seed(123)

# + no effect on productivity (zero restriction)
# + positive effect on stock prices (positive sign restriction) 
sign_irf       = matrix(c(0, 1, rep(NA, 23)), 5, 5)
specification  = specify_bsvarSIGN$new(optimism, sign_irf = sign_irf)

# estimate the model
posterior      = estimate(specification, 10)
#> **************************************************|
#>  bsvarSIGNs: Bayesian Structural VAR with sign,   |
#>              zero and narrative restrictions      |
#> **************************************************|
#>  Progress of simulation for 10 independent draws
#>  Press Esc to interrupt the computations
#> **************************************************|

# compute historical decompositions
hd            = compute_historical_decompositions(posterior)
#> **************************************************|
#> bsvars: Bayesian Structural Vector Autoregressions|
#> **************************************************|
#>  Computing historical decomposition               |
#> **************************************************|
#>  This might take a little while :)                
#> **************************************************|

# workflow with the pipe |>
############################################################
set.seed(123)
optimism |>
  specify_bsvarSIGN$new(sign_irf = sign_irf) |> 
  estimate(S = 10) |> 
  compute_historical_decompositions() -> hd
#> **************************************************|
#>  bsvarSIGNs: Bayesian Structural VAR with sign,   |
#>              zero and narrative restrictions      |
#> **************************************************|
#>  Progress of simulation for 10 independent draws
#>  Press Esc to interrupt the computations
#> **************************************************|
#> **************************************************|
#> bsvars: Bayesian Structural Vector Autoregressions|
#> **************************************************|
#>  Computing historical decomposition               |
#> **************************************************|
#>  This might take a little while :)                
#> **************************************************|
  
```
