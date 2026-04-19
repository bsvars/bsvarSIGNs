# Computes posterior draws from data predictive density

Each of the draws from the posterior estimation of models from packages
bsvars or bsvarSIGNs is transformed into a draw from the data predictive
density.

## Usage

``` r
# S3 method for class 'PosteriorBSVARSIGN'
compute_fitted_values(posterior)
```

## Arguments

- posterior:

  posterior estimation outcome - an object of class `PosteriorBSVARSIGN`
  obtained by running the `estimate` function.

## Value

An object of class `PosteriorFitted`, that is, an `NxTxS` array with
attribute `PosteriorFitted` containing `S` draws from the data
predictive density.

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

# compute draws from in-sample predictive density
fitted         = compute_fitted_values(posterior)

# workflow with the pipe |>
############################################################
set.seed(123)
optimism |>
  specify_bsvarSIGN$new(sign_irf = sign_irf) |> 
  estimate(S = 20) |> 
  compute_fitted_values() -> fitted
#> **************************************************|
#>  bsvarSIGNs: Bayesian Structural VAR with sign,   |
#>              zero and narrative restrictions      |
#> **************************************************|
#>  Progress of simulation for 20 independent draws
#>  Press Esc to interrupt the computations
#> **************************************************|
```
