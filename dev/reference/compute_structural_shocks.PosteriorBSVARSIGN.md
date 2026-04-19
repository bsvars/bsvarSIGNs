# Computes posterior draws of structural shocks

Each of the draws from the posterior estimation of models from packages
bsvars or bsvarSIGNs is transformed into a draw from the posterior
distribution of the structural shocks.

## Usage

``` r
# S3 method for class 'PosteriorBSVARSIGN'
compute_structural_shocks(posterior)
```

## Arguments

- posterior:

  posterior estimation outcome - an object of class `PosteriorBSVARSIGN`
  obtained by running the `estimate` function.

## Value

An object of class `PosteriorShocks`, that is, an `NxTxS` array with
attribute `PosteriorShocks` containing `S` draws of the structural
shocks.

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

# compute structural shocks
shocks         = compute_structural_shocks(posterior)

# workflow with the pipe |>
############################################################
set.seed(123)
optimism |>
  specify_bsvarSIGN$new(sign_irf = sign_irf) |> 
  estimate(S = 20) |> 
  compute_structural_shocks() -> ss
#> **************************************************|
#>  bsvarSIGNs: Bayesian Structural VAR with sign,   |
#>              zero and narrative restrictions      |
#> **************************************************|
#>  Progress of simulation for 20 independent draws
#>  Press Esc to interrupt the computations
#> **************************************************|
```
