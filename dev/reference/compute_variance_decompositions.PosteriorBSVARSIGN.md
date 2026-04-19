# Computes posterior draws of the forecast error variance decomposition

Each of the draws from the posterior estimation of the model is
transformed into a draw from the posterior distribution of the forecast
error variance decomposition.

## Usage

``` r
# S3 method for class 'PosteriorBSVARSIGN'
compute_variance_decompositions(posterior, horizon)
```

## Arguments

- posterior:

  posterior estimation outcome - an object of class `PosteriorBSVARSIGN`
  obtained by running the `estimate` function.

- horizon:

  a positive integer number denoting the forecast horizon for the
  impulse responses computations.

## Value

An object of class `PosteriorFEVD`, that is, an `NxNx(horizon+1)xS`
array with attribute `PosteriorFEVD` containing `S` draws of the
forecast error variance decomposition.

## References

Kilian, L., & Lütkepohl, H. (2017). Structural VAR Tools, Chapter 4, In:
Structural vector autoregressive analysis. Cambridge University Press.

## See also

[`compute_impulse_responses.PosteriorBSVARSIGN`](https://bsvars.org/bsvarSIGNs/dev/reference/compute_impulse_responses.PosteriorBSVARSIGN.md),
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

# compute forecast error variance decomposition 2 years ahead
fevd           = compute_variance_decompositions(posterior, horizon = 8)

# workflow with the pipe |>
############################################################
set.seed(123)
optimism |>
  specify_bsvarSIGN$new(sign_irf = sign_irf) |> 
  estimate(S = 10) |> 
  compute_variance_decompositions(horizon = 8) -> fevd
#> **************************************************|
#>  bsvarSIGNs: Bayesian Structural VAR with sign,   |
#>              zero and narrative restrictions      |
#> **************************************************|
#>  Progress of simulation for 10 independent draws
#>  Press Esc to interrupt the computations
#> **************************************************|
  
```
