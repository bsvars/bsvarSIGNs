# Computes posterior draws of impulse responses

Each of the draws from the posterior estimation of models from packages
bsvars or bsvarSIGNs is transformed into a draw from the posterior
distribution of the impulse responses.

## Usage

``` r
# S3 method for class 'PosteriorBSVARSIGN'
compute_impulse_responses(posterior, horizon, standardise = FALSE)
```

## Arguments

- posterior:

  posterior estimation outcome - an object of class `PosteriorBSVARSIGN`
  obtained by running the `estimate` function.

- horizon:

  a positive integer number denoting the forecast horizon for the
  impulse responses computations.

- standardise:

  a logical value. If `TRUE`, the impulse responses are standardised so
  that the variables' own shocks at horizon 0 are equal to 1. Otherwise,
  the parameter estimates determine this magnitude.

## Value

An object of class PosteriorIR, that is, an `NxNx(horizon+1)xS` array
with attribute PosteriorIR containing `S` draws of the impulse
responses.

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

# compute impulse responses 2 years ahead
irf           = compute_impulse_responses(posterior, horizon = 8)
#> Argument standardise is forcibly set to FALSE due to zero restrictions imposed on the diagonal element(s) of the on-impact impulse response matrix.

# workflow with the pipe |>
############################################################
set.seed(123)
optimism |>
  specify_bsvarSIGN$new(sign_irf = sign_irf) |> 
  estimate(S = 10) |> 
  compute_impulse_responses(horizon = 8) -> ir
#> **************************************************|
#>  bsvarSIGNs: Bayesian Structural VAR with sign,   |
#>              zero and narrative restrictions      |
#> **************************************************|
#>  Progress of simulation for 10 independent draws
#>  Press Esc to interrupt the computations
#> **************************************************|
#> Argument standardise is forcibly set to FALSE due to zero restrictions imposed on the diagonal element(s) of the on-impact impulse response matrix.

```
