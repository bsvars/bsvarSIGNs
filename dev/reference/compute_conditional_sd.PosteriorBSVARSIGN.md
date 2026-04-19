# Computes posterior draws of structural shock conditional standard deviations

Each of the draws from the posterior estimation of models is transformed
into a draw from the posterior distribution of the structural shock
conditional standard deviations.

## Usage

``` r
# S3 method for class 'PosteriorBSVARSIGN'
compute_conditional_sd(posterior)
```

## Arguments

- posterior:

  posterior estimation outcome - an object of class `PosteriorBSVARSIGN`
  obtained by running the `estimate` function.

## Value

An object of class `PosteriorSigma`, that is, an `NxTxS` array with
attribute `PosteriorSigma` containing `S` draws of the structural shock
conditional standard deviations.

## See also

[`estimate.BSVARSIGN`](https://bsvars.org/bsvarSIGNs/dev/reference/estimate.BSVARSIGN.md)

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

# compute structural shocks' conditional standard deviations
sigma          = compute_conditional_sd(posterior)
#> The model is homoskedastic. Returning an NxTxS matrix of conditional sd all equal to 1.

# workflow with the pipe |>
############################################################
set.seed(123)
optimism |>
  specify_bsvarSIGN$new(sign_irf = sign_irf) |> 
  estimate(S = 10) |> 
  compute_conditional_sd() -> csd
#> **************************************************|
#>  bsvarSIGNs: Bayesian Structural VAR with sign,   |
#>              zero and narrative restrictions      |
#> **************************************************|
#>  Progress of simulation for 10 independent draws
#>  Press Esc to interrupt the computations
#> **************************************************|
#> The model is homoskedastic. Returning an NxTxS matrix of conditional sd all equal to 1.
```
