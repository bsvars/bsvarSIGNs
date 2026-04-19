# Forecasting using Structural Vector Autoregression

Samples from the joint predictive density of all of the dependent
variables for models from packages bsvars, bsvarSIGNs or bvarPANELs at
forecast horizons from 1 to `horizon` specified as an argument of the
function. Also facilitates forecasting using models with exogenous
variables and conditional forecasting given projected future
trajcetories of (some of the) variables.

## Usage

``` r
# S3 method for class 'PosteriorBSVARSIGN'
forecast(
  posterior,
  horizon = 1,
  exogenous_forecast = NULL,
  conditional_forecast = NULL
)
```

## Arguments

- posterior:

  posterior estimation outcome - an object of class `PosteriorBSVARSIGN`
  obtained by running the `estimate` function.

- horizon:

  a positive integer, specifying the forecasting horizon.

- exogenous_forecast:

  a matrix of dimension `horizon x d` containing forecasted values of
  the exogenous variables.

- conditional_forecast:

  a `horizon x N` matrix with forecasted values for selected variables.
  It should only contain `numeric` or `NA` values. The entries with `NA`
  values correspond to the values that are forecasted conditionally on
  the realisations provided as `numeric` values.

## Value

A list of class `Forecasts` containing the draws from the predictive
density and data. The output list includes element:

- forecasts:

  an `NxhorizonxS` array with the draws from predictive density

- Y:

  an \\NxT\\ matrix with the data on dependent variables

## See also

[`estimate.BSVARSIGN`](https://bsvars.org/bsvarSIGNs/dev/reference/estimate.BSVARSIGN.md),
[`summary`](https://rdrr.io/r/base/summary.html),
[`plot`](https://rdrr.io/r/graphics/plot.default.html)

## Author

Tomasz Woźniak <wozniak.tom@pm.me> and Xiaolei Wang
<adamwang15@gmail.com>

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

# sample from predictive density 1 year ahead
predictive     = forecast(posterior, 4)

# workflow with the pipe |>
############################################################
set.seed(123)
optimism |>
  specify_bsvarSIGN$new(sign_irf = sign_irf) |>
  estimate(S = 20) |> 
  forecast(horizon = 4) -> predictive
#> **************************************************|
#>  bsvarSIGNs: Bayesian Structural VAR with sign,   |
#>              zero and narrative restrictions      |
#> **************************************************|
#>  Progress of simulation for 20 independent draws
#>  Press Esc to interrupt the computations
#> **************************************************|

# conditional forecasting 2 quarters ahead conditioning on 
#  provided future values for the Gross Domestic Product 
############################################################
cf         = matrix(NA , 2, 5)
# # conditional forecasts equal to the last consumption observation
cf[,3]     = tail(optimism, 1)[3]
predictive = forecast(posterior, 2, conditional_forecast = cf)

# workflow with the pipe |>
############################################################
set.seed(123)
optimism |>
  specify_bsvarSIGN$new(sign_irf = sign_irf) |>
  estimate(S = 10) |> 
  forecast(horizon = 2, conditional_forecast = cf) -> predictive
#> **************************************************|
#>  bsvarSIGNs: Bayesian Structural VAR with sign,   |
#>              zero and narrative restrictions      |
#> **************************************************|
#>  Progress of simulation for 10 independent draws
#>  Press Esc to interrupt the computations
#> **************************************************|
```
