# R6 Class Representing Forecasts

R6 class representing draws from the predictive density of a Bayesian
Structural Vector Autoregression model.

## Value

An object of class `Forecasts`.

## Details

The class contains the following objects:

- `forecasts`:

  An `N x horizon x S` array containing draws from the predictive
  density.

- `forecast_mean`:

  An `N x horizon x S` array containing the conditional means of the
  predictive density.

- `forecast_covariance`:

  An `N x N x horizon x S` array containing the conditional covariance
  matrices of the predictive density.

- `Y`:

  An `N x T` matrix containing the data on the dependent variables used
  for estimation.

The method `as_list()` returns the contents of the `Forecasts` object as
a list.

## Public fields

- `forecasts`:

  An `N x horizon x S` numeric array containing draws from the
  predictive density.

- `forecast_mean`:

  An `N x horizon x S` numeric array containing the conditional means of
  the predictive density.

- `forecast_covariance`:

  An `N x N x horizon x S` numeric array containing the conditional
  covariance matrices of the predictive density.

- `Y`:

  An `N x T` numeric matrix containing the data on the dependent
  variables used for estimation.

## Methods

### Public methods

- [`Forecasts$new()`](#method-Forecasts-initialize)

- [`Forecasts$get_forecasts()`](#method-Forecasts-get_forecasts)

- [`Forecasts$clone()`](#method-Forecasts-clone)

------------------------------------------------------------------------

### `Forecasts$new()`

Creates a new `Forecasts` object from the output of the forecasting
procedure.

#### Usage

    Forecasts$new(output, Y)

#### Arguments

- `output`:

  A list containing the forecasting output, including `forecasts`,
  `forecast_mean`, and `forecast_cov`.

- `Y`:

  An `N x T` matrix containing the data on the dependent variables.

#### Returns

An object of class `Forecasts`.

------------------------------------------------------------------------

### `Forecasts$get_forecasts()`

Converts the `Forecasts` object to a list.

#### Usage

    Forecasts$get_forecasts()

#### Returns

A list containing `forecasts`, `forecast_mean`, `forecast_covariance`,
and `Y`.

------------------------------------------------------------------------

### `Forecasts$clone()`

The objects of this class are cloneable with this method.

#### Usage

    Forecasts$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
spec = specify_bsvarSIGN$new(optimism * 100)
post = estimate(spec, 5)
#> **************************************************|
#>  bsvarSIGNs: Bayesian Structural VAR with sign,   |
#>              zero and narrative restrictions      |
#> **************************************************|
#>  Progress of simulation for 5 independent draws
#>  Press Esc to interrupt the computations
#> **************************************************|
fore = forecast(post, 4)
apply(fore$forecasts, 1:2, mean) # compute mean forecasts 
#>               [,1]          [,2]          [,3]          [,4]
#> [1,]    84.5451248    85.1549035    84.8759131  8.445742e+01
#> [2,] -1070.4948614 -1069.6004138 -1070.8006405 -1.069942e+03
#> [3,]  -338.6833530  -338.5817412  -338.7693864 -3.387349e+02
#> [4,]     0.5422752     0.5167357     0.1908743  1.367001e-02
#> [5,]  -784.2297911  -784.9333610  -785.2907366 -7.851993e+02
```
