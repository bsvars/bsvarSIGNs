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
#>              [,1]          [,2]         [,3]         [,4]
#> [1,]    84.482668  8.462837e+01    84.725136    84.313684
#> [2,] -1068.752173 -1.074758e+03 -1074.418499 -1079.844384
#> [3,]  -338.421656 -3.386165e+02  -338.746547  -338.466008
#> [4,]    -1.057242  3.738026e-02     1.401288     1.249104
#> [5,]  -784.616510 -7.848578e+02  -786.006054  -786.205879
```
