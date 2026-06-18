# R6 Class representing the specification of the BSVARSIGN model

The class BSVARSIGN presents complete specification for the Bayesian
Structural VAR model with sign and narrative restrictions.

## See also

[`estimate.BSVARSIGN`](https://bsvars.org/bsvarSIGNs/dev/reference/estimate.BSVARSIGN.md),
[`specify_posterior_bsvarSIGN`](https://bsvars.org/bsvarSIGNs/dev/reference/specify_posterior_bsvarSIGN.md)

## Public fields

- `p`:

  a non-negative integer specifying the autoregressive lag order of the
  model.

- `identification`:

  an object IdentificationBSVARSIGN with the identifying restrictions.

- `prior`:

  an object PriorBSVARSIGN with the prior specification.

- `data_matrices`:

  an object DataMatricesBSVARSIGN with the data matrices.

- `starting_values`:

  an object StartingValuesBSVARSIGN with the starting values.

## Methods

### Public methods

- [`specify_bsvarSIGN$new()`](#method-BSVARSIGN-new)

- [`specify_bsvarSIGN$get_data_matrices()`](#method-BSVARSIGN-get_data_matrices)

- [`specify_bsvarSIGN$no_dummy_observations()`](#method-BSVARSIGN-no_dummy_observations)

- [`specify_bsvarSIGN$estimate_hyper()`](#method-BSVARSIGN-estimate_hyper)

- [`specify_bsvarSIGN$get_identification()`](#method-BSVARSIGN-get_identification)

- [`specify_bsvarSIGN$get_prior()`](#method-BSVARSIGN-get_prior)

- [`specify_bsvarSIGN$get_starting_values()`](#method-BSVARSIGN-get_starting_values)

- [`specify_bsvarSIGN$clone()`](#method-BSVARSIGN-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new specification of the Bayesian Structural VAR model with
sign and narrative restrictions BSVARSIGN.

#### Usage

    specify_bsvarSIGN$new(
      data,
      p = 1L,
      sign_irf,
      sign_narrative,
      sign_structural,
      max_tries = Inf,
      exogenous = NULL,
      stationary = rep(FALSE, ncol(data)),
      hyper_mu = TRUE,
      hyper_delta = TRUE,
      hyper_lambda = TRUE,
      hyper_psi = TRUE,
      hyper_covid = NULL
    )

#### Arguments

- `data`:

  a `(T+p)xN` matrix with time series data.

- `p`:

  a positive integer providing model's autoregressive lag order.

- `sign_irf`:

  a `NxNxH` array - sign and zero restrictions on the impulse response
  functions, ±1 for positive/negative sign restriction 0 for zero
  restrictions and NA for no restrictions, the `h`-th slice `NxN` matrix
  contains the restrictions on the `h-1` horizon.

- `sign_narrative`:

  a list of objects of class "narrative" - narrative sign restrictions.

- `sign_structural`:

  a `NxN` matrix with entries ±1 or NA - sign restrictions on the
  contemporaneous relations `B` between reduced-form errors `E` and
  structural shocks `U` where `BE=U`.

- `max_tries`:

  a positive integer with the maximum number of iterations for finding a
  rotation matrix \\Q\\ that would satisfy sign restrictions

- `exogenous`:

  a `(T+p)xd` matrix of exogenous variables.

- `stationary`:

  an `N` logical vector - its element set to `FALSE` sets the prior mean
  for the autoregressive parameters of the `N`th equation to the white
  noise process, otherwise to random walk.

- `hyper_mu`:

  whether to estimate the hyper-parameter in the sum-of-coefficients
  dummy prior.

- `hyper_delta`:

  whether to estimate the hyper-parameter in the single-unit-root dummy
  prior.

- `hyper_lambda`:

  whether to estimate the hyper-parameter of the shrinkage in the
  Minnesota prior.

- `hyper_psi`:

  whether to estimate the hyper-parameter of the variances in the
  Minnesota prior.

- `hyper_covid`:

  NULL or positive integer indicating the start of the COVID-19
  pandemic.

#### Returns

A new complete specification for the Bayesian Structural VAR model
BSVARSIGN.

------------------------------------------------------------------------

### Method `get_data_matrices()`

Returns the data matrices as the DataMatricesBSVAR object.

#### Usage

    specify_bsvarSIGN$get_data_matrices()

#### Examples

    # specify a model with the optimism data and 4 lags

    data(optimism)
    spec = specify_bsvarSIGN$new(
       data = optimism,
       p = 4
    )

    # get the data matrices
    spec$get_data_matrices()

------------------------------------------------------------------------

### Method `no_dummy_observations()`

Sets the sum-of-coefficients and single-unit-root dummy observations to
zero (removes the dummy observation prior).

#### Usage

    specify_bsvarSIGN$no_dummy_observations()

#### Examples

    # specify the model
    data(optimism)
    spec = specify_bsvarSIGN$new(optimism, p = 4)
    spec$no_dummy_observations() # remove dummy observations

------------------------------------------------------------------------

### Method `estimate_hyper()`

Estimates hyper-parameters with adaptive Metropolis algorithm.

#### Usage

    specify_bsvarSIGN$estimate_hyper(S = 10000, burn_in = S/2)

#### Arguments

- `S`:

  number of MCMC draws.

- `burn_in`:

  number of burn-in draws.

#### Examples

    # specify the model and set seed
    set.seed(123)
    data(optimism)
    spec = specify_bsvarSIGN$new(optimism, p = 4)

    # estimate hyper parameters with adaptive Metropolis algorithm
    spec$estimate_hyper(S = 10)

    # trace plot
    hyper = t(spec$prior$hyper)[, 4:8]
    colnames(hyper) = paste("psi", 1:5, sep = "")
    plot.ts(hyper)

------------------------------------------------------------------------

### Method `get_identification()`

Returns the identifying restrictions as the IdentificationBSVARSIGN
object.

#### Usage

    specify_bsvarSIGN$get_identification()

#### Examples

    # specify a model with the optimism data and 4 lags
    data(optimism)
    spec = specify_bsvarSIGN$new(
       data = optimism,
       p = 4
    )

    # get the identifying restrictions
    spec$get_identification()

------------------------------------------------------------------------

### Method `get_prior()`

Returns the prior specification as the PriorBSVAR object.

#### Usage

    specify_bsvarSIGN$get_prior()

#### Examples

    # specify a model with the optimism data and 4 lags

    data(optimism)
    spec = specify_bsvarSIGN$new(
       data = optimism,
       p = 4
    )

    # get the prior specification
    spec$get_prior()

------------------------------------------------------------------------

### Method `get_starting_values()`

Returns the starting values as the StartingValuesBSVAR object.

#### Usage

    specify_bsvarSIGN$get_starting_values()

#### Examples

    # specify a model with the optimism data and 4 lags

    data(optimism)
    spec = specify_bsvarSIGN$new(
       data = optimism,
       p = 4
    )

    # get the starting values
    spec$get_starting_values()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    specify_bsvarSIGN$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# specify a model with the optimism data and 4 lags

data(optimism)
specification = specify_bsvarSIGN$new(
   data = optimism,
   p = 4
)


## ------------------------------------------------
## Method `specify_bsvarSIGN$get_data_matrices`
## ------------------------------------------------

# specify a model with the optimism data and 4 lags

data(optimism)
spec = specify_bsvarSIGN$new(
   data = optimism,
   p = 4
)

# get the data matrices
spec$get_data_matrices()
#> <DataMatricesBSVAR>
#>   Public:
#>     X: 0.200390767 -11.08316654 -4.301382328 0.012867114 -7.572 ...
#>     Y: 0.194524324 -11.02210681 -4.293947793 0.024503568 -7.577 ...
#>     clone: function (deep = FALSE) 
#>     get_data_matrices: function () 
#>     initialize: function (data, p = 1L, exogenous = NULL) 


## ------------------------------------------------
## Method `specify_bsvarSIGN$no_dummy_observations`
## ------------------------------------------------

# specify the model
data(optimism)
spec = specify_bsvarSIGN$new(optimism, p = 4)
spec$no_dummy_observations() # remove dummy observations


## ------------------------------------------------
## Method `specify_bsvarSIGN$estimate_hyper`
## ------------------------------------------------

# specify the model and set seed
set.seed(123)
data(optimism)
spec = specify_bsvarSIGN$new(optimism, p = 4)

# estimate hyper parameters with adaptive Metropolis algorithm
spec$estimate_hyper(S = 10)
#> **************************************************|
#>  Adaptive Metropolis MCMC: hyper parameters       |
#> **************************************************|

# trace plot
hyper = t(spec$prior$hyper)[, 4:8]
colnames(hyper) = paste("psi", 1:5, sep = "")
plot.ts(hyper)



## ------------------------------------------------
## Method `specify_bsvarSIGN$get_identification`
## ------------------------------------------------

# specify a model with the optimism data and 4 lags
data(optimism)
spec = specify_bsvarSIGN$new(
   data = optimism,
   p = 4
)

# get the identifying restrictions
spec$get_identification()
#> <IdentificationBSVARSIGN>
#>   Public:
#>     VB: list
#>     clone: function (deep = FALSE) 
#>     get_identification: function () 
#>     initialize: function (N, sign_irf, sign_narrative, sign_structural, max_tries = Inf) 
#>     max_tries: Inf
#>     set_identification: function (N, sign_irf, sign_narrative, sign_structural) 
#>     sign_irf: NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA ...
#>     sign_narrative: list
#>     sign_structural: 1 NA NA NA NA NA 1 NA NA NA NA NA 1 NA NA NA NA NA 1 NA  ...


## ------------------------------------------------
## Method `specify_bsvarSIGN$get_prior`
## ------------------------------------------------

# specify a model with the optimism data and 4 lags

data(optimism)
spec = specify_bsvarSIGN$new(
   data = optimism,
   p = 4
)

# get the prior specification
spec$get_prior()
#> <PriorBSVARSIGN>
#>   Public:
#>     A: 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0  ...
#>     S: 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1
#>     V: 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0  ...
#>     X: 0.200390767 -11.08316654 -4.301382328 0.012867114 -7.572 ...
#>     Xsoc: 0.201757572 0 0 0 0 0.201757572 0 0 0 0 0.201757572 0 0  ...
#>     Xsur: 0.201757572 -11.0726542625 -4.2879895105 0.00226154775 - ...
#>     Y: 0.194524324 -11.02210681 -4.293947793 0.024503568 -7.577 ...
#>     Ysoc: 0.201757572 0 0 0 0 0 -11.0726542625 0 0 0 0 0 -4.287989 ...
#>     Ysur: 0.201757572 -11.0726542625 -4.2879895105 0.00226154775 - ...
#>     clone: function (deep = FALSE) 
#>     covid: NULL
#>     data: NA
#>     delta.scale: 0.618033988749895
#>     delta.shape: 2.61803398874989
#>     get_prior: function () 
#>     hyper: 1 1 0.2 6.63317874235789e-05 0.00655120654453227 1.68010 ...
#>     initialize: function (data, p, exogenous = NULL, stationary = rep(FALSE, 
#>     lambda.scale: 0.540312423743285
#>     lambda.shape: 1.37015621187164
#>     mu.scale: 0.618033988749895
#>     mu.shape: 2.61803398874989
#>     nu: 7
#>     p: 4
#>     psi.scale: 0.000799362037821841
#>     psi.shape: 0.998405094554603


## ------------------------------------------------
## Method `specify_bsvarSIGN$get_starting_values`
## ------------------------------------------------

# specify a model with the optimism data and 4 lags

data(optimism)
spec = specify_bsvarSIGN$new(
   data = optimism,
   p = 4
)

# get the starting values
spec$get_starting_values()
#> NULL
```
