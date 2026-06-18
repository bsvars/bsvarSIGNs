# R6 Class Representing PriorBSVAR

The class PriorBSVARSIGN presents a prior specification for the
homoskedastic bsvar model.

## Public fields

- `p`:

  a positive integer - the number of lags.

- `hyper`:

  a `(N+3)xS` matrix of hyper-parameters \\\mu, \delta, \lambda, \psi\\.

- `A`:

  a `NxK` normal prior mean matrix for the autoregressive parameters.

- `V`:

  a `KxK` matrix determining the normal prior column-specific covariance
  for the autoregressive parameters.

- `S`:

  an `NxN` matrix determining the inverted-Wishart prior scale of error
  terms covariance matrix.

- `nu`:

  a positive scalar greater than `N+1` - the shape of the
  inverted-Wishart prior for error terms covariance matrix.

- `data`:

  an `TxN` matrix of observations.

- `Y`:

  an `NxT` matrix of dependent variables.

- `X`:

  an `KxT` matrix of independent variables.

- `Ysoc`:

  an `NxN` matrix with the sum-of-coefficients dummy observations.

- `Xsoc`:

  an `KxN` matrix with the sum-of-coefficients dummy observations.

- `Ysur`:

  an `NxN` matrix with the single-unit-root dummy observations.

- `Xsur`:

  an `KxN` matrix with the single-unit-root dummy observations.

- `mu.scale`:

  a positive scalar - the shape of the gamma prior for \\\mu\\.

- `mu.shape`:

  a positive scalar - the shape of the gamma prior for \\\mu\\.

- `delta.scale`:

  a positive scalar - the shape of the gamma prior for \\\delta\\.

- `delta.shape`:

  a positive scalar - the shape of the gamma prior for \\\delta\\.

- `lambda.scale`:

  a positive scalar - the shape of the gamma prior for \\\lambda\\.

- `lambda.shape`:

  a positive scalar - the shape of the gamma prior for \\\lambda\\.

- `psi.scale`:

  a positive scalar - the shape of the inverted gamma prior for
  \\\psi\\.

- `psi.shape`:

  a positive scalar - the shape of the inverted gamma prior for
  \\\psi\\.

- `covid`:

  NULL or a positive integer indicating the start of the COVID-19
  pandemic.

## Methods

### Public methods

- [`specify_prior_bsvarSIGN$new()`](#method-PriorBSVARSIGN-new)

- [`specify_prior_bsvarSIGN$get_prior()`](#method-PriorBSVARSIGN-get_prior)

- [`specify_prior_bsvarSIGN$no_dummy()`](#method-PriorBSVARSIGN-no_dummy)

- [`specify_prior_bsvarSIGN$estimate_hyper()`](#method-PriorBSVARSIGN-estimate_hyper)

- [`specify_prior_bsvarSIGN$clone()`](#method-PriorBSVARSIGN-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new prior specification PriorBSVAR.

#### Usage

    specify_prior_bsvarSIGN$new(
      data,
      p,
      exogenous = NULL,
      stationary = rep(FALSE, ncol(data))
    )

#### Arguments

- `data`:

  the `TxN` data matrix of observations.

- `p`:

  a positive integer - the autoregressive lag order of the SVAR model.

- `exogenous`:

  a `Txd` matrix of exogenous variables.

- `stationary`:

  an `N` logical vector - its element set to `FALSE` sets the prior mean
  for the autoregressive parameters of the `N`th equation to the white
  noise process, otherwise to random walk.

#### Returns

A new prior specification PriorBSVARSIGN.

#### Examples

    # a prior for 5-variable example with one lag and stationary data
    data(optimism)
    prior = specify_prior_bsvarSIGN$new(optimism, p = 1)
    prior$B # show autoregressive prior mean

------------------------------------------------------------------------

### Method `get_prior()`

Returns the elements of the prior specification PriorBSVAR as a `list`.

#### Usage

    specify_prior_bsvarSIGN$get_prior()

#### Examples

    # a prior for 5-variable example with four lags
    prior = specify_prior_bsvar$new(N = 5, p = 4)
    prior$get_prior() # show the prior as list

------------------------------------------------------------------------

### Method `no_dummy()`

Sets the sum-of-coefficients and single-unit-root dummy observations to
zero (removes the dummy observation prior).

#### Usage

    specify_prior_bsvarSIGN$no_dummy()

#### Examples

    # a prior for 5-variable example with four lags
    data(optimism)
    prior = specify_prior_bsvarSIGN$new(optimism, p = 4)
    prior$no_dummy() # remove dummy observations

------------------------------------------------------------------------

### Method `estimate_hyper()`

Estimates hyper-parameters with adaptive Metropolis algorithm.

#### Usage

    specify_prior_bsvarSIGN$estimate_hyper(
      S = 10000,
      burn_in = S/2,
      mu = TRUE,
      delta = TRUE,
      lambda = TRUE,
      psi = TRUE,
      covid = NULL
    )

#### Arguments

- `S`:

  number of MCMC draws.

- `burn_in`:

  number of burn-in draws.

- `mu`:

  whether to estimate the hyper-parameter in the sum-of-coefficients
  dummy prior.

- `delta`:

  whether to estimate the hyper-parameter in the single-unit-root dummy
  prior.

- `lambda`:

  whether to estimate the hyper-parameter of the shrinkage in the
  Minnesota prior.

- `psi`:

  whether to estimate the hyper-parameter of the variances in the
  Minnesota prior.

- `covid`:

  NULL or positive integer indicating the start of the COVID-19 pandemic

#### Examples

    # specify the model and set seed
    set.seed(123)
    data(optimism)
    prior = specify_prior_bsvarSIGN$new(optimism, p = 4)

    # estimate hyper parameters with adaptive Metropolis algorithm
    prior$estimate_hyper(S = 10, psi = TRUE)

    # trace plot
    hyper = t(prior$hyper)[, 4:8]
    colnames(hyper) = paste("psi", 1:5, sep = "")
    plot.ts(hyper)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    specify_prior_bsvarSIGN$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# a prior for 5-variable example with one lag
data(optimism)
prior = specify_prior_bsvarSIGN$new(optimism, p = 1)
prior$A  # show autoregressive prior mean
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    0    0    0    0    0
#> [2,]    0    1    0    0    0    0
#> [3,]    0    0    1    0    0    0
#> [4,]    0    0    0    1    0    0
#> [5,]    0    0    0    0    1    0


## ------------------------------------------------
## Method `specify_prior_bsvarSIGN$new`
## ------------------------------------------------

# a prior for 5-variable example with one lag and stationary data
data(optimism)
prior = specify_prior_bsvarSIGN$new(optimism, p = 1)
prior$B # show autoregressive prior mean
#> NULL


## ------------------------------------------------
## Method `specify_prior_bsvarSIGN$get_prior`
## ------------------------------------------------

# a prior for 5-variable example with four lags
prior = specify_prior_bsvar$new(N = 5, p = 4)
prior$get_prior() # show the prior as list
#> $A
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
#> [1,]    1    0    0    0    0    0    0    0    0     0     0     0     0     0
#> [2,]    0    1    0    0    0    0    0    0    0     0     0     0     0     0
#> [3,]    0    0    1    0    0    0    0    0    0     0     0     0     0     0
#> [4,]    0    0    0    1    0    0    0    0    0     0     0     0     0     0
#> [5,]    0    0    0    0    1    0    0    0    0     0     0     0     0     0
#>      [,15] [,16] [,17] [,18] [,19] [,20] [,21]
#> [1,]     0     0     0     0     0     0     0
#> [2,]     0     0     0     0     0     0     0
#> [3,]     0     0     0     0     0     0     0
#> [4,]     0     0     0     0     0     0     0
#> [5,]     0     0     0     0     0     0     0
#> 
#> $A_V_inv
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#>  [1,]    1    0    0    0    0    0    0    0    0     0     0     0     0
#>  [2,]    0    1    0    0    0    0    0    0    0     0     0     0     0
#>  [3,]    0    0    1    0    0    0    0    0    0     0     0     0     0
#>  [4,]    0    0    0    1    0    0    0    0    0     0     0     0     0
#>  [5,]    0    0    0    0    1    0    0    0    0     0     0     0     0
#>  [6,]    0    0    0    0    0    4    0    0    0     0     0     0     0
#>  [7,]    0    0    0    0    0    0    4    0    0     0     0     0     0
#>  [8,]    0    0    0    0    0    0    0    4    0     0     0     0     0
#>  [9,]    0    0    0    0    0    0    0    0    4     0     0     0     0
#> [10,]    0    0    0    0    0    0    0    0    0     4     0     0     0
#> [11,]    0    0    0    0    0    0    0    0    0     0     9     0     0
#> [12,]    0    0    0    0    0    0    0    0    0     0     0     9     0
#> [13,]    0    0    0    0    0    0    0    0    0     0     0     0     9
#> [14,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [15,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [16,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [17,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [18,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [19,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [20,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#> [21,]    0    0    0    0    0    0    0    0    0     0     0     0     0
#>       [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21]
#>  [1,]     0     0     0     0     0     0     0     0
#>  [2,]     0     0     0     0     0     0     0     0
#>  [3,]     0     0     0     0     0     0     0     0
#>  [4,]     0     0     0     0     0     0     0     0
#>  [5,]     0     0     0     0     0     0     0     0
#>  [6,]     0     0     0     0     0     0     0     0
#>  [7,]     0     0     0     0     0     0     0     0
#>  [8,]     0     0     0     0     0     0     0     0
#>  [9,]     0     0     0     0     0     0     0     0
#> [10,]     0     0     0     0     0     0     0     0
#> [11,]     0     0     0     0     0     0     0     0
#> [12,]     0     0     0     0     0     0     0     0
#> [13,]     0     0     0     0     0     0     0     0
#> [14,]     9     0     0     0     0     0     0     0
#> [15,]     0     9     0     0     0     0     0     0
#> [16,]     0     0    16     0     0     0     0     0
#> [17,]     0     0     0    16     0     0     0     0
#> [18,]     0     0     0     0    16     0     0     0
#> [19,]     0     0     0     0     0    16     0     0
#> [20,]     0     0     0     0     0     0    16     0
#> [21,]     0     0     0     0     0     0     0     1
#> 
#> $B_V_inv
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    0    0    0    0
#> [2,]    0    1    0    0    0
#> [3,]    0    0    1    0    0
#> [4,]    0    0    0    1    0
#> [5,]    0    0    0    0    1
#> 
#> $B_nu
#> [1] 5
#> 
#> $hyper_nu_B
#> [1] 10
#> 
#> $hyper_a_B
#> [1] 10
#> 
#> $hyper_s_BB
#> [1] 100
#> 
#> $hyper_nu_BB
#> [1] 1
#> 
#> $hyper_nu_A
#> [1] 10
#> 
#> $hyper_a_A
#> [1] 10
#> 
#> $hyper_s_AA
#> [1] 10
#> 
#> $hyper_nu_AA
#> [1] 10
#> 


## ------------------------------------------------
## Method `specify_prior_bsvarSIGN$no_dummy`
## ------------------------------------------------

# a prior for 5-variable example with four lags
data(optimism)
prior = specify_prior_bsvarSIGN$new(optimism, p = 4)
prior$no_dummy() # remove dummy observations


## ------------------------------------------------
## Method `specify_prior_bsvarSIGN$estimate_hyper`
## ------------------------------------------------

# specify the model and set seed
set.seed(123)
data(optimism)
prior = specify_prior_bsvarSIGN$new(optimism, p = 4)

# estimate hyper parameters with adaptive Metropolis algorithm
prior$estimate_hyper(S = 10, psi = TRUE)
#> **************************************************|
#>  Adaptive Metropolis MCMC: hyper parameters       |
#> **************************************************|

# trace plot
hyper = t(prior$hyper)[, 4:8]
colnames(hyper) = paste("psi", 1:5, sep = "")
plot.ts(hyper)

```
