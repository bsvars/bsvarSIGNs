# Bayesian estimation of a Structural Vector Autoregression with traditional and narrative sign restrictions via Gibbs sampler

Estimates Bayesian Structural Vector Autoregression model using the
Gibbs sampler proposed by Waggoner & Zha (2003) with traditional sign
restrictions following Rubio-Ramírez, Waggoner & Zha (2010) and
narrative sign restrictions following Antolín-Díaz & Rubio-Ramírez
(2018). Additionally, the parameter matrices \\A\\ and \\B\\ follow a
Minnesota prior and generalised-normal prior distributions respectively
with the matrix-specific overall shrinkage parameters estimated using a
hierarchical prior distribution.

Given sign restrictions, in each Gibbs sampler iteration, the sampler
draws rotation matrix \\Q\\ uniformly from the space of `NxN` orthogonal
matrices and checks if the sign restrictions are satisfied. If a valid
\\Q\\ is found within `max_tries` (defined in `specify_bsvarSIGN`), the
sampler saves the current \\A\\ and \\B\\ draw and proceeds to the next
iteration. Otherwise, the sampler then proceeds to next iteration
without saving the current \\A\\ and \\B\\ draw. If a narrative sign
restriction is given, the posterior draws are resampled with
`algorithm 1` in Antolín-Díaz & Rubio-Ramírez (2018).

See section **Details** for the model equations.

## Usage

``` r
# S3 method for class 'BSVARSIGN'
estimate(specification, S, thin = 1, show_progress = TRUE)
```

## Arguments

- specification:

  an object of class BSVARSIGN generated using the
  `specify_bsvarSIGN$new()` function.

- S:

  a positive integer, the number of posterior draws to be generated

- thin:

  a positive integer, specifying the frequency of MCMC output thinning

- show_progress:

  a logical value, if `TRUE` the estimation progress bar is visible

## Value

An object of class `PosteriorBSVARSIGN` containing the Bayesian
estimation output and containing two elements:

`posterior` a list with a collection of `S` draws from the posterior
distribution generated via Gibbs sampler containing:

- A:

  an `NxKxS` array with the posterior draws for matrix \\A\\

- B:

  an `NxNxS` array with the posterior draws for matrix \\B\\

- hyper:

  a `5xS` matrix with the posterior draws for the hyper-parameters of
  the hierarchical prior distribution

- skipped:

  an integer of the total skipped iterations, the Gibbs sampler performs
  a total of S+skipped iterations, when the sampler does not find a
  valid rotation matrix `Q` within `max_tries`, the current iteration is
  skipped (i.e. the current draw of `A,B` is not saved). A message is
  shown when skipped/(skipped+S/thin) \> 0.05, where S/thin is the total
  number of draws returned.

`last_draw` an object of class BSVARSIGN with the last draw of the
current MCMC run as the starting value to be passed to the continuation
of the MCMC estimation using `estimate()`.

## Details

The Structural VAR model is given by the reduced form equation: \$\$Y =
AX + E\$\$ where \\Y\\ is an `NxT` matrix of dependent variables, \\X\\
is a `KxT` matrix of explanatory variables, \\E\\ is an `NxT` matrix of
reduced form error terms, and \\A\\ is an `NxK` matrix of autoregressive
slope coefficients and parameters on deterministic terms in \\X\\.

The structural equation is given by \$\$BE = U\$\$ where \\U\\ is an
`NxT` matrix of structural form error terms, and \\B\\ is an `NxN`
matrix of contemporaneous relationships. More specifically, \$\$B =
Q'P\$\$ where \\Q\\ is an `NxN` rotation matrix and \\P\\ is an `NxN`
lower triangular matrix.

Finally, the structural shocks, `U`, are temporally and
contemporaneously independent and jointly normally distributed with zero
mean and unit variances.

## References

Antolín-Díaz & Rubio-Ramírez (2018) Narrative Sign Restrictions for
SVARs, American Economic Review, 108(10), 2802-29,
\<doi:10.1257/aer.20161852\>.

Arias, Rubio-Ramírez, & Waggoner (2018), Inference Based on Structural
Vector Autoregressions Identified With Sign and Zero Restrictions:
Theory and Applications, Econometrica, 86(2), 685-720,
\<doi:10.3982/ECTA14468\>.

Giannone, Lenza, Primiceri (2015) Prior Selection for Vector
Autoregressions, Review of Economics and Statistics, 97(2), 436-451
\<doi:10.1162/REST_a_00483\>.

Rubio-Ramírez, Waggoner & Zha (2010) Structural Vector Autoregressions:
Theory of Identification and Algorithms for Inference, The Review of
Economic Studies, 77(2), 665-696,
\<doi:10.1111/j.1467-937X.2009.00578.x\>.

## Author

Tomasz Woźniak <wozniak.tom@pm.me>, Xiaolei Wang <adamwang15@gmail.com>

## Examples

``` r
# investigate the effects of the optimism shock
data(optimism)

# specify identifying restrictions:
# + no effect on productivity (zero restriction)
# + positive effect on stock prices (positive sign restriction) 
sign_irf       = matrix(c(0, 1, rep(NA, 23)), 5, 5)

# specify the model and set seed
set.seed(123)
specification  = specify_bsvarSIGN$new(optimism * 100,
                                       p        = 12,
                                       sign_irf = sign_irf)
                                       
# estimate the model
posterior      = estimate(specification, S = 10)
#> **************************************************|
#>  bsvarSIGNs: Bayesian Structural VAR with sign,   |
#>              zero and narrative restrictions      |
#> **************************************************|
#>  Progress of simulation for 10 independent draws
#>  Press Esc to interrupt the computations
#> **************************************************|
```
