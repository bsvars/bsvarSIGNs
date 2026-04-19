# Package index

## bsvarSIGNs: Bayesian Estimation of Structural Vector Autoregressions Identified by Sign, Zero, and Narrative Restrictions

Browse package information

- [`bsvarSIGNs-package`](https://bsvars.org/bsvarSIGNs/dev/reference/bsvarSIGNs-package.md)
  [`bsvarSIGNs`](https://bsvars.org/bsvarSIGNs/dev/reference/bsvarSIGNs-package.md)
  : Bayesian Estimation of Structural Vector Autoregressions Identified
  by Sign, Zero, and Narrative Restrictions

## Data

Upload sample data set

- [`optimism`](https://bsvars.org/bsvarSIGNs/dev/reference/optimism.md)
  : A 5-variable US business cycle data, from 1955 Q1 to 2004 Q4
- [`monetary`](https://bsvars.org/bsvarSIGNs/dev/reference/monetary.md)
  : A 6-variable US monetary policy data, from 1965 Jan to 2007 Aug

## Model specification

Choose a model to work with

- [`specify_bsvarSIGN`](https://bsvars.org/bsvarSIGNs/dev/reference/specify_bsvarSIGN.md)
  : R6 Class representing the specification of the BSVARSIGN model
- [`specify_identification_bsvarSIGN`](https://bsvars.org/bsvarSIGNs/dev/reference/specify_identification_bsvarSIGN.md)
  : R6 Class Representing IdentificationBSVARSIGN
- [`specify_narrative()`](https://bsvars.org/bsvarSIGNs/dev/reference/specify_narrative.md)
  : vector specifying one narrative restriction
- [`specify_posterior_bsvarSIGN`](https://bsvars.org/bsvarSIGNs/dev/reference/specify_posterior_bsvarSIGN.md)
  : R6 Class Representing PosteriorBSVARSIGN
- [`specify_prior_bsvarSIGN`](https://bsvars.org/bsvarSIGNs/dev/reference/specify_prior_bsvarSIGN.md)
  : R6 Class Representing PriorBSVAR

## More detailed model specification

Adjust or inspect the specified model

- [`specify_identification_bsvarSIGN`](https://bsvars.org/bsvarSIGNs/dev/reference/specify_identification_bsvarSIGN.md)
  : R6 Class Representing IdentificationBSVARSIGN
- [`specify_prior_bsvarSIGN`](https://bsvars.org/bsvarSIGNs/dev/reference/specify_prior_bsvarSIGN.md)
  : R6 Class Representing PriorBSVAR

## Estimation

Run Bayesian estimation of your model and inspect the outputs

- [`estimate(`*`<BSVARSIGN>`*`)`](https://bsvars.org/bsvarSIGNs/dev/reference/estimate.BSVARSIGN.md)
  : Bayesian estimation of a Structural Vector Autoregression with
  traditional and narrative sign restrictions via Gibbs sampler
- [`specify_posterior_bsvarSIGN`](https://bsvars.org/bsvarSIGNs/dev/reference/specify_posterior_bsvarSIGN.md)
  : R6 Class Representing PosteriorBSVARSIGN

## Forecasting

Predict future values of your variables

- [`forecast(`*`<PosteriorBSVARSIGN>`*`)`](https://bsvars.org/bsvarSIGNs/dev/reference/forecast.PosteriorBSVARSIGN.md)
  : Forecasting using Structural Vector Autoregression

## Structural analyses

Compute interpretable outcomes

- [`compute_conditional_sd(`*`<PosteriorBSVARSIGN>`*`)`](https://bsvars.org/bsvarSIGNs/dev/reference/compute_conditional_sd.PosteriorBSVARSIGN.md)
  : Computes posterior draws of structural shock conditional standard
  deviations
- [`compute_fitted_values(`*`<PosteriorBSVARSIGN>`*`)`](https://bsvars.org/bsvarSIGNs/dev/reference/compute_fitted_values.PosteriorBSVARSIGN.md)
  : Computes posterior draws from data predictive density
- [`compute_historical_decompositions(`*`<PosteriorBSVARSIGN>`*`)`](https://bsvars.org/bsvarSIGNs/dev/reference/compute_historical_decompositions.PosteriorBSVARSIGN.md)
  : Computes posterior draws of historical decompositions
- [`compute_impulse_responses(`*`<PosteriorBSVARSIGN>`*`)`](https://bsvars.org/bsvarSIGNs/dev/reference/compute_impulse_responses.PosteriorBSVARSIGN.md)
  : Computes posterior draws of impulse responses
- [`compute_structural_shocks(`*`<PosteriorBSVARSIGN>`*`)`](https://bsvars.org/bsvarSIGNs/dev/reference/compute_structural_shocks.PosteriorBSVARSIGN.md)
  : Computes posterior draws of structural shocks
- [`compute_variance_decompositions(`*`<PosteriorBSVARSIGN>`*`)`](https://bsvars.org/bsvarSIGNs/dev/reference/compute_variance_decompositions.PosteriorBSVARSIGN.md)
  : Computes posterior draws of the forecast error variance
  decomposition

## Posterior summaries

Analyse the posterior summaries of the posterior estimation outcomes
using function [`summary()`](https://rdrr.io/r/base/summary.html)

## Plot your results

Prepare beautiful and informative plots for your analyses using function
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)
