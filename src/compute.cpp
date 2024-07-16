
#include <RcppArmadillo.h>
#include "bsvars.h"

using namespace Rcpp;
using namespace arma;

// cpp wrappers for bsvars functions

// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube bsvarSIGNs_structural_shocks (
    const arma::cube&     posterior_B,    // (N, N, S)
    const arma::cube&     posterior_A,    // (N, K, S)
    const arma::mat&      Y,              // NxT dependent variables
    const arma::mat&      X               // KxT dependent variables
) {
  cube structural_shocks = bsvars::bsvars_structural_shocks (posterior_B, posterior_A, Y, X);
  return structural_shocks;
} // END bsvarSIGNs_structural_shocks



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube bsvarSIGNs_fitted_values (
    arma::cube&     posterior_A,        // NxKxS
    arma::cube&     posterior_B,        // NxNxS
    arma::cube&     posterior_sigma,    // NxTxS
    arma::mat&      X                   // KxT
) {
  cube fitted_values = bsvars::bsvars_fitted_values (posterior_A, posterior_B, posterior_sigma, X);
  return fitted_values;
} // END bsvarSIGNs_fitted_values