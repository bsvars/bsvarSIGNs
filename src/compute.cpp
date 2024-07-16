
#include <RcppArmadillo.h>
#include "bsvars.h"

using namespace Rcpp;
using namespace arma;

// cpp wrappers for bsvars functions

// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube structural_shocks (
    const arma::cube&     posterior_B,    // (N, N, S)
    const arma::cube&     posterior_A,    // (N, K, S)
    const arma::mat&      Y,              // NxT dependent variables
    const arma::mat&      X               // KxT dependent variables
) {
  cube structural_shocks = bsvars::bsvars_structural_shocks (posterior_B, posterior_A, Y, X);
  return structural_shocks;
} // END structural_shocks
