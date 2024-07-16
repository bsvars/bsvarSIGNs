
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




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube ir1_cpp (
    const arma::mat& B,           // KxN
    const arma::mat& Theta0,      // NxN
    int              horizon,
    const int&       p
) {
  
  horizon++;
  int N         = Theta0.n_cols;
  
  cube irf      = zeros(N, N, horizon);
  irf.slice(0)  = Theta0;
  
  for(int t = 2; t <= horizon; t++) {
    for(int j = 1; j <= std::min(p, t - 1); j++) {
      mat B_j           = B.rows((j - 1) * N, j * N - 1).t();
      irf.slice(t - 1) += B_j * irf.slice(t - j - 1);
    } // END j loop
  } // END t loop
  
  return irf;
} // END ir1_cpp



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarSIGNs_ir (
    arma::cube&   posterior_B,        // (K, N, S)
    arma::cube&   posterior_Theta0,   // (N, N, S)
    const int     horizon,
    const int     p,
    const bool    standardise = false
) {
  
  const int       N = posterior_B.n_cols;
  const int       S = posterior_B.n_slices;
  
  cube            aux_irfs(N, N, horizon + 1);
  field<cube>     irfs(S);
  
  for (int s=0; s<S; s++) {
    mat irf_0           = posterior_Theta0.slice(s);
    if ( standardise ) {
      irf_0             = irf_0 * diagmat(pow(diagvec(irf_0), -1));
    }
    aux_irfs            = ir1_cpp( posterior_B.slice(s), irf_0, horizon, p );
    irfs(s)             = aux_irfs;
  } // END s loop
  
  return irfs;
} // END bsvarSIGNs_ir
