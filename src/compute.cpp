
#include <RcppArmadillo.h>
#include "progress.hpp"
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



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarSIGNs_hd (
    arma::field<arma::cube>&    posterior_irf_T,    // output of bsvars_irf with irfs at T horizons
    arma::cube&                 structural_shocks,  // NxTxS output bsvars_structural_shocks
    const bool                  show_progress = true
) {
  field<cube> hds = bsvars::bsvars_hd (posterior_irf_T, structural_shocks, show_progress);
  return hds;
} // END bsvarSIGNs_hd



// compute historical decomposition from t to t+h of the i-th variable
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat hd1_cpp(
    const int&        var_i,   // i-th variable
    const int&        t,       // start at period t
    const int&        h,       // number of horizons
    const arma::mat&  Epsilon, // structural shocks, NxT
    const arma::cube& irf
) {
  
  int ii = var_i - 1;
  int tt = t - 1;
  int N  = Epsilon.n_rows;
  
  mat hd(N, h + 1);
  
  // for each shock jj
  for(int jj = 0; jj < N; jj++) {
    // for each horizon hh
    for(int hh = 0; hh <= h; hh++) {
      // for each lag ll
      for(int ll = 0; ll <= hh; ll++) {
        hd(jj, hh) += irf(ii, jj, ll) * Epsilon(jj, tt + hh - ll);
      }
    }
  }
  return hd;
} // END hd1_cpp



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::cube> bsvarSIGNs_fevd (
    arma::field<arma::cube>&    posterior_irf   // output of bsvars_irf
) {
  field<cube>     fevds = bsvars::bsvars_fevd_homosk (posterior_irf);
  return fevds;
} // END bsvarSIGNs_fevd
