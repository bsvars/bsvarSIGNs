
#include <RcppArmadillo.h>

#include <bsvars.h>

#include "utils.h"

using namespace Rcpp;
using namespace arma;


// If matches signs
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
bool match_sign(const arma::mat& A, const arma::mat sign) {
  return accu(((A % sign) > 0)) >= accu(abs(sign));
}


// If matches traditional sign restrictions on impulse response functions
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
bool match_sign_irf(const arma::mat&  Q,
                    const arma::cube& irf,
                    const arma::cube& sign_irf) {
  int horizons = irf.n_slices;
  for (int h=0; h<horizons; h++) {
    if ( !match_sign(irf.slice(h)*Q, sign_irf.slice(h)) ) {
      return false;
    }
  }
  return true;
}


// If matches narrative sign restrictions on historical decompositions
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
bool match_sign_hd(const arma::mat& Q,
                   const arma::mat& sign_hd) {
  // TODO
  return true;
}


// Sample rotation matrix Q
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_Q(arma::mat aux_B,
                   arma::mat aux_A,
                   const int& lags,
                   const arma::cube& sign_irf,
                   const arma::mat& sign_hd) {
  int N = aux_B.n_rows;
  int h = sign_irf.n_slices;
  
  int n_tries   = 0;
  int max_tries = pow(10, 6);
  bool success  = false;
  
  mat Q(N, N);
  
  cube irf = bsvars::bsvars_ir1(aux_B, aux_A, h-1, lags);
  
  while (n_tries < max_tries && !success) {
    Q = rortho_cpp(N);
    
    if (match_sign_irf(Q, irf, sign_irf) && match_sign_hd(Q, sign_hd)) {
      success = true;
    }
    
    n_tries++;
  }
  
  if (!success) {
    Rcout << "Warning: could not find a valid Q matrix after " << max_tries << " tries." << endl;
  }
  
  return Q;
}