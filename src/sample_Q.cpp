
#include <RcppArmadillo.h>

#include "utils.h"
#include "restrictions_narrative.h"
#include "restrictions_zero.h"
#include "compute.h"

using namespace Rcpp;
using namespace arma;


// If matches traditional sign restrictions on impulse response functions
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
bool match_sign_irf(
    const arma::mat&  Q,
    const arma::cube& sign_irf,
    const arma::cube& irf
) {
  
  int h = sign_irf.n_slices;
  for (int t=0; t<h; t++) {
    if ( !match_sign(irf.slice(t) * Q, sign_irf.slice(t)) ) {
      return false;
    }
  }
  return true;
}


// Sample rotation matrix Q and updates importance weight by reference
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::mat> sample_Q(    
    const int&                    p,
    const arma::mat&              Y,
    const arma::mat&              X,
    arma::mat&                    B,
    arma::mat&                    h_invp,
    arma::mat&                    chol_Sigma,
    const Rcpp::List&             prior,
    const arma::cube&             sign_irf,
    const arma::mat&              sign_narrative,
    const arma::mat&              sign_B,
    const arma::field<arma::mat>& Z,
    const int&                    max_tries
) {
  
  const int N          = Y.n_cols;
  const int T          = Y.n_rows;
  
  int    n_tries       = 0;
  int    h             = sign_narrative.col(5).max(); // maximum horizon
  if (h < sign_irf.n_slices - 1) {h = sign_irf.n_slices - 1;}
  
  bool   has_narrative = sign_narrative(0, 0) != 0;
  bool   has_zero      = Z.n_elem > 0;
  
  mat    Q(N, N);
  mat    lt_shocks     = h_invp * (Y - X * B).t();  // lower triangular identified shocks
  
  cube   irf           = ir1_cpp(B, chol_Sigma, h, p);  // reduced-form irf
  
  bool   success       = false;
  mat    shocks;
  while (n_tries < max_tries && !success) {
    if (has_zero) {
      Q = rzeroQ(Z, irf.slice(0));
    } else {
      Q = rortho_cpp(N);
    }
    
    shocks = Q.t() * lt_shocks;
    
    if (match_sign_irf(Q, sign_irf, irf) && match_sign(Q.t() * h_invp, sign_B)) {
      if (!has_narrative) {
        success = true;
      } else {
        success = match_sign_narrative(shocks, sign_narrative, irf);
      }
    }
    
    n_tries++;
  }
  
  double w = 1;
  if (!success) {
    w = 0;
  } else {
    if (has_narrative) {
      w *= weight_narrative(T, sign_narrative, irf);  
    }
    if (has_zero) {
      w *= weight_zero(Z, B, h_invp, Q);
    }
  }
  
  field<mat> result(3);
  result(0) = Q;
  result(1) = shocks;
  result(2) = w;
  
  return result;
}

