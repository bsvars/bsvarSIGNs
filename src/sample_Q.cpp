
#include <RcppArmadillo.h>

#include "utils.h"
#include "restrictions_narrative.h"
#include "restrictions_zero.h"
#include "bsvarTOOLs.h"

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
    if ( !match_sign(irf.slice(t)*Q, sign_irf.slice(t)) ) {
      return false;
    }
  }
  return true;
}


// Sample rotation matrix Q and updates importance weight by reference
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_Q(    
    const int&                    p,
    const arma::mat&              Y,
    const arma::mat&              X,
    double&                       aux_w,
    arma::mat&                    aux_A,
    arma::mat&                    aux_B,
    arma::mat&                    chol_SIGMA,
    const Rcpp::List&             prior,
    const arma::field<arma::mat>& VB,
    const arma::cube&             sign_irf,
    const arma::mat&              sign_narrative,
    const arma::mat&              sign_B,
    const arma::field<arma::mat>& Z,
    const int&                    max_tries,
    bool&                         success
) {

  const int N          = Y.n_cols;
  const int T          = Y.n_rows;
  
  int    n_tries       = 0;
  int    h             = sign_narrative.col(5).max(); // maximum horizon
  if (h < sign_irf.n_slices - 1) {h = sign_irf.n_slices - 1;}
  
  bool   has_narrative = sign_narrative(0, 0) != 0;
  bool   has_zero      = Z.n_elem > 0;
  
  mat    Q(N, N);
  mat    U             = aux_B * (Y.t() - aux_A.t() * X.t());
  
  cube   irf           = ir1_cpp(aux_A, chol_SIGMA, h, p);
  
  while (n_tries < max_tries && !success) {
    if (has_zero) {
      Q = rzeroQ(Z, irf.slice(0));
    } else {
      Q = rortho_cpp(N);
    }
    
    if (match_sign_irf(Q, sign_irf, irf) && match_sign(Q.t() * aux_B, sign_B)) {
      if (!has_narrative) {
        success = true;
      } else {
        success = match_sign_narrative(Q.t() * U, sign_narrative, irf);
      }
    }
    
    n_tries++;
  }
  
  aux_w = 1;
  if (!success) {
    aux_w = 0;
  } else {
    if (has_narrative) {
      aux_w *= weight_narrative(T, sign_narrative, irf);  
    }
    if (has_zero) {
      aux_w *= weight_zero(Z, aux_A, aux_B, Q);
    }
  }
  
  return Q;
}