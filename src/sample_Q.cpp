
#include <RcppArmadillo.h>

#include <bsvars.h>

#include "utils.h"

using namespace Rcpp;
using namespace arma;


// If matches signs
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
bool match_sign(
    const arma::mat& A, 
    const arma::mat& sign
) {
  return accu(((A % sign) > 0)) == accu(sign != 0);
}


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


// compute historical decomposition from t to t+h of the i-th variable
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat hd1(
    const int&        var_i,  // i-th variable
    const int&        t,      // start at period t
    const int&        h,      // number of horizons
    const arma::mat&  U,      // structural shocks 
    const arma::cube& irf
) {
  
   int ii = var_i - 1;
   int tt = t - 1;
   int N  = U.n_rows;
   
   mat hd(N, h + 1);
   
   // for each shock jj
   for(int jj = 0; jj < N; jj++) {
     // for each horizon hh
     for(int hh = 0; hh <= h; hh++) {
       // for each lag ll
       for(int ll = 0; ll <= hh; ll++) {
         hd(jj, hh) += irf(ii, jj, ll) * U(jj, tt + hh - ll);
       }
     }
   }
   return hd;
 }


// If matches narrative sign restrictions
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
bool match_sign_narrative(
    const arma::mat&  U,
    const arma::mat&  sign_narrative,
    const arma::cube& irf
) {

  bool is_greater;
  
  int type, var_i, shock_j, t, h;
  
  mat one, u, hd, hd_j;
  
  for (int k=0; k<sign_narrative.n_rows; k++) {
    
    // sign_narrative: | type | sign | var_i | shock_j | start | horizons |
    type        = sign_narrative(k, 0);
    is_greater  = sign_narrative(k, 1) == 0;
    var_i       = sign_narrative(k, 2);
    shock_j     = sign_narrative(k, 3) - 1;
    t           = sign_narrative(k, 4);
    h           = sign_narrative(k, 5);
    
    one         = ones(1, h + 1);
    if (!is_greater) {
      // if restriction is smaller than zero
      one       = -one;
    }
    
    if (type == 1) {
      // narrative sign restrictions on structural shocks
      u = U.submat(shock_j, t - 1, shock_j, t + h - 1);
      
      if (!match_sign(u, one)) {
        return false;
      }
    } else {
      // narrative sign restrictions on historical decomposition
      hd   = abs(hd1(var_i, t, h, U, irf));
      hd_j = hd.row(shock_j);
      hd.shed_row(shock_j);
      
      if (type == 2) {
        // type A
        if ( is_greater && !match_sign(hd_j - max(hd), one)) {
          return false;
        }
        if (!is_greater && !match_sign(hd_j - min(hd), one)) {
          return false;
        }
      } else if (          !match_sign(hd_j - sum(hd), one)) {
        // type B
        return false;
      }
    }
  }
  return true;
}


// approximate importance weight
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
double approximate_w(
    const int&                    T,
    arma::mat                     sign_narrative,
    const arma::cube&             irf
) {
  
  const int M         = 1e+04;  // number of draws to approximate normal distribution
  
  double    n_success = 1.0e-15;
  
  cube      Z(irf.n_rows, sign_narrative.col(5).max() + 1, M, fill::randn);
  // cube      Z(irf.n_rows, T, M, fill::randn);
  
  // change all starting period to the first period
  // since we use the same M draws for all narrative restrictions
  sign_narrative.col(4) = ones(sign_narrative.n_rows, 1);
  
  for (int m=0; m<M; m++) {
    if (match_sign_narrative(Z.slice(m), sign_narrative, irf)) {
      n_success++;
    }
  }
  return M / n_success;
}


// Sample rotation matrix Q and updates importance weight by reference
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat sample_Q(    
    const int&                    lags,
    const arma::mat&              Y,
    const arma::mat&              X,
    double&                       aux_w,
    arma::mat&                    aux_A,
    arma::mat&                    aux_B,
    arma::mat&                    aux_hyper,
    const Rcpp::List&             prior,
    const arma::field<arma::mat>& VB,
    const arma::cube&             sign_irf,
    const arma::mat&              sign_narrative,
    const arma::mat&              sign_B,
    const int&                    max_tries,
    int&                          n_fails
) {
  
  const int N          = Y.n_rows;
  const int T          = Y.n_cols;
  
  int    n_tries       = 0;
  int    h             = sign_narrative.col(5).max(); // maximum horizon
  if (h < sign_irf.n_slices - 1) {h = sign_irf.n_slices - 1;}
  
  bool   success       = false;
  bool   has_narrative = sign_narrative(0, 0) != 0;
  
  mat    Q(N, N);
  mat    U             = aux_B * (Y - aux_A * X);
  
  cube   irf           = bsvars::bsvars_ir1(aux_B, aux_A, h, lags);
  
  while (n_tries < max_tries && !success) {
    Q = rortho_cpp(N);
    
    if (match_sign_irf(Q, sign_irf, irf) && match_sign(aux_B, sign_B)) {
      if (!has_narrative) {
        success = true;
      } else {
        success = match_sign_narrative(Q.t() * U, sign_narrative, irf);
      }
    }
    
    n_tries++;
  }
  
  if (!success) {
    n_fails++;
    
    aux_A         = bsvars::sample_A_homosk1(aux_A, aux_B, aux_hyper, Y, X, prior);
    aux_B         = bsvars::sample_B_homosk1(aux_B, aux_A, aux_hyper, Y, X, prior, VB);
    
    return sample_Q(lags, Y, X,
                    aux_w, aux_A, aux_B, aux_hyper, 
                    prior, VB,
                    sign_irf, sign_narrative, sign_B,
                    max_tries, n_fails);
  }
  
  if (has_narrative) {
    aux_w = approximate_w(T, sign_narrative, irf);
  }
  
  return Q;
}