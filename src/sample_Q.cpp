
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
    const arma::mat sign
) {
  
  return accu(((A % sign) > 0)) == accu(sign != 0);
  
}


// If matches traditional sign restrictions on impulse response functions
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
bool match_sign_irf(
    const arma::mat&  Q,
    const arma::cube& irf,
    const arma::cube& sign_irf
) {
  
  int h = sign_irf.n_slices;
  for (int t=0; t<h; t++) {
    if ( !match_sign(irf.slice(t)*Q, sign_irf.slice(t)) ) {
      return false;
    }
  }
  return true;
}


// compute historical decomposition once
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat hd1(
    const int&        i,  // i-th variable
    const int&        t,  // start at period t
    const int&        h,  // number of horizons
    const arma::mat&  U,  // structural shocks 
    const arma::cube& irf
) {
  
   int ii = i - 1;
   int tt = t - 1;
   int N  = U.n_rows;
   
   mat hd = zeros(N, h);
   
   // for each shock jj
   for(int jj = 0; jj < N; jj++) {
     // for each horizon hh
     for(int hh = 0; hh < h; hh++) {
       // for each lag ll
       for(int ll = 0; ll <= hh; ll++) {
         hd(jj, hh) += irf(ii, jj, ll) * U(jj, tt + hh - ll);
       }
     }
   }
   return hd;
 }


// If matches narrative sign restrictions on historical decomposition
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
arma::mat sample_Q(    
    const int&        lags,
    const arma::mat&  Y,
    const arma::mat&  X,
    arma::mat         aux_A,
    arma::mat         aux_B,
    const arma::cube& sign_irf,
    const arma::mat&  sign_hd,
    const arma::mat&  sign_B
) {
  
  int N          = aux_B.n_rows;
  int n_tries    = 0;
  int max_tries  = pow(10, 6);
  int h          = sign_hd.col(5).max(); // maximum horizon
  if (h < sign_irf.n_slices-1) {
    h = sign_irf.n_slices-1;
  }
  
  bool success   = false;
  
  mat Q(N, N);
  
  cube irf       = bsvars::bsvars_ir1(aux_B, aux_A, h, lags);
  
  while (n_tries < max_tries && !success) {
    Q = rortho_cpp(N);
    
    if (match_sign_hd(Q, sign_hd) &&
        match_sign_irf(Q, irf, sign_irf) &&
        match_sign(aux_B, sign_B)) {
      success = true;
    }
    
    n_tries++;
  }
  
  if (!success) {
    Rcout << "Warning: could not find a valid Q matrix after " << max_tries << " tries." << endl;
    Q.zeros(); // return a matrix of zeros if no valid Q matrix is found
  }
  
  return Q;
}