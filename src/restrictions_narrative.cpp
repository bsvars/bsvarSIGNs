
#include <RcppArmadillo.h>

#include "utils.h"
#include "bsvarTOOLs.h"

using namespace Rcpp;
using namespace arma;


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


// approximate importance weight for narrative restrictions
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
double weight_narrative(
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
