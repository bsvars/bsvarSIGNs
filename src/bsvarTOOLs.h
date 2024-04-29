#ifndef _BSVARTOOLS_H_
#define _BSVARTOOLS_H_

#include <RcppArmadillo.h>

using namespace arma;

using namespace Rcpp;
using namespace arma;


arma::mat hd1_cpp(
    const int&        var_i,  // i-th variable
    const int&        t,      // start at period t
    const int&        h,      // number of horizons
    const arma::mat&  U,      // structural shocks 
    const arma::cube& irf
);
  
  
arma::cube ir1_cpp(
    const arma::mat& B,           // KxN
    const arma::mat& Theta0,      // NxN
    int              horizon,
    const int&       p
);

#endif  // _BSVARTOOLS_H_