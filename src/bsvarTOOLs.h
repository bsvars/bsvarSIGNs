#ifndef _BSVARTOOLS_H_
#define _BSVARTOOLS_H_

#include <RcppArmadillo.h>

using namespace arma;

using namespace Rcpp;
using namespace arma;


arma::mat hd1_cpp(
    const int&        var_i,
    const int&        t,
    const int&        h,
    const arma::mat&  Epsilon,
    const arma::cube& irf
);
  
  
arma::cube ir1_cpp(
    const arma::mat& B,           // KxN
    const arma::mat& Theta0,      // NxN
    int              horizon,
    const int&       p
);


arma::field<arma::cube> bsvarSIGNs_ir (
    arma::cube&   posterior_B,        // (K, N, S)
    arma::cube&   posterior_Theta0,   // (N, N, S)
    const int     horizon,
    const int     p
);


#endif  // _BSVARTOOLS_H_