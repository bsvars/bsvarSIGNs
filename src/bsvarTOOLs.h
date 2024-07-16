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
  

#endif  // _BSVARTOOLS_H_