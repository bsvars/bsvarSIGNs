#ifndef _BSVARTOOLS_H_
#define _BSVARTOOLS_H_

#include <RcppArmadillo.h>

using namespace arma;

using namespace Rcpp;
using namespace arma;

arma::cube ir1_cpp(
    const arma::mat& A, 
    const arma::mat& B, 
    const int& periods,
    const int& p
);

#endif  // _BSVARTOOLS_H_