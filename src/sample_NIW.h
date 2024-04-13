#ifndef _SAMPLE_NIW_H_
#define _SAMPLE_NIW_H_

#include <RcppArmadillo.h>

#include "utils.h"

using namespace Rcpp;
using namespace arma;

void niw_cpp(
        arma::cube&      A,
        arma::cube&      SIGMA,
        const arma::mat& Y,       // (T, N)
        const arma::mat& X,       // (T, K)
        const int&       S,
        const arma::mat& prior_A,
        const arma::mat& prior_V,
        const arma::mat& prior_S,
        const int&       prior_nu
);

#endif  // _SAMPLE_NIW_H_
