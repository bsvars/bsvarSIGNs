#ifndef _SAMPLE_NIW_H_
#define _SAMPLE_NIW_H_

#include <RcppArmadillo.h>

#include "utils.h"

using namespace Rcpp;
using namespace arma;

void niw_cpp(
                arma::cube&      A,
                arma::cube&      SIGMA,
                const arma::mat& Y,
                const arma::mat& X,
                const int&       S,
                const Rcpp::List prior
);

#endif  // _SAMPLE_NIW_H_
