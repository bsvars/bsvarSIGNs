#ifndef _SAMPLE_NIW_H_
#define _SAMPLE_NIW_H_

#include <RcppArmadillo.h>

#include "utils.h"

using namespace Rcpp;
using namespace arma;

void niw_cpp(arma::mat&       post_A,
             arma::mat&       post_V,
             arma::mat&       post_S,
             int&             post_nu,
             const arma::mat& Y,
             const arma::mat& X,
             const Rcpp::List prior);

#endif  // _SAMPLE_NIW_H_
