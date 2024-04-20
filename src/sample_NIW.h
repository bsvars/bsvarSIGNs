#ifndef _SAMPLE_NIW_H_
#define _SAMPLE_NIW_H_

#include <RcppArmadillo.h>

#include "utils.h"

using namespace Rcpp;
using namespace arma;


arma::mat rmatnorm_cpp(const arma::mat& M,
                       const arma::mat& U,
                       const arma::mat& V);

arma::mat riwish_cpp (
    const arma::mat&  S, 
    const double&     nu
);

Rcpp::List niw_cpp(
    const arma::mat& Y,
    const arma::mat& X,
    const Rcpp::List prior
);

#endif  // _SAMPLE_NIW_H_
