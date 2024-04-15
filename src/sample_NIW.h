#ifndef _SAMPLE_NIW_H_
#define _SAMPLE_NIW_H_

#include <RcppArmadillo.h>

#include "utils.h"

using namespace Rcpp;
using namespace arma;


arma::mat rmnorm_cpp(const arma::mat& M,
                      const arma::mat& U,
                      const arma::mat& V);

arma::mat riwish_cpp (
    const arma::mat&  S, 
    const double&     nu
);

void niw_cpp(arma::mat&       post_A,
             arma::mat&       post_V,
             arma::mat&       post_S,
             int&             post_nu,
             const arma::mat& Y,
             const arma::mat& X,
             const Rcpp::List prior);

#endif  // _SAMPLE_NIW_H_
