#ifndef _SAMPLE_HYPER_H_
#define _SAMPLE_HYPER_H_

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


Rcpp::List mn_prior(
    const int&       p,
    const double&    lambda,
    const arma::vec& psi
);


Rcpp::List extend_dummy(
    const int&       p,
    const arma::vec& hyper,
    const arma::mat& Y,
    const arma::mat& X
);


#endif  // _SAMPLE_HYPER_H_
