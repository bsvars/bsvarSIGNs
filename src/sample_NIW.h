#ifndef _SAMPLE_NIW_H_
#define _SAMPLE_NIW_H_

#include <RcppArmadillo.h>

arma::mat rmatnorm_cpp(
    const arma::mat& M,
    const arma::mat& U,
    const arma::mat& V
);

arma::mat riwish_cpp (
    const arma::mat&  S, 
    const double&     nu
);

arma::field<arma::mat> niw_cpp(
    const arma::mat& Y,
    const arma::mat& X,
    const arma::mat& prior_B,
    const arma::mat& prior_V,
    const arma::mat& prior_S,
    const int&       prior_nu
);

#endif  // _SAMPLE_NIW_H_
