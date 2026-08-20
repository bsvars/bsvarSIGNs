#ifndef _SAMPLE_SOE_H_
#define _SAMPLE_SOE_H_

#include <RcppArmadillo.h>

arma::field<arma::mat> sample_restricted_B_cpp(
    const arma::mat& post_B,
    const arma::mat& post_V,
    const arma::mat& Sigma,
    const int& p,
    const int& N,
    const int& Nf,
    const int& K
);

#endif  // _SAMPLE_SOE_H_
