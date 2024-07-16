
#ifndef _COMPUTE_H_
#define _COMPUTE_H_

#include <RcppArmadillo.h>


arma::cube structural_shocks (
    const arma::cube&     posterior_B,    // (N, N, S)
    const arma::cube&     posterior_A,    // (N, K, S)
    const arma::mat&      Y,              // NxT dependent variables
    const arma::mat&      X               // KxT dependent variables
);


#endif  // _COMPUTE_H_