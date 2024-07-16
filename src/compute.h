
#ifndef _COMPUTE_H_
#define _COMPUTE_H_

#include <RcppArmadillo.h>


arma::cube bsvarSIGNs_structural_shocks (
    const arma::cube&     posterior_B,    // (N, N, S)
    const arma::cube&     posterior_A,    // (N, K, S)
    const arma::mat&      Y,              // NxT dependent variables
    const arma::mat&      X               // KxT dependent variables
);


arma::cube bsvarSIGNs_fitted_values (
    arma::cube&     posterior_A,        // NxKxS
    arma::cube&     posterior_B,        // NxNxS
    arma::cube&     posterior_sigma,    // NxTxS
    arma::mat&      X                   // KxT
);


#endif  // _COMPUTE_H_