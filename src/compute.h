
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


arma::cube ir1_cpp (
    const arma::mat& B,           // KxN
    const arma::mat& Theta0,      // NxN
    int              horizon,
    const int&       p
);


arma::field<arma::cube> bsvarSIGNs_ir (
    arma::cube&   posterior_B,        // (K, N, S)
    arma::cube&   posterior_Theta0,   // (N, N, S)
    const int     horizon,
    const int     p,
    const bool    standardise = false
);


#endif  // _COMPUTE_H_