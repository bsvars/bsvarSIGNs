#ifndef _BSVARS_SIGN_H_
#define _BSVARS_SIGN_H_


#include <RcppArmadillo.h>

Rcpp::List bsvar_sign_cpp(
    const int&  S,                        // number of draws from the posterior
    const int&  lags,                     // number of lags
    const arma::mat&  Y,                  // NxT dependent variables
    const arma::mat&  X,                  // KxT dependent variables
    const arma::field<arma::mat>& VB,     // N-list
    const arma::cube& sign_irf,           // NxNxh cube of signs for impulse response function
    const arma::mat&  sign_hd,            // Mx6 matrix of signs for historical decomposition
    const Rcpp::List& prior,              // a list of priors
    const Rcpp::List& starting_values,    // a list of starting values
    const int         thin = 100,         // introduce thinning
    const bool        show_progress = true
);

#endif  // _BSVARS_SIGN_H_