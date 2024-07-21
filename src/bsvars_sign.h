#ifndef _BSVARS_SIGN_H_
#define _BSVARS_SIGN_H_

#include <RcppArmadillo.h>

Rcpp::List bsvar_sign_cpp(
        const int&        S,                  // number of draws from the posterior
        const int&        p,                  // number of lags
        const arma::mat&  Y,                  // TxN dependent variables
        const arma::mat&  X,                  // TxK dependent variables
        const arma::cube& sign_irf,           // NxNxh cube of signs for impulse response function
        const arma::mat&  sign_narrative,     // ANYx6 matrix of signs for historical decomposition
        const arma::mat&  sign_B,             // NxN matrix of signs for B
        const arma::field<arma::mat>& Z,      // a list of zero restrictions
        const Rcpp::List& prior,              // a list of priors
        const bool        show_progress,
        const bool        parallel,
        const int&        max_tries           // maximum tries for Q draw
);

#endif  // _BSVARS_SIGN_H_