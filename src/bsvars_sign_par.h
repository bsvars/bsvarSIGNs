#ifndef _BSVARS_SIGN_PAR_H_
#define _BSVARS_SIGN_PAR_H_

#include <RcppArmadillo.h>

Rcpp::List bsvar_sign_par_cpp(
    const int &p,                    // number of lags
    const arma::mat &Y,              // TxN dependent variables
    const arma::mat &X,              // TxK dependent variables
    const arma::cube &sign_irf,      // NxNxh cube of signs for impulse response function
    const arma::mat &sign_narrative, // ANYx6 matrix of signs for historical decomposition
    const arma::mat &sign_B,         // NxN matrix of signs for B
    const arma::field<arma::mat> &Z, // a list of zero restrictions
    const int &Nf,                   // number of foreign variables for SOE
    const Rcpp::List &prior,         // a list of priors
    const int &max_tries = 10000,    // maximum tries for Q draw
    const int idx = 0                // specific hyperparameter column index
);

#endif // _BSVARS_SIGN_PAR_H_