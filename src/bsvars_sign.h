#ifndef _BSVARS_SIGN_H_
#define _BSVARS_SIGN_H_

#include <RcppArmadillo.h>

arma::field<arma::mat> bsvar_sign_single_draw_cpp(
    const int &p,
    const arma::mat &Y,
    const arma::mat &X,
    const arma::cube &sign_irf,
    const arma::mat &sign_narrative,
    const arma::mat &sign_B,
    const arma::field<arma::mat> &Z,
    const int &Nf,
    const Rcpp::List &prior,
    const int &max_tries,
    const int idx
);

Rcpp::List bsvar_sign_par_cpp(
    const int &p,
    const arma::mat &Y,
    const arma::mat &X,
    const arma::cube &sign_irf,
    const arma::mat &sign_narrative,
    const arma::mat &sign_B,
    const arma::field<arma::mat> &Z,
    const int &Nf,
    const Rcpp::List &prior,
    const int &max_tries = 10000,
    const int idx = 0
);

Rcpp::List bsvar_sign_cpp(
        const int&        S,                  // number of draws from the posterior
        const int&        p,                  // number of lags
        const arma::mat&  Y,                  // TxN dependent variables
        const arma::mat&  X,                  // TxK dependent variables
        const arma::cube& sign_irf,           // NxNxh cube of signs for impulse response function
        const arma::mat&  sign_narrative,     // ANYx6 matrix of signs for historical decomposition
        const arma::mat&  sign_B,             // NxN matrix of signs for B
        const arma::field<arma::mat>& Z,      // a list of zero restrictions
        const int&        Nf,                 // number of foreign variables for SOE
        const Rcpp::List& prior,              // a list of priors
        const bool        show_progress = true,
        const int         thin = 100,         // introduce thinning
        const int&        max_tries = 10000   // maximum tries for Q draw
);

#endif  // _BSVARS_SIGN_H_