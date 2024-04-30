#ifndef _SAMPLE_Q_H_
#define _SAMPLE_Q_H_

#include <RcppArmadillo.h>

#include "utils.h"

using namespace Rcpp;
using namespace arma;

Rcpp::List sample_Q(    
    const int&                    p,
    const arma::mat&              Y,
    const arma::mat&              X,
    arma::mat&                    B,
    arma::mat&                    h_invp,
    arma::mat&                    chol_Sigma,
    const Rcpp::List&             prior,
    const arma::field<arma::mat>& VB,
    const arma::cube&             sign_irf,
    const arma::mat&              sign_narrative,
    const arma::mat&              sign_B,
    const arma::field<arma::mat>& Z,
    const int&                    max_tries
);

#endif  // _SAMPLE_Q_H_
