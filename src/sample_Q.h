#ifndef _SAMPLE_Q_H_
#define _SAMPLE_Q_H_

#include <RcppArmadillo.h>

#include <bsvars.h>

#include "utils.h"

using namespace Rcpp;
using namespace arma;

arma::mat sample_Q(    
    const int&                    lags,
    const arma::mat&              Y,
    const arma::mat&              X,
    double&                       aux_w,
    arma::mat&                    aux_A,
    arma::mat&                    aux_B,
    arma::mat&                    chol_SIGMA,
    const Rcpp::List&             prior,
    const arma::field<arma::mat>& VB,
    const arma::cube&             sign_irf,
    const arma::mat&              sign_narrative,
    const arma::mat&              sign_B,
    const arma::cube&             zero_irf,
    const int&                    max_tries,
    bool&                         success
);

#endif  // _SAMPLE_Q_H_
