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
    arma::mat&                    aux_A,
    arma::mat&                    aux_B,
    arma::mat&                    aux_hyper,
    const Rcpp::List&             prior,              // a list of priors
    const arma::field<arma::mat>& VB,     // N-list
    const arma::cube&             sign_irf,
    const arma::mat&              sign_narrative,
    const arma::mat&              sign_B,
    const arma::cube&             Z,
    double&                       aux_w
);

#endif  // _SAMPLE_Q_H_
