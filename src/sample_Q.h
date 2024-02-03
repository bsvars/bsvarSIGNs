#ifndef _SAMPLE_Q_H_
#define _SAMPLE_Q_H_

#include <RcppArmadillo.h>

#include <bsvars.h>

#include "utils.h"

using namespace Rcpp;
using namespace arma;

arma::mat sample_Q(    
    const int&        lags,               // number of lags
    const arma::mat&  Y,                  // NxT dependent variables
    const arma::mat&  X,                  // KxT dependent variables
    arma::mat         aux_A,
    arma::mat         aux_B,
    const arma::cube& sign_irf,
    const arma::mat&  sign_hd,
    const arma::mat&  sign_B
);

#endif  // _SAMPLE_Q_H_
