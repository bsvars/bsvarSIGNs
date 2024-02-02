#ifndef _SAMPLE_Q_H_
#define _SAMPLE_Q_H_

#include <RcppArmadillo.h>

#include <bsvars.h>

#include "utils.h"

using namespace Rcpp;
using namespace arma;

bool match_sign_cpp(const arma::mat& A, const arma::mat sign);

bool match_sign_irf(const arma::mat&  Q,
                    const arma::cube& irf,
                    const arma::cube& sign_irf);

arma::mat sample_Q(arma::mat aux_B,
                   arma::mat aux_A,
                   const int& lags,
                   const arma::cube& sign_irf,
                   const arma::mat& sign_hd);

#endif  // _SAMPLE_Q_H_
