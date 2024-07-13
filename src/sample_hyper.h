#ifndef _SAMPLE_HYPER_H_
#define _SAMPLE_HYPER_H_

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


arma::mat sample_hyper(
    const int&        S,
    const int&        start,
    const arma::vec&  init,
    const arma::vec&  model,
    const arma::mat&  Y,
    const arma::mat&  X,
    const arma::mat&  W,
    const Rcpp::List& prior
);

#endif  // _SAMPLE_HYPER_H_
