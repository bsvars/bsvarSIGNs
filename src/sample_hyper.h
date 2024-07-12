#ifndef _SAMPLE_HYPER_H_
#define _SAMPLE_HYPER_H_

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


Rcpp::List mn_prior(
    const int&       p,
    const double&    lambda,
    const arma::vec& psi
);


Rcpp::List extend_dummy(
    const int&       p,
    const arma::vec& hyper,
    const arma::mat& Y,
    const arma::mat& X
);


arma::mat sample_hyper(
    const int&        S,
    const int&        start,
    const int&        p,
    const arma::vec&  init,
    const arma::vec&  model,
    const arma::mat&  Y,
    const arma::mat&  X,
    const arma::mat&  W
);

#endif  // _SAMPLE_HYPER_H_
