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


Rcpp::List update_prior(
    const int&        p,
    const arma::vec&  hyper,
    const Rcpp::List& prior
);


Rcpp::List extend_dummy(
    const int&       p,
    const arma::vec& hyper,
    const arma::mat& Y,
    const arma::mat& X
);


arma::mat sample_hyper(
    const int&        S,
    const int&        p,
    const double&     c,
    const arma::mat&  Y,
    const arma::mat&  X,
    const Rcpp::List& prior
);


#endif  // _SAMPLE_HYPER_H_
