#ifndef _MCMC_H_
#define _MCMC_H_


#include <functional>
#include <iostream>
#include <RcppArmadillo.h>
#include "Rcpp/Rmath.h"
#include "progress.hpp"

using namespace arma;


arma::mat metropolis(
    const int& T,
    const int& t0,
    arma::vec  x,
    arma::mat  Sigma,
    const std::function<double(const arma::vec&)>& log_target
);


#endif  // _MCMC_H_