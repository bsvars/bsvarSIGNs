#ifndef _MCMC_H_
#define _MCMC_H_

#include <RcppArmadillo.h>


arma::mat metropolis(
    const int& T,
    const int& t0,
    arma::vec  x,
    arma::mat  Sigma,
    const std::function<double(const arma::vec&)>& log_target
);


#endif  // _MCMC_H_