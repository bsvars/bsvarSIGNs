#ifndef _UTILS_BSVARSIGNS_H_
#define _UTILS_BSVARSIGNS_H_

#include <RcppArmadillo.h>


arma::mat qr_sign_cpp(const arma::mat& A);

arma::mat rortho_cpp(const int& N);

bool match_sign(
    const arma::mat& A, 
    const arma::mat& sign
);

arma::mat Df(
    const  std::function<arma::colvec(const arma::colvec&)>& f,
    const  arma::colvec& x,
    double               h = 1e-10
);

arma::mat metropolis(
    const int& T,
    const int& t0,
    arma::vec  x,
    arma::mat  Sigma,
    const std::function<double(const arma::vec&)>& log_target
);

#endif  // _UTILS_BSVARSIGNS_H_
