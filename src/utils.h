#ifndef _UTILS_H_
#define _UTILS_H_

#include <RcppArmadillo.h>


arma::mat qr_sign_cpp(const arma::mat& A);

arma::mat rortho_cpp(const int& N);

bool match_sign(
    const arma::mat& A, 
    const arma::mat& sign
);

arma::mat Df(
    const std::function<arma::vec(const arma::vec&)>& f,
    const arma::vec& x
);

arma::mat metropolis(
    const int& T,
    const int& t0,
    arma::vec  x,
    arma::mat  Sigma,
    const std::function<double(const arma::vec&)>& log_target
);

#endif  // _UTILS_H_
