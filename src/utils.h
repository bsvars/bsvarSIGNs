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
    const  std::function<arma::colvec(const arma::colvec&)>& f,
    const  arma::colvec& x,
    double               h = 1e-10
);

#endif  // _UTILS_H_
