#ifndef _UTILS_H_
#define _UTILS_H_

#include <RcppArmadillo.h>

arma::mat qr_sign_cpp(const arma::mat& A);

arma::mat rortho_cpp(const int& N);

#endif  // _UTILS_H_
