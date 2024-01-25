#ifndef _UTILS_H_
#define _UTILS_H_

#include <RcppArmadillo.h>

double log_det_lower_cpp(const arma::mat& L);

arma::mat qr_sign_cpp(const arma::mat& A);

#endif  // _UTILS_H_
