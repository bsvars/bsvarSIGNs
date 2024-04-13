#ifndef _UTILS_H_
#define _UTILS_H_

#include <RcppArmadillo.h>


arma::mat matnrnd_cpp(const arma::mat& M,
                      const arma::mat& U,
                      const arma::mat& V);

arma::mat qr_sign_cpp(const arma::mat& A);

arma::mat rortho_cpp(const int& N);

arma::field<arma::cube> irf_cpp(
    arma::cube&         posterior_B,        // (N, N, S)
    arma::cube&         posterior_A,        // (N, K, S)
    const int           horizon,
    const int           p
);

#endif  // _UTILS_H_
