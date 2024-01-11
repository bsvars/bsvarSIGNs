#ifndef _UTILS_H_
#define _UTILS_H_

#include <RcppArmadillo.h>

arma::mat matnrnd_cpp(const arma::mat& M,
                      const arma::mat& U,
                      const arma::mat& V);

arma::mat mvnrnd_inverse_cpp(const arma::mat& mu, const arma::mat& inv_Sigma);

double log_det_lower_cpp(const arma::mat& L);

arma::mat inv_chol_cpp(const arma::mat& L);

double log_mvnpdf_cpp(const arma::mat& x,
                      const arma::mat& mu,
                      const arma::mat& inv_Sigma,
                      const arma::mat& L);

arma::mat qr_sign_cpp(const arma::mat& A);

#endif  // _UTILS_H_
