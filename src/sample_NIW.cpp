
#include <RcppArmadillo.h>

#include "utils.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
void niw_cpp(
    arma::cube&      A,
    arma::cube&      SIGMA,
    const arma::mat& Y,       // (T, N)
    const arma::mat& X,       // (T, K)
    const int&       S,
    const arma::mat& prior_A,
    const arma::mat& prior_V,
    const arma::mat& prior_S,
    const int&       prior_nu
) {

  int T = Y.n_rows;
  int N = Y.n_cols;
  int K = X.n_cols;
  
  // analytic solutions
  mat prior_V_inv = inv_sympd(prior_V);
  mat post_V_inv  = prior_V_inv + X.t() * X;
  mat post_V      = inv_sympd(post_V_inv);
  mat post_A      = post_V * (prior_V_inv * prior_A + X.t() * Y);
  
  // marginal posterior of Sigma
  mat post_S  = prior_S + Y.t() * Y + prior_A.t() * prior_V_inv * prior_A - post_A.t() * post_V_inv * post_A;
  post_S      = symmatu(post_S);
  int post_nu = prior_nu + T;
  
  mat a     = zeros(K, N);
  mat sigma = zeros(N, N);
  
  for(int s = 0; s < S; s++) {
    sigma          = iwishrnd(post_S, post_nu);
    a              = matnrnd_cpp(post_A, post_V, sigma);
    SIGMA.slice(s) = sigma;
    A.slice(s)     = a.t();
  }
}

