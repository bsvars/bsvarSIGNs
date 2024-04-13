
#include <RcppArmadillo.h>

#include "utils.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
void niw_cpp(
    arma::cube&      A,
    arma::cube&      SIGMA,
    const arma::mat& Y,
    const arma::mat& X,
    const int&       S,
    const Rcpp::List prior
) {

  const int T  = Y.n_cols;
  
  mat prior_A  = as<mat>(prior["A"]);
  mat prior_V  = as<mat>(prior["V"]);
  mat prior_S  = as<mat>(prior["S"]);
  int prior_nu = as<int>(prior["nu"]);
  
  // analytic solutions
  mat prior_V_inv = inv_sympd(prior_V);
  mat post_V_inv  = prior_V_inv + X * X.t();
  mat post_V      = inv_sympd(post_V_inv);
  mat post_A      = (prior_A * prior_V_inv + Y * X.t()) * post_V;
  
  // marginal posterior of Sigma
  mat post_S  = prior_S + Y * Y.t() + prior_A * prior_V_inv * prior_A.t() - post_A * post_V_inv * post_A.t();
  post_S      = symmatu(post_S);
  int post_nu = prior_nu + T;
  
  mat sigma, a;
  for(int s = 0; s < S; s++) {
    sigma          = iwishrnd(post_S, post_nu);
    a              = matnrnd_cpp(post_A, sigma, post_V);
    SIGMA.slice(s) = sigma;
    A.slice(s)     = a;
  }
}

