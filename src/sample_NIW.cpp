
#include <RcppArmadillo.h>

#include "utils_bsvarsigns.h"

using namespace Rcpp;
using namespace arma;


// sample from matrix normal, U = var between rows, V = var between columns
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat rmatnorm_cpp(const arma::mat& M,
                      const arma::mat& U,
                      const arma::mat& V) {
  
  mat X = mat(size(M), fill::randn);
  return M + chol(U).t() * X * chol(V);
}


// sample from inverse Wishart
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat riwish_cpp (
    const arma::mat&  S, 
    const double&     nu
) {
  // Based on algorithm B.4.4. from Appendinx B by Bauwens, Lubrano, Richard (1999) Bayesian Inference in Dynamic Econometric Models, Oxford Uni Press
  
  int N           = S.n_cols;
  
  mat s_chol      = chol(S, "lower");
  
  mat Q(N, N, fill::zeros);
  Q.diag()        = sqrt(pow(chi2rnd(nu - N + 1, N), -1));
  
  for (int i = 0; i < (N - 1); i++) {
    Q.submat(i + 1, i, N - 1, i) = randn(N - i - 1);
  }
  mat Q_inv       = inv(trimatu(Q));
  
  return s_chol * Q_inv.t() * Q_inv * s_chol.t();
}


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::mat> niw_cpp(
    const arma::mat& Y,
    const arma::mat& X,
    const arma::mat& prior_B,
    const arma::mat& prior_V,
    const arma::mat& prior_S,
    const int&       prior_nu
) {

  const int T  = Y.n_rows;

  // analytic solutions
  mat prior_V_inv = inv_sympd(prior_V);
  mat post_V_inv  = prior_V_inv + X.t() * X;
  mat post_V      = inv_sympd(post_V_inv);
  mat post_B      = post_V * (X.t() * Y + prior_V_inv * prior_B);
  
  // marginal posterior of Sigma
  mat post_S  = prior_S + Y.t() * Y + 
                prior_B.t() * prior_V_inv * prior_B - 
                post_B.t() * post_V_inv * post_B;
  int post_nu = prior_nu + T;
  
  field<mat> post(4);
  post(0) = post_B;
  post(1) = post_V;
  post(2) = symmatu(post_S);
  post(3) = mat(1, 1, fill::ones) * post_nu;
  
  return post;
}

