
#include <RcppArmadillo.h>

#include "utils.h"

using namespace Rcpp;
using namespace arma;


// sample from matrix normal, U = var between rows, V = var between columns
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
arma::mat rmatnorm_cpp(const arma::mat& M,
                      const arma::mat& U,
                      const arma::mat& V) {
  
  mat X = mat(size(M), fill::randn);
  return M + chol(U).t() * X * chol(V);
}


// sample from inverse Wishart
// [[Rcpp:interface(cpp)]]
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


// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
void niw_cpp(
    arma::mat&       post_A,
    arma::mat&       post_V,
    arma::mat&       post_S,
    int&             post_nu,
    const arma::mat& Y,
    const arma::mat& X,
    const Rcpp::List prior
) {

  const int T  = Y.n_cols;
  
  mat prior_A  = as<mat>(prior["B"]);
  mat prior_V  = as<mat>(prior["V"]);
  mat prior_S  = as<mat>(prior["S"]);
  int prior_nu = as<int>(prior["nu"]);
  
  prior_A = prior_A.t();
  
  // analytic solutions
  mat prior_V_inv = inv_sympd(prior_V);
  mat post_V_inv  = prior_V_inv + X * X.t();
  post_V      = inv_sympd(post_V_inv);
  post_A      = (prior_A * prior_V_inv + Y * X.t()) * post_V;
  
  // marginal posterior of Sigma
  post_S  = prior_S + Y * Y.t() + prior_A * prior_V_inv * prior_A.t() - post_A * post_V_inv * post_A.t();
  post_S      = symmatu(post_S);
  post_nu = prior_nu + T;
}

