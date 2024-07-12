
#include <RcppArmadillo.h>

#include "utils.h"

using namespace Rcpp;
using namespace arma;


// get prior from hyper-parameters
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
Rcpp::List mn_prior(
  const int&       p,
  const double&    lambda,
  const arma::vec& psi
) {
  
  int N = psi.n_elem;
  int K = 1 + N * p;
  
  mat B(K, N, fill::eye);
  
  vec v(K);
  v.rows(0, K - 2) = kron(pow(linspace(1, p, p), -2), 1 / psi);
  mat V            = lambda * lambda * diagmat(v);
  V(K - 1, K - 1)  = 1e+6;
  
  mat S  = diagmat(psi);
  int nu = N + 2;
  
  return List::create(
    _["B"]   = B,
    _["V"]   = V,
    _["S"]   = S,
    _["nu"]  = nu
  );
}


// extend data with dummy observations
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
Rcpp::List extend_dummy(
    const int&       p,
    const arma::vec& hyper,
    const arma::vec& model,
    const arma::mat& Y,
    const arma::mat& X
) {
  
  int N        = Y.n_cols;
  int K        = X.n_cols;
  
  mat ybar0    = mean(X.submat(0, 0, p - 1, N - 1), 0);
  
  double mu, delta;
  mat Ystar(0, N), Xstar(0, K);
  
  if (model(0)) {
    mu         = hyper(0);
    mat yp     = diagmat(ybar0 / mu);
    mat xp     = join_horiz(repmat(yp, 1, p), zeros(N, 1));
    
    Ystar      = join_vert(Ystar, yp);
    Xstar      = join_vert(Xstar, xp);
  }
  
  if (model(1)) {
    delta      = hyper(1);
    mat ypp    = ybar0 / delta;
    mat xpp    = join_horiz(repmat(ypp, 1, p), mat(1, 1, fill::value(1 / delta)));
    
    Ystar      = join_vert(Ystar, ypp);
    Xstar      = join_vert(Xstar, xpp);
  }
  
  mat Yplus    = join_vert(Ystar, Y);
  mat Xplus    = join_vert(Xstar, X);
  
  return List::create(
    _["Yplus"] = Yplus,
    _["Xplus"] = Xplus,
    _["Ystar"] = Ystar,
    _["Xstar"] = Xstar
  );
}


// log density of gamma distribution
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_dgamma(
    const double& x,
    const double& k,
    const double& theta
) {
  
  return (k - 1) * log(x) - x / theta - k * log(theta) - lgamma(k);
}


// log density of inverse gamma distribution
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_dinvgamma(
    const double& x,
    const double& alpha,
    const double& beta
) {
  
  return alpha * log(beta) - (alpha + 1) * log(x) - beta / x - lgamma(alpha);
}


// log prior density of hyper-parameters
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_prior_hyper(
    const arma::vec&  hyper,
    const arma::vec&  model,
    const Rcpp::List& prior
) {
  
  double log_prior   = 0;
  
  mat    prior_hyper = prior["hyper"];
  
  if (model(0)) {
    log_prior += log_dgamma(hyper(0), prior_hyper(0, 0), prior_hyper(0, 1));
  }
  
  if (model(1)) {
    log_prior += log_dgamma(hyper(1), prior_hyper(1, 0), prior_hyper(1, 1));
  }
  
  if (model(2)) {
    log_prior += log_dgamma(hyper(2), prior_hyper(2, 0), prior_hyper(2, 1));
  }
  
  if (model(3)) {
    for (int i = 3; i < hyper.n_elem; i++) {
      log_prior += log_dinvgamma(hyper(i), prior_hyper(i, 0), prior_hyper(i, 1));
    }
  }
  
  return log_prior;
}


// log multivariate gamma function
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_mvgamma(
    const int&    n,
    const double& x
) {
  
  if (n == 1) {
    return lgamma(x);
  }
  
  double c = (n - 1) / 2;
  return c * log(M_PI) + lgamma(x - c) + log_mvgamma(n - 1, x);
}


// log marginal likelihood, notation as in Giannone, Lenza & Primiceri (2014)
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_ml(
    const int&       p,
    const arma::mat& b,
    const arma::mat& Omega,
    const arma::mat& Psi,
    const int&       d,
    const arma::mat& Y,
    const arma::mat& X
) {
  
  int T = Y.n_rows;
  int N = Y.n_cols;
  
  double log_ml    = 0;
  mat    inv_Omega = diagmat(1 / Omega.diag());
  
  try {
    mat Bhat = inv_sympd(X.t() * X + inv_Omega) * (X.t() * Y + inv_Omega * b);
    mat ehat = Y - X * Bhat;
    
    log_ml += - N * T / 2.0 * log(M_PI);
    log_ml += log_mvgamma(N, (T + d) / 2.0) - log_mvgamma(N, d / 2.0);
    log_ml += - N / 2.0 * log_det_sympd(Omega);
    log_ml += d / 2.0 * log_det_sympd(Psi);
    log_ml += - N / 2.0 * log_det_sympd(X.t() * X + inv_Omega);
    mat A   = Psi + ehat.t() * ehat + (Bhat - b).t() * inv_Omega * (Bhat - b);
    log_ml += - (T + d) / 2.0 * log_det_sympd(A);
    
    if (!std::isfinite(log_ml)) {
      log_ml = -1e+10;
    }
    
  } catch(...) {
    log_ml = -1e+10;
  }
  
  return log_ml;
}


// log marginal likelihood with dummy observations
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_ml_dummy(
    const int&        p,
    const arma::vec&  hyper,
    const arma::vec&  model,
    const arma::mat&  Y,
    const arma::mat&  X
) {
  
  int N         = Y.n_cols;
  
  List extended = extend_dummy(p, hyper, model, Y, X);
  mat  Yplus    = as<mat>(extended["Yplus"]);
  mat  Xplus    = as<mat>(extended["Xplus"]);
  mat  Ystar    = as<mat>(extended["Ystar"]);
  mat  Xstar    = as<mat>(extended["Xstar"]);
  
  double lambda   = hyper(2);
  vec    psi      = hyper.rows(3, N + 2);
  List   prior    = mn_prior(p, lambda, psi);
  mat    prior_B  = as<mat>(prior["B"]);
  mat    prior_V  = as<mat>(prior["V"]);
  mat    prior_S  = as<mat>(prior["S"]);
  int    prior_nu = as<int>(prior["nu"]);
  
  double log_ml_plus = log_ml(p, prior_B, prior_V, prior_S, prior_nu, Yplus, Xplus);
  double log_ml_star = log_ml(p, prior_B, prior_V, prior_S, prior_nu, Ystar, Xstar);
  
  return log_ml_plus - log_ml_star;
}


// log posterior of hyper-parameters (up to a constant)
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_posterior_hyper(
    const int&        p,
    const arma::vec&  hyper,
    const arma::vec&  model,
    const arma::mat&  Y,
    const arma::mat&  X
) {
  
  // double log_prior = log_prior_hyper(hyper, model, prior);
  double log_prior = 0;
  double log_ml    = log_ml_dummy(p, hyper, model, Y, X);
  
  return log_prior + log_ml;
}


// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
arma::mat extend_hyper(
    const arma::vec& init,
    const arma::vec& model,
    const arma::mat& hypers
) {
  
  int i = 0;
  
  mat extended = repmat(init, 1, hypers.n_cols);
  
  if (model(0)) {
    extended.row(0) = hypers.row(i);
    i++;  
  }
  
  if (model(1)) {
    extended.row(1) = hypers.row(i);
    i++;
  }
  
  if (model(2)) {
    extended.row(2) = hypers.row(i);
    i++;
  }
  
  if (model(3)) {
    extended.rows(3, extended.n_rows - 1) = hypers.rows(i, hypers.n_rows - 1);
  }
  
  return extended;
}


// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
arma::mat narrow_hyper(
    const arma::vec& model,
    arma::mat        hypers
) {
  
  uvec indices;
  
  if (!model(0)) {
    indices = join_vert(indices, uvec({0}));
  }
  
  if (!model(1)) {
    indices = join_vert(indices, uvec({1}));
  }
  
  if (!model(2)) {
    indices = join_vert(indices, uvec({2}));
  }
  
  if (!model(3)) {
    indices = join_vert(indices, regspace<uvec>(3, hypers.n_rows - 1));
  }
  
  hypers.shed_rows(indices);
  
  return hypers;
}


// sample hyper-parameters
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
arma::mat sample_hyper(
    const int&        S,
    const int&        start,
    const int&        p,
    const arma::vec&  init,
    const arma::vec&  model,
    const arma::mat&  Y,
    const arma::mat&  X,
    const arma::mat&  W
) {
  
  mat hypers = metropolis(
    S, start, narrow_hyper(model, init), W,
    [p, init, model, Y, X](const vec& x) {
      vec extended = extend_hyper(init, model, x);
      return log_posterior_hyper(p, extended, model, Y, X);
    }
  );
  
  hypers = extend_hyper(init, model, hypers);
  
  return hypers;
}









