
#include <RcppArmadillo.h>

#include "utils.h"
#include "mcmc.h"

using namespace Rcpp;
using namespace arma;


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
  
  double log_prior = 0, shape, scale;
  
  if (model(0)) {
    shape      = as<double>(prior["mu.shape"]);
    scale      = as<double>(prior["mu.scale"]);
    log_prior += log_dgamma(hyper(0), shape, scale);
  }
  
  if (model(1)) {
    shape      = as<double>(prior["delta.shape"]);
    scale      = as<double>(prior["delta.scale"]);
    log_prior += log_dgamma(hyper(1), shape, scale);
  }
  
  if (model(2)) {
    shape      = as<double>(prior["lambda.shape"]);
    scale      = as<double>(prior["lambda.scale"]);
    log_prior += log_dgamma(hyper(2), shape, scale);
  }
  
  if (model(3)) {
    shape      = as<double>(prior["psi.shape"]);
    scale      = as<double>(prior["psi.scale"]);
    for (int i = 3; i < hyper.n_elem; i++) {
      log_prior += log_dinvgamma(hyper(i), shape, scale);
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
    log_ml += log_mvgamma(N, (T + d) / 2.0);
    log_ml += -log_mvgamma(N, d / 2.0);
    log_ml += - N / 2.0 * log_det_sympd(Omega);
    log_ml += d / 2.0 * log_det_sympd(Psi);
    log_ml += - N / 2.0 * log_det_sympd(X.t() * X + inv_Omega);
    mat A   = Psi + ehat.t() * ehat + (Bhat - b).t() * inv_Omega * (Bhat - b);
    log_ml += - (T + d) / 2.0 * log_det_sympd(A);
    
  } catch(...) {
    log_ml = -1e+10;
  }
  
  return log_ml;
}


// log marginal likelihood with dummy observations
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_ml_dummy(
    const arma::vec&  hyper,
    const arma::vec&  model,
    const arma::mat&  Y,
    const arma::mat&  X,
    const Rcpp::List& prior
) {
  
  int    N           = Y.n_cols;
  double mu          = hyper(0);
  double delta       = hyper(1);
  double lambda      = hyper(2);
  vec    psi         = hyper.rows(3, N + 2);
  
  // update Minnesota prior
  mat    prior_B     = as<mat>(prior["B"]);
  mat    prior_V     = diagmat(join_vert(
                               lambda*lambda * kron(as<vec>(prior["Vp"]), 1 / psi),
                               as<vec>(prior["Vd"])));
  mat    prior_S     = diagmat(psi);
  int    prior_nu    = as<int>(prior["nu"]);
  
  // update dummy observation prior
  mat    Ystar       = join_vert(as<mat>(prior["Ysoc"]) / mu, 
                                 as<mat>(prior["Ysur"]) / delta);
  mat    Xstar       = join_vert(as<mat>(prior["Xsoc"]) / mu, 
                                 as<mat>(prior["Xsur"]) / delta);
  mat    Yplus       = join_vert(Ystar, Y);
  mat    Xplus       = join_vert(Xstar, X);
  
  double log_ml_plus = log_ml(prior_B, prior_V, prior_S, prior_nu, Yplus, Xplus);
  double log_ml_star = log_ml(prior_B, prior_V, prior_S, prior_nu, Ystar, Xstar);
  
  return log_ml_plus - log_ml_star;
}


// log posterior of hyper-parameters (up to a constant)
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_posterior_hyper(
    const arma::vec&  hyper,
    const arma::vec&  model,
    const arma::mat&  Y,
    const arma::mat&  X,
    const Rcpp::List& prior
) {
  
  double log_prior = log_prior_hyper(hyper, model, prior);
  double log_ml    = log_ml_dummy(hyper, model, Y, X, prior);
  double log_post  = log_prior + log_ml;
  
  if (!std::isfinite(log_post)) {
    log_post       = -1e+10;
  }
  
  return log_post;
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
    const arma::vec&  init,
    const arma::vec&  model,
    const arma::mat&  Y,
    const arma::mat&  X,
    const arma::mat&  W,
    const Rcpp::List& prior
) {
  
  mat hypers = metropolis(
    
    S, start, narrow_hyper(model, init), W,
    
    [init, model, Y, X, prior](const vec& x) {
      
      vec extended = extend_hyper(init, model, x);
      return log_posterior_hyper(extended, model, Y, X, prior);
    }
  );
  
  hypers = extend_hyper(init, model, hypers);
  
  return hypers;
}

