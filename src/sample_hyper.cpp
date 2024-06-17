
#include <RcppArmadillo.h>

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


// update prior with hyper-parameters
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
Rcpp::List update_prior(
    const int&        p,
    const arma::vec&  hyper,
    const arma::vec&  model,
    const Rcpp::List& prior
) {
  
  mat prior_B   = as<mat>(prior["B"]);
  int prior_nu  = as<int>(prior["nu"]);
  
  mat prior_V, prior_S;
  if (model(3)) {
    int K       = prior_B.n_rows;
    int n_hyper = hyper.n_elem;
    
    vec psi     = hyper.rows(3, n_hyper - 1);
    prior_S     = diagmat(psi);
    
    vec v(K);
    v(K - 1)         = 1e+6;
    v.rows(0, K - 2) = kron(pow(linspace(1, p, p), -2), 1 / psi);
    prior_V          = diagmat(v);
  } else {
    prior_V = as<mat>(prior["V"]);
    prior_S = as<mat>(prior["S"]);
  }
  
  double lambda;
  if (model(2)) {
    lambda = hyper(2);
  } else {
    lambda = 0.2;
  }
  prior_V *= lambda * lambda;
  
  return List::create(
    _["B"]  = prior_B,
    _["V"]  = prior_V,
    _["S"]  = prior_S,
    _["nu"] = prior_nu
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


// log marginal likelihood
// notation as in Giannone, Lenza & Primiceri (2014)
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_ml(
    const int&       p,
    const arma::mat& b,
    const arma::mat& Omega,
    const arma::mat& Psi,
    const int&       d,
    const arma::mat& inv_Omega,
    const arma::mat& Y,
    const arma::mat& X
) {
  
  int T = Y.n_rows;
  int N = Y.n_cols;
  
  mat Bhat = solve(X.t() * X + inv_Omega, X.t() * Y + inv_Omega * b);
  mat ehat = Y - X * Bhat;
  
  double log_ml = 0;
  
  log_ml += - N * T / 2 * log(M_PI);
  log_ml += log_mvgamma(N, (T + d) / 2) - log_mvgamma(N, d / 2);
  log_ml += - N / 2 * log_det_sympd(Omega);
  log_ml += d / 2 * log_det_sympd(Psi);
  log_ml += - N / 2 * log_det_sympd(X.t() * X + inv_Omega);
  mat A   = ehat.t() * ehat + (Bhat - b).t() * inv_Omega * (Bhat - b);
  log_ml += - (T + d) / 2 * log_det_sympd(Psi + A);
  
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
    const arma::mat&  X,
    Rcpp::List        prior
) {
  
  int N         = Y.n_cols;
  int n_hyper   = hyper.n_elem;
  
  List extended = extend_dummy(p, model, hyper, Y, X);
  mat  Yplus    = as<mat>(extended["Yplus"]);
  mat  Xplus    = as<mat>(extended["Xplus"]);
  mat  Ystar    = as<mat>(extended["Ystar"]);
  mat  Xstar    = as<mat>(extended["Xstar"]);
  
  prior         = update_prior(p, hyper, model, prior);
  mat prior_B   = as<mat>(prior["B"]);
  mat prior_V   = as<mat>(prior["V"]);
  mat prior_S   = as<mat>(prior["S"]);
  int prior_nu  = as<int>(prior["nu"]);
  
  mat inv_V     = diagmat(1 / prior_V.diag());
  
  double log_ml_extended = log_ml(p, prior_B, prior_V, prior_S, prior_nu, inv_V, Yplus, Xplus);
  double log_ml_dummy    = log_ml(p, prior_B, prior_V, prior_S, prior_nu, inv_V, Ystar, Xstar);
  
  return log_ml_extended - log_ml_dummy;
}


// log posterior of hyper-parameters (up to a constant)
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_posterior_hyper(
    const int&        p,
    const arma::vec&  hyper,
    const arma::vec&  model,
    const arma::mat&  Y,
    const arma::mat&  X,
    const Rcpp::List& prior
) {
  
  double log_prior = log_prior_hyper(hyper, model, prior);
  double log_ml    = log_ml_dummy(p, hyper, model, Y, X, prior);
  
  return log_prior + log_ml;
}


