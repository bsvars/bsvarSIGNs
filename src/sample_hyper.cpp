
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// update prior with hyper-parameters
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
Rcpp::List update_prior(
    const int&        p,
    const arma::vec&  hyper,
    const Rcpp::List& prior
) {
  
  int n_hyper   = hyper.n_elem;
  
  double lambda = hyper(2);
  // vec    psi    = hyper.rows(3, n_hyper - 1);
  
  mat prior_B   = as<mat>(prior["B"]);
  
  int K            = prior_B.n_rows;
  // vec v(K);
  // v(K - 1)         = 1e+6;
  // v.rows(0, K - 2) = lambda * lambda * kron(pow(linspace(1, p, p), -2), 1 / psi);
  mat prior_V      = lambda * lambda * as<mat>(prior["V"]);
  // mat prior_V      = diagmat(v);
  
  mat prior_S   = as<mat>(prior["S"]);
  // mat prior_S   = diagmat(psi);
  int prior_nu  = as<int>(prior["nu"]);
  
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
    const arma::mat& Y,
    const arma::mat& X
) {
  
  int N        = Y.n_cols;
  
  double mu    = hyper(0);
  double delta = hyper(1);
  
  mat ybar0    = mean(X.submat(0, 0, p - 1, N - 1), 0);
  mat yplus    = diagmat(ybar0 / mu);
  mat ypplus   = ybar0 / delta;
  
  mat Ystar    = join_vert(ypplus, yplus);
  mat Yplus    = join_vert(Ystar, Y);
  
  mat Xstar    = zeros(N + 1, 1);
  Xstar(0, 0)  = 1 / delta;
  Xstar        = join_horiz(repmat(Ystar, 1, p), Xstar);
  mat Xplus    = join_vert(Xstar, X);
  
  return List::create(
    _["Y"] = Yplus,
    _["X"] = Xplus
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


// log prior density of hyper-parameters (up to a constant)
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_prior_hyper(
    const arma::vec&  hyper,
    const Rcpp::List& prior
) {
  
  double log_prior   = 0;
  
  mat    prior_hyper = prior["hyper"];
  
  for (int i = 0; i < hyper.n_elem; i++) {
    // first 3 hyper-parameters ~ gamma
    if (i < 3) {
      log_prior += log_dgamma(hyper(i), prior_hyper(i, 0), prior_hyper(i, 1));
    } else {
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


// log marginal likelihood with dummy observations (up to a constant)
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_ml_dummy(
    const int&        p,
    const arma::vec&  hyper,
    const arma::mat&  Y,
    const arma::mat&  X,
    Rcpp::List        prior
) {
  
  int N         = Y.n_cols;
  int n_hyper   = hyper.n_elem;
  
  List extended = extend_dummy(p, hyper, Y, X);
  mat  Yplus    = as<mat>(extended["Y"]);
  mat  Xplus    = as<mat>(extended["X"]);
  
  prior         = update_prior(p, hyper, prior);
  mat prior_B   = as<mat>(prior["B"]);
  mat prior_V   = as<mat>(prior["V"]);
  mat prior_S   = as<mat>(prior["S"]);
  int prior_nu  = as<int>(prior["nu"]);
  
  mat inv_V     = diagmat(1 / prior_V.diag());
  
  mat Ystar     = Yplus.rows(0, N);
  mat Xstar     = Xplus.rows(0, N);
  
  double log_ml_extended = log_ml(p, prior_B, prior_V, prior_S, prior_nu, inv_V, Yplus, Xplus);
  double log_ml_dummy    = log_ml(p, prior_B, prior_V, prior_S, prior_nu, inv_V, Ystar, Xstar);
  
  return log_ml_extended - log_ml_dummy;
}


// log posterior of hyper-parameters
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_posterior_hyper(
    const int&        p,
    const arma::vec&  hyper,
    const arma::mat&  Y,
    const arma::mat&  X,
    const Rcpp::List& prior
) {
  
  double log_prior = log_prior_hyper(hyper, prior);
  double log_ml    = log_ml_dummy(p, hyper, Y, X, prior);
  
  return log_prior + log_ml;
}


// sample hyper-parameters with Metropolis-Hastings
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
arma::mat sample_hyper(
    const int&        S,
    const int&        p,
    const double&     c,
    const arma::mat&  Y,
    const arma::mat&  X,
    const Rcpp::List& prior
) {
  
  int    N       = Y.n_cols;
  int    n_hyper = N + 3;
  
  vec    hyper   = as<vec>(prior["map"]);
  
  double logp    = log_posterior_hyper(p, hyper, Y, X, prior);
  
  mat    W       = c * as<mat>(prior["W"]);
  
  mat posterior_hyper(n_hyper, S);
  
  int success = 0;
  for (int s = 0; s < S; s++) {
    
    // sample from proposal distribution
    vec    hyperp = mvnrnd(hyper, c * W);
    double logpp  = log_posterior_hyper(p, hyperp, Y, X, prior);
    double r      = exp(logpp - logp);
    
    if (randu() < r) {
      hyper = hyperp;
      logp  = logpp;
      success++;
    }
    
    posterior_hyper.col(s) = hyper;
  }
  
  cout << "Acceptance rate: " << (double) success / S << endl;
  
  return posterior_hyper;
}

