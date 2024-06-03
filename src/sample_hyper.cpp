
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
  vec    phi    = hyper.rows(3, n_hyper - 1);
  
  mat prior_B   = as<mat>(prior["B"]);
  mat prior_V   = lambda * lambda * as<mat>(prior["V"]);
  mat prior_S   = diagmat(phi);
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


// log density of gamma distribution (up to a constant)
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_dgamma(
    const double& x,
    const double& alpha,  // shape
    const double& beta    // rate
) {
  
  return (alpha - 1) * log(x) - beta * x;
}


// log density of inverse gamma distribution (up to a constant)
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double log_dinvgamma(
    const double& x,
    const double& alpha,  // shape
    const double& beta    // scale
) {
  
  return -(alpha + 1) * log(x) - beta / x;
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


// log marginal likelihood (up to a constant)
// notation as in Domenico, Lenza & Primiceri (2014)
double log_ml(
    const int& p,
    const mat& b,
    const mat& Omega,
    const mat& Phi,
    const int& d,
    const mat& D_Omega,
    const mat& D_Phi,
    const mat& inv_Omega,
    const mat& Y,
    const mat& X
) {
  
  int T = Y.n_rows + p;
  int N = Y.n_cols;
  int K = X.n_cols;
  
  mat Bhat = solve(X.t() * X + inv_Omega, X.t() * Y + inv_Omega * b);
  mat ehat = Y - X * Bhat;
  
  double log_ml = 0;
  
  // cancels out in p(Yplus) / p(Ystar)
  // log_ml -= (T - p) / 2 * log_det_sympd(Phi);
  
  log_ml -= N / 2 * log_det(D_Omega.t() * X.t() * X * D_Omega + eye(K, K)).real();
  // log_ml -= N / 2 * log_det(X.t() * X + inv_Omega).real();
  
  mat A   = ehat.t() * ehat + (Bhat - b).t() * inv_Omega * (Bhat - b);
  log_ml -= (T - p + d) / 2 * log_det(D_Phi.t() * A * D_Phi + eye(N, N)).real();
  // log_ml -= (T - p + d) / 2 * log_det(Phi + A).real();
  
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
  
  vec    phi    = hyper.rows(3, n_hyper - 1);  
  
  List extended = extend_dummy(p, hyper, Y, X);
  mat  Yplus    = as<mat>(extended["Y"]);
  mat  Xplus    = as<mat>(extended["X"]);
  
  prior         = update_prior(p, hyper, prior);
  mat prior_B   = as<mat>(prior["B"]);
  mat prior_V   = as<mat>(prior["V"]);
  mat prior_S   = as<mat>(prior["S"]);
  int prior_nu  = as<int>(prior["nu"]);
  
  // mat chol_V    = chol(prior_V, "lower");
  // mat inv_V     = inv_sympd(prior_V);
  mat chol_V    = diagmat(sqrt(prior_V.diag()));
  mat inv_V     = diagmat(1 / prior_V.diag());
  mat cholinv_S = diagmat(sqrt(1 / phi));
  
  mat Ystar     = Yplus.rows(0, N);
  mat Xstar     = Xplus.rows(0, N);
  
  double log_ml_extended = log_ml(p, prior_B, prior_V, prior_S, prior_nu, 
                                  chol_V, cholinv_S, inv_V, Yplus, Xplus);
  double log_ml_dummy    = log_ml(p, prior_B, prior_V, prior_S, prior_nu, 
                                  chol_V, cholinv_S, inv_V, Ystar, Xstar);
  
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
    const arma::mat&  Y,
    const arma::mat&  X,
    const Rcpp::List& prior
) {
  
  int    N       = Y.n_cols;
  int    n_hyper = N + 3;
  
  vec    hyper   = as<vec>(prior["map"]);
  
  double logp    = log_posterior_hyper(p, hyper, Y, X, prior);
  
  double c       = 1;
  mat    W       = as<mat>(prior["W"]);
  
  mat posterior_hyper(n_hyper, S);
  
  for (int s = 0; s < S; s++) {
    
    // sample from proposal distribution
    vec    hyperp = mvnrnd(hyper, c * W);
    double logpp  = log_posterior_hyper(p, hyperp, Y, X, prior);
    double r      = exp(logpp - logp);
    
    if (r < 0.2) {
      c /= 10;
    } else if (r > 0.5) {
      c *= 2;
    }
    
    if (randu() < r) {
      hyper = hyperp;
      logp  = logpp;
    }
    
    posterior_hyper.col(s) = hyper;
  }
  
  return posterior_hyper;
}

