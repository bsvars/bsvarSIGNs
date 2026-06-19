#include <RcppArmadillo.h>
#include "Rcpp/Rmath.h"
#include <bsvars.h>

#include "sample_hyper.h"
#include "sample_Q.h"
#include "sample_NIW.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List bsvar_sign_par_cpp(
    const int&        p,                  // number of lags
    const arma::mat&  Y,                  // TxN dependent variables
    const arma::mat&  X,                  // TxK dependent variables
    const arma::cube& sign_irf,           // NxNxh cube of signs for impulse response function
    const arma::mat&  sign_narrative,     // ANYx6 matrix of signs for historical decomposition
    const arma::mat&  sign_B,             // NxN matrix of signs for B
    const arma::field<arma::mat>& Z,      // a list of zero restrictions
    const int&        Nf,                 // number of foreign variables for SOE
    const Rcpp::List& prior,              // a list of priors
    const int&        max_tries = 10000   // maximum tries for Q draw
) {
  
  const int  T = Y.n_rows;
  const int  N = Y.n_cols;
  const int  K = X.n_cols;
  
  mat        hypers = as<mat>(prior["hyper"]);
  int        S_hyper  = hypers.n_cols - 1;
  int        prior_nu = as<int>(prior["nu"]);
  int        post_nu  = prior_nu + T;
  int        n_tries;
  
  double     w, mu, delta, lambda;
  
  vec        hyper, psi;
  vec        prior_v = as<mat>(prior["V"]).diag();
  
  mat        B, Sigma, chol_Sigma, h_invp, Q, shocks;
  mat        prior_V, prior_S, post_B, post_V, post_S;
  mat        Ystar, Xstar, Yplus, Xplus;
  mat        prior_B = as<mat>(prior["B"]);
  mat        Ysoc    = as<mat>(prior["Ysoc"]);
  mat        Xsoc    = as<mat>(prior["Xsoc"]);
  mat        Ysur    = as<mat>(prior["Ysur"]);
  mat        Xsur    = as<mat>(prior["Xsur"]);
  
  field<mat> result;
  
  hyper        = hypers.col(randi(distr_param(0, S_hyper)));
  mu           = hyper(0);
  delta        = hyper(1);
  lambda       = hyper(2);
  psi          = hyper.rows(3, N + 2);
  
  // update Minnesota prior
  prior_V      = diagmat(prior_v % join_vert(lambda * lambda * repmat(1 / psi, p, 1),
                                             ones<vec>(K - N * p)));
  prior_S      = diagmat(psi);
  
  // update dummy observation prior
  Ystar        = join_vert(Ysoc / mu, Ysur / delta);
  Xstar        = join_vert(Xsoc / mu, Xsur / delta);
  
  mat Y_scaled = Y;
  mat X_scaled = X;
  int covid = as<int>(prior["covid"]);
  if (covid > 0 && covid <= T) {
    int c_idx = covid - 1;
    double s0 = hyper(N + 3);
    double s1 = hyper(N + 4);
    double s2 = hyper(N + 5);
    double rho = hyper(N + 6);
    
    vec s = ones<vec>(T);
    if (c_idx < T) s(c_idx) = s0;
    if (c_idx + 1 < T) s(c_idx + 1) = s1;
    if (c_idx + 2 < T) s(c_idx + 2) = s2;
    for (int t = c_idx + 3; t < T; t++) {
      s(t) = 1.0 + (s2 - 1.0) * std::pow(rho, t - c_idx - 2);
    }
    
    Y_scaled.each_col() /= s;
    X_scaled.each_col() /= s;
  }
  
  Yplus        = join_vert(Ystar, Y_scaled);
  Xplus        = join_vert(Xstar, X_scaled);
  
  // posterior parameters
  result       = niw_cpp(Yplus, Xplus, prior_B, prior_V, prior_S, prior_nu);
  post_B       = result(0);
  post_V       = result(1);
  post_S       = result(2);
  post_nu      = as_scalar(result(3));
  
  w            = 0;
  n_tries      = 0;
  
  while (w == 0 and (n_tries < max_tries or max_tries == 0)) {
    
    checkUserInterrupt();
    
    // sample reduced-form parameters
    Sigma      = iwishrnd(post_S, post_nu);
    chol_Sigma = chol(Sigma, "lower");
    B          = rmatnorm_cpp(post_B, post_V, Sigma);
    h_invp     = inv(trimatl(chol_Sigma)); // lower tri, h(Sigma) is upper tri
    
    result     = sample_Q(p, Y_scaled, X_scaled, B, h_invp, chol_Sigma, prior, 
                          sign_irf, sign_narrative, sign_B, Z, Nf, 1);
    Q          = result(0);
    shocks     = result(1);
    w          = as_scalar(result(2));
    n_tries++;
  }
  
  return List::create(
    _["w"]        = w,
    _["hyper"]    = hyper,
    _["A"]        = B.t(),
    _["B"]        = Q.t() * h_invp,
    _["Q"]        = Q,
    _["Sigma"]    = Sigma,
    _["Theta0"]   = chol_Sigma * Q,
    _["shocks"]   = shocks
  );
} // END bsvar_sign_par_cpp
