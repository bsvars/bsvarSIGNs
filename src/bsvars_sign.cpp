
#include <RcppArmadillo.h>
#include "progress.hpp"
#include "Rcpp/Rmath.h"

#include <bsvars.h>

#include "sample_hyper.h"
#include "sample_Q.h"
#include "sample_NIW.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List bsvar_sign_cpp(
    const int&        S,                  // number of draws from the posterior
    const int&        lags,               // number of lags
    const arma::mat&  Y,                  // NxT dependent variables
    const arma::mat&  X,                  // KxT dependent variables
    const arma::field<arma::mat>& VB,     // N-list
    const arma::cube& sign_irf,           // NxNxh cube of signs for impulse response function
    const arma::mat&  sign_narrative,     // Mx6 matrix of signs for historical decomposition
    const arma::mat&  sign_B,             // Mx6 matrix of signs for B
    const arma::field<arma::mat>& Z,      // a list of zero restrictions
    const Rcpp::List& prior,              // a list of priors
    const Rcpp::List& starting_values,    // a list of starting values
    const bool        show_progress = true,
    const int         thin = 100,         // introduce thinning
    const int&        max_tries = 10000   // maximum tries for Q draw
) {
  
  std::string oo = "";
  if ( thin != 1 ) {
    oo      = bsvars::ordinal(thin) + " ";
  }
  
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, S, 50));
  if (show_progress) {
    Rcout << "**************************************************|" << endl;
    Rcout << "bsvars: Bayesian Structural Vector Autoregressions|" << endl;
    Rcout << "**************************************************|" << endl;
    Rcout << " Gibbs sampler for the SVAR model                 |" << endl;
    Rcout << "**************************************************|" << endl;
    Rcout << " Progress of the MCMC simulation for " << S << " draws" << endl;
    Rcout << "    Every " << oo << "draw is saved via MCMC thinning" << endl;
    Rcout << " Press Esc to interrupt the computations" << endl;
    Rcout << "**************************************************|" << endl;
  }
  Progress p(50, show_progress);
  
  const int T = Y.n_rows;
  const int N = Y.n_cols;
  const int K = X.n_cols;
  
  vec   posterior_w(S);
  mat   posterior_hyper(N + 3, S);
  cube  posterior_A(N, K, S);
  cube  posterior_B(N, N, S);
  cube  posterior_Sigma(N, N, S);
  cube  posterior_Theta0(N, N, S);
  cube  posterior_shocks(N, T, S);
  
  mat   hypers = as<mat>(prior["hyper"]);
  
  int    s = 0, S_hyper = hypers.n_cols - 1, post_nu;
  double w = 1, lambda;
  vec    hyper, psi;
  mat    B, Sigma, chol_Sigma, h_invp, Q, Epsilon, post_B, post_V, post_S;
  List   result;
  
  mat        prior_B  = as<mat>(prior["B"]), prior_V, prior_S;
  int        prior_nu = as<int>(prior["nu"]);
  field<mat> post;
  
  while (s < S) {

    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    // update Minnesota prior
    hyper      = hypers.col(randi(distr_param(0, S_hyper)));
    lambda     = hyper(2);
    psi        = hyper.rows(3, N + 2);
    
    prior_V    = diagmat(join_vert(
      lambda*lambda * kron(as<vec>(prior["Vp"]), 1 / psi),
      as<vec>(prior["Vd"])  
    ));
    prior_S    = diagmat(psi);
    
    // update dummy observation prior
    mat    Ystar(0, N), Xstar(0, K);
    double mu       = hyper(0);
    Ystar           = join_vert(Ystar, as<mat>(prior["Ysoc"]) / mu);
    Xstar           = join_vert(Xstar, as<mat>(prior["Xsoc"]) / mu);
    
    double delta    = hyper(1);
    Ystar           = join_vert(Ystar, as<mat>(prior["Ysur"]) / delta);
    Xstar           = join_vert(Xstar, as<mat>(prior["Xsur"]) / delta);
    
    mat    Yplus    = join_vert(Ystar, Y);
    mat    Xplus    = join_vert(Xstar, X);
  
    // posterior parameters
    post       = niw_cpp(Yplus, Xplus, prior_B, prior_V, prior_S, prior_nu);
    post_B     = post(0);
    post_V     = post(1);
    post_S     = post(2);
    post_nu    = as_scalar(post(3));
    
    // sample reduced-form parameters (B, Sigma, Q)
    Sigma      = iwishrnd(post_S, post_nu);
    chol_Sigma = chol(Sigma, "lower");
    B          = rmatnorm_cpp(post_B, post_V, Sigma);
    h_invp     = inv(trimatl(chol_Sigma)); // lower tri, h(Sigma) is upper tri
    
    result     = sample_Q(lags, Y, X, 
                          B, h_invp, chol_Sigma,
                          prior, VB,
                          sign_irf, sign_narrative, sign_B, Z,
                          max_tries);
    Q          = as<mat>(result["Q"]);
    w          = as<double>(result["w"]);
    Epsilon    = as<mat>(result["Epsilon"]);
    
    if (w > 0) {
      // Increment progress bar
      if (any(prog_rep_points == s)) p.increment();
      
      posterior_w(s)            = w;
      posterior_A.slice(s)      = B.t();
      posterior_B.slice(s)      = Q.t() * h_invp;
      posterior_Sigma.slice(s)  = Sigma;
      posterior_Theta0.slice(s) = chol_Sigma * Q;
      posterior_shocks.slice(s) = Epsilon;
      posterior_hyper.col(s)    = hyper;
      
      s++;
    }
  } // END s loop
  
  return List::create(
    _["posterior"]  = List::create(
      _["w"]        = posterior_w,
      _["hyper"]    = posterior_hyper,
      _["A"]        = posterior_A,
      _["B"]        = posterior_B,
      _["Sigma"]    = posterior_Sigma,
      _["Theta0"]   = posterior_Theta0,
      _["shocks"]   = posterior_shocks
    )
  );
} // END bsvar_sign_cpp



