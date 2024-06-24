
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
  
  while (s < S) {

    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    hyper      = hypers.col(randi(distr_param(0, S_hyper)));
    lambda     = hyper(2);
    psi        = hyper.rows(3, N + 2);
    
    result     = niw_cpp(Y, X, mn_prior(lags, lambda, psi));
    post_B     = as<mat>(result["B"]);
    post_V     = as<mat>(result["V"]);
    post_S     = as<mat>(result["S"]);
    post_nu    = as<int>(result["nu"]);
    
    Sigma      = iwishrnd(post_S, post_nu);
    chol_Sigma = chol(Sigma, "lower");
    B          = rmatnorm_cpp(post_B, post_V, Sigma);
    h_invp     = inv(trimatl(chol_Sigma));    // lower triangular, h(Sigma) is upper triangular
    
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



