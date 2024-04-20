
#include <RcppArmadillo.h>
#include "progress.hpp"
#include "Rcpp/Rmath.h"

#include <bsvars.h>

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
    const int         thin = 100,         // introduce thinning
    const bool        show_progress = true,
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
  
  mat Yt = Y.t();
  mat Xt = X.t();
  
  const int N       = Y.n_cols;
  const int K       = X.n_cols;
  
  mat   aux_B       = as<mat>(starting_values["B"]);
  mat   aux_A       = as<mat>(starting_values["A"]);
  mat   aux_hyper   = as<mat>(starting_values["hyper"]);
  
  cube  posterior_A(N, K, S);
  cube  posterior_B(N, N, S);
  cube  posterior_hyper(2 * N + 1, 2, S);
  
  double aux_w = 1;
  vec    w     = ones(S);
  
  List post    = niw_cpp(Y, X, prior);
  mat  post_B  = as<mat>(post["B"]);
  mat  post_V  = as<mat>(post["V"]);
  mat  post_S  = as<mat>(post["S"]);
  int  post_nu = as<int>(post["nu"]);
  
  mat B, Sigma, chol_Sigma, h_invp, Q;
  
  bool success;
  int  s = 0;
  while (s < S) {

    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    Sigma      = iwishrnd(post_S, post_nu);
    chol_Sigma = chol(Sigma, "lower");
    B          = rmatnorm_cpp(post_B, post_V, Sigma);
    h_invp     = inv(trimatl(chol_Sigma));    // lower triangular, h(Sigma) is upper triangular
    
    success = false;
    Q       = sample_Q(lags, Y, X, 
                       aux_w, B, h_invp, chol_Sigma, 
                       prior, VB,
                       sign_irf, sign_narrative, sign_B, Z,
                       max_tries, success);
    
    if (success) {
      // Increment progress bar
      if (any(prog_rep_points == s)) p.increment();
      
      posterior_A.slice(s) = B.t();
      posterior_B.slice(s) = Q.t() * h_invp;
      w(s)                 = aux_w;
      s++;
    }
  } // END s loop
  
  return List::create(
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      _["A"]        = posterior_A,
      _["hyper"]    = posterior_hyper,
      _["w"]        = w
    )
  );
} // END bsvar_sign_cpp



