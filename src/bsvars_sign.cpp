
#include <RcppArmadillo.h>
#include "progress.hpp"
#include "Rcpp/Rmath.h"

#include <bsvars.h>

#include "sample_Q.h"

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
  
  const int N       = Y.n_rows;
  const int K       = X.n_rows;
  
  mat   aux_B       = as<mat>(starting_values["B"]);
  mat   aux_A       = as<mat>(starting_values["A"]);
  mat   aux_hyper   = as<mat>(starting_values["hyper"]);
  mat   Q(N, N);
  
  const int   SS    = floor(S / thin);
  
  cube  posterior_B(N, N, SS);
  cube  posterior_A(N, K, SS);
  cube  posterior_hyper(2 * N + 1, 2, SS);
  cube  posterior_Q(N, N, SS);
  
  int    n_fails    = 0;
  int    ss         = 0;
  double aux_w      = 1;
  vec    w          = ones(SS);
  
  for (int s=0; s<S; s++) {
    
    // Increment progress bar
    if (any(prog_rep_points == s)) p.increment();
    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    aux_hyper     = bsvars::sample_hyperparameters(aux_hyper, aux_B, aux_A, VB, prior);
    aux_A         = bsvars::sample_A_homosk1(aux_A, aux_B, aux_hyper, Y, X, prior);
    aux_B         = bsvars::sample_B_homosk1(aux_B, aux_A, aux_hyper, Y, X, prior, VB);
    
    if (s % thin == 0) {
      Q                          = sample_Q(lags, Y, X, 
                                            aux_w, aux_A, aux_B, aux_hyper, 
                                            prior, VB,
                                            sign_irf, sign_narrative, sign_B,
                                            max_tries, n_fails);
      
      posterior_B.slice(ss)      = aux_B;
      posterior_A.slice(ss)      = aux_A;
      posterior_hyper.slice(ss)  = aux_hyper;
      posterior_Q.slice(ss)      = Q;
      w(ss)                      = aux_w;
      ss++;
    }
  } // END s loop
  
  double fail = static_cast<double>(n_fails) / SS;
  
  return List::create(
    _["last_draw"]  = List::create(
      _["B"]        = aux_B,
      _["A"]        = aux_A,
      _["hyper"]    = aux_hyper,
      _["Q"]        = Q
    ),
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      _["A"]        = posterior_A,
      _["hyper"]    = posterior_hyper,
      _["Q"]        = posterior_Q,
      _["w"]        = w,
      _["fail"]     = fail
    )
  );
} // END bsvar_sign_cpp
