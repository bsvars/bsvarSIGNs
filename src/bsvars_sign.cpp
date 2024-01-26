
#include <RcppArmadillo.h>
#include "progress.hpp"
#include "Rcpp/Rmath.h"

#include <bsvars.h>
#include "utils.h"

using namespace Rcpp;
using namespace arma;


// If matches traditional sign restriction
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
bool match_sign_cpp(const arma::mat& A, const arma::mat sign) {
 return accu(((A % sign) > 0)) == accu(abs(sign));
}


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List bsvar_sign_cpp(
    const int&  S,                        // number of draws from the posterior
    const int&  lags,                     // number of lags
    const arma::mat&  Y,                  // NxT dependent variables
    const arma::mat&  X,                  // KxT dependent variables
    const arma::field<arma::mat>& VB,     // N-list
    const arma::cube& sign_irf,           // NxNxh cube of signs for impulse response function
    const arma::mat&  sign_hd,            // Mx6 matrix of signs for historical decomposition
    const Rcpp::List& prior,              // a list of priors
    const Rcpp::List& starting_values,    // a list of starting values
    const int         thin = 100,         // introduce thinning
    const bool        show_progress = true
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
  const int h       = sign_irf.n_slices;
  
  mat   aux_B       = as<mat>(starting_values["B"]);
  mat   aux_A       = as<mat>(starting_values["A"]);
  mat   aux_hyper   = as<mat>(starting_values["hyper"]);
  mat   Q(N, N);
  
  const int   SS    = floor(S / thin);
  
  cube  posterior_B(N, N, SS);
  cube  posterior_A(N, K, SS);
  cube  posterior_hyper(2 * N + 1, 2, SS);
  cube  irf(N, N, h);
  
  int  ss        = 0;
  
  int  max_tries = pow(10, 6);
  int  n_tries;
  bool success;
  
  for (int s=0; s<S; s++) {
    
    // Increment progress bar
    if (any(prog_rep_points == s)) p.increment();
    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    aux_hyper     = bsvars::sample_hyperparameters(aux_hyper, aux_B, aux_A, VB, prior);
    aux_A         = bsvars::sample_A_homosk1(aux_A, aux_B, aux_hyper, Y, X, prior);
    aux_B         = bsvars::sample_B_homosk1(aux_B, aux_A, aux_hyper, Y, X, prior, VB);
    
    if (s % thin == 0) {
      ////////////////////////////////////////////////////////////////////////////////////
      // sign restrictions
      
      irf     = bsvars::bsvars_ir1(aux_B, aux_A, h-1, lags);
      n_tries = 0;
      success = true;
      
      while (n_tries < max_tries) {
        Q = rortho_cpp(N);
        
        // traditional sign restrictions on impulse response functions
        for (int hh=0; hh<h; hh++) {
          if ( !match_sign_cpp(irf.slice(hh) * Q, sign_irf.slice(hh)) ) {
            success = false;
            break;
          }
          
        }
        // END traditional sign restrictions
        

        // narrative sign restrictions on historical decomposition
        // TODO
        
        // END narrative sign restrictions
        
        if (success) break;
        n_tries++;
      }
      
      if (!success) {
        Rcout << "Warning: could not find a valid Q matrix after " << max_tries << " tries." << endl;
        continue;
      }
      ////////////////////////////////////////////////////////////////////////////////////
      
      posterior_B.slice(ss)      = Q.t() * aux_B;
      posterior_A.slice(ss)      = aux_A;
      posterior_hyper.slice(ss)  = aux_hyper;
      ss++;
    }
  } // END s loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["B"]        = aux_B,
      _["A"]        = aux_A,
      _["hyper"]    = aux_hyper
    ),
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      _["A"]        = posterior_A,
      _["hyper"]    = posterior_hyper
    )
  );
} // END bsvar_sign_cpp
