
#include <RcppArmadillo.h>
#include "progress.hpp"
#include "Rcpp/Rmath.h"
#include <bsvars.h>

#ifdef _OPENMP
  #include <omp.h>
#else
  // for machines with compilers void of openmp support
  #define omp_get_num_threads()  1
  #define omp_get_thread_num()   0
  #define omp_get_max_threads()  1
  #define omp_get_thread_limit() 1
  #define omp_get_num_procs()    1
  #define omp_set_nested(a)   // empty statement to remove the call
  #define omp_get_wtime()        0
#endif

#include "sample_hyper.h"
#include "sample_Q.h"
#include "sample_NIW.h"

using namespace Rcpp;
using namespace arma;


/***
All notations in the C++ code except for compute.cpp and forecast_bsvarSIGNs.cpp
 are consistent with the notations in the papers:
 Antolín-Díaz and Rubio-Ramírez (2018) and Arias, Rubio-Ramírez and Waggoner (2018)
 which are different from the notations in the R code.
***/


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List bsvar_sign_cpp(
    const int&        S,                  // number of draws from the posterior
    const int&        p,                  // number of lags
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
  double num_threads;
  #pragma omp parallel
  {
    num_threads = omp_get_num_threads();
  }
  vec prog_rep_points = arma::round(arma::linspace(0, S / num_threads, 50));
  if (show_progress) {
    Rcout << "**************************************************|" << endl;
    Rcout << " bsvarSIGNs: Bayesian Structural VAR with sign,   |" << endl;
    Rcout << "             zero and narrative restrictions      |" << endl;
    Rcout << "**************************************************|" << endl;
    // Rcout << " Gibbs sampler for the SVAR model                 |" << endl;
    // Rcout << "**************************************************|" << endl;
    Rcout << " Progress of simulation for " << S << " independent draws" << endl;
    // Rcout << "    Every " << oo << "draw is saved via MCMC thinning" << endl;
    Rcout << " Press Esc to interrupt the computations" << endl;
    Rcout << "**************************************************|" << endl;
  }
  Progress bar(50, show_progress);
  
  const int  T = Y.n_rows;
  const int  N = Y.n_cols;
  const int  K = X.n_cols;
  
  vec        posterior_w(S);
  mat        posterior_hyper(N + 3, S);
  cube       posterior_A(N, K, S);
  cube       posterior_B(N, N, S);
  cube       posterior_Q(N, N, S);
  cube       posterior_Sigma(N, N, S);
  cube       posterior_Theta0(N, N, S);
  cube       posterior_shocks(N, T, S);
  
  mat        hypers = as<mat>(prior["hyper"]);
  
  int        S_hyper  = hypers.n_cols - 1;
  int        prior_nu = as<int>(prior["nu"]);
  
  vec        prior_v = as<mat>(prior["V"]).diag();
  
  mat        prior_B = as<mat>(prior["B"]);
  mat        Ysoc    = as<mat>(prior["Ysoc"]);
  mat        Xsoc    = as<mat>(prior["Xsoc"]);
  mat        Ysur    = as<mat>(prior["Ysur"]);
  mat        Xsur    = as<mat>(prior["Xsur"]);
  
  #pragma omp parallel for
  for (int s = 0; s < S; s++) {
    
    vec    hyper      = hypers.col(randi(distr_param(0, S_hyper)));
    double mu         = hyper(0);
    double delta      = hyper(1);
    double lambda     = hyper(2);
    vec    psi        = hyper.rows(3, N + 2);
    
    // update Minnesota prior
    mat    prior_V    = diagmat(prior_v % join_vert(lambda * lambda * repmat(1 / psi, p, 1),
                                                    ones<vec>(K - N * p)));
    mat    prior_S    = diagmat(psi);
    
    // update dummy observation prior
    mat    Ystar      = join_vert(Ysoc / mu, Ysur / delta);
    mat    Xstar      = join_vert(Xsoc / mu, Xsur / delta);
    mat    Yplus      = join_vert(Ystar, Y);
    mat    Xplus      = join_vert(Xstar, X);
    
    // posterior parameters
    field<mat> result = niw_cpp(Yplus, Xplus, prior_B, prior_V, prior_S, prior_nu);
    mat    post_B     = result(0);
    mat    post_V     = result(1);
    mat    post_S     = result(2);
    int    post_nu    = as_scalar(result(3));
    
    double w          = 0;
    
    mat    Sigma, chol_Sigma, B, h_invp, Q, shocks;
    
    while (w == 0) {
      
      checkUserInterrupt();
      
      // sample reduced-form parameters
      Sigma      = iwishrnd(post_S, post_nu);
      chol_Sigma = chol(Sigma, "lower");
      B          = rmatnorm_cpp(post_B, post_V, Sigma);
      h_invp     = inv(trimatl(chol_Sigma)); // lower tri, h(Sigma) is upper tri
      
      result     = sample_Q(p, Y, X, B, h_invp, chol_Sigma, prior, VB, 
                            sign_irf, sign_narrative, sign_B, Z, max_tries);
      Q          = result(0);
      shocks     = result(1);
      w          = as_scalar(result(2));
    }
    
    posterior_w(s)            = w;
    posterior_hyper.col(s)    = hyper;
    posterior_A.slice(s)      = B.t();
    posterior_B.slice(s)      = Q.t() * h_invp;
    posterior_Q.slice(s)      = Q;
    posterior_Sigma.slice(s)  = Sigma;
    posterior_Theta0.slice(s) = chol_Sigma * Q;
    posterior_shocks.slice(s) = shocks;
    
    if (omp_get_thread_num() == 0 and any(prog_rep_points == s)) {
      bar.increment();
    }
    
  } // END s loop
  
  return List::create(
    _["posterior"]  = List::create(
      _["w"]        = posterior_w,
      _["hyper"]    = posterior_hyper,
      _["A"]        = posterior_A,
      _["B"]        = posterior_B,
      _["Q"]        = posterior_Q,
      _["Sigma"]    = posterior_Sigma,
      _["Theta0"]   = posterior_Theta0,
      _["shocks"]   = posterior_shocks
    )
  );
} // END bsvar_sign_cpp

