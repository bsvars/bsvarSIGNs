
#include <RcppArmadillo.h>
#include "bsvars.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube forecast_bsvarSIGNs (
    arma::cube&   posterior_Sigma,    // (N, N, S)
    arma::cube&   posterior_A,        // (N, K, S)
    arma::vec&    X_T,                // (K)
    arma::mat&    exogenous_forecast, // (horizon, d)
    arma::mat&    cond_forecast,      // (horizon, N)
    const int&    horizon
) {
  
  const int   N = posterior_Sigma.n_rows;
  const int   S = posterior_Sigma.n_slices;
  const int   K = posterior_A.n_cols;
  const int   d = exogenous_forecast.n_cols;
  
  bool        do_exog = exogenous_forecast.is_finite();
  vec         x_t;
  if ( do_exog ) {
    x_t       = X_T.rows(0, K - 1 - d);
  } else {
    x_t       = X_T.rows(0, K - 1);
  } // END if do_exog
  
  vec         Xt(K);
  cube        forecasts(N, horizon, S);
  
  for (int s=0; s<S; s++) {
    
    if ( do_exog ) {
      Xt          = join_cols(x_t, trans(exogenous_forecast.row(0)));
    } else {
      Xt          = x_t;
    } // END if do_exog
    
    for (int h=0; h<horizon; h++) {
      
      mat   Sigma             = posterior_Sigma.slice(s);
      vec   cond_forecast_h   = trans(cond_forecast.row(h));
      uvec  nonf_el           = find_nonfinite( cond_forecast_h );
      int   nonf_no           = nonf_el.n_elem;
      
      if ( nonf_no == N ) {
        forecasts.slice(s).col(h) = mvnrnd( posterior_A.slice(s) * Xt, Sigma );
      } else {
        forecasts.slice(s).col(h) = bsvars::mvnrnd_cond( cond_forecast_h, posterior_A.slice(s) * Xt, Sigma );   // does not work if cond_fc_h is all nan
      } // END if nonf_no
      
      if ( h != horizon - 1 ) {
        if ( do_exog ) {
          Xt          = join_cols( forecasts.slice(s).col(h), Xt.subvec(N, K - 1 - d), trans(cond_forecast.row(h + 1)) );
        } else {
          Xt          = join_cols( forecasts.slice(s).col(h), Xt.subvec(N, K - 1) );
        }
      } // END if h
      
    } // END h loop
  } // END s loop
  
  return forecasts;
} // END forecast_bsvarSIGNs

