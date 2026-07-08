
#include <RcppArmadillo.h>
#include "bsvars.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List forecast_bsvarSIGNs (
    arma::cube&   posterior_Sigma,    // (N, N, S)
    arma::cube&   posterior_A,        // (N, K, S)
    arma::mat&    posterior_hyper,    // (N+7, S)
    arma::vec&    X_T,                // (K)
    arma::mat&    exogenous_forecast, // (horizon, d)
    arma::mat&    cond_forecast,      // (horizon, N)
    const int&    covid,
    const int&    T,
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
  cube        forecast_mean(N, horizon, S);
  cube        cov_s(N, N, horizon);
  field<cube> forecast_cov(S);

  for (int s=0; s<S; s++) {
    
    if ( do_exog ) {
      Xt          = join_cols(x_t, trans(exogenous_forecast.row(0)));
    } else {
      Xt          = x_t;
    } // END if do_exog
    
    // rescale the covariance matrix if covid is specified
    vec scale = ones<vec>(T+horizon);
    mat hyper = posterior_hyper.col(s);
    if (covid > 0 && covid <= T) {
      int c_idx = covid - 1;
      double s0 = hyper(N + 3);
      double s1 = hyper(N + 4);
      double s2 = hyper(N + 5);
      double rho = hyper(N + 6);

      if (c_idx < T)
        scale(c_idx) = s0;
      if (c_idx + 1 < T)
        scale(c_idx + 1) = s1;
      if (c_idx + 2 < T)
        scale(c_idx + 2) = s2;
      for (int t = c_idx + 3; t < T+horizon; t++)
      {
        scale(t) = 1.0 + (s2 - 1.0) * std::pow(rho, t - c_idx - 2);
      }
    }

    for (int h=0; h<horizon; h++) {
      
      mat   Sigma             = posterior_Sigma.slice(s);
      mat   mean              = posterior_A.slice(s) * Xt;
      mat   cov               = scale(T + h) * Sigma;
      vec   cond_forecast_h   = trans(cond_forecast.row(h));
      uvec  nonf_el           = find_nonfinite( cond_forecast_h );
      int   nonf_no           = nonf_el.n_elem;
      
      if ( nonf_no == N ) {
        forecasts.slice(s).col(h) = mvnrnd( mean, cov );
      } else {
        forecasts.slice(s).col(h) = bsvars::mvnrnd_cond( cond_forecast_h, mean, cov );   // does not work if cond_fc_h is all nan
      } // END if nonf_no
      
      if ( h != horizon - 1 ) {
        if ( do_exog ) {
          Xt          = join_cols( forecasts.slice(s).col(h), Xt.subvec(N, K - 1 - d), trans(exogenous_forecast.row(h + 1)) );
        } else {
          Xt          = join_cols( forecasts.slice(s).col(h), Xt.subvec(N, K - 1) );
        }
      } // END if h
      
      forecast_mean.slice(s).col(h) = mean;
      cov_s.slice(h) = cov;

    } // END h loop
    
    forecast_cov(s) = cov_s;

  } // END s loop
  
  return List::create(
    _["forecasts"]     = forecasts,
    _["forecast_mean"] = forecast_mean,
    _["forecast_cov"]  = forecast_cov
  );
} // END forecast_bsvarSIGNs

