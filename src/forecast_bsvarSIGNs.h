
#ifndef _FORECAST_BSVARSIGNS_H_
#define _FORECAST_BSVARSIGNS_H_

#include <RcppArmadillo.h>


arma::cube forecast_bsvarSIGNs (
    arma::cube&   posterior_Sigma,    // (N, N, S)
    arma::cube&   posterior_A,        // (N, K, S)
    arma::vec&    X_T,                // (K)
    arma::mat&    exogenous_forecast, // (horizon, d)
    arma::mat&    cond_forecast,      // (horizon, N)
    const int&    horizon
);


#endif  // _FORECAST_BSVARSIGNS_H_