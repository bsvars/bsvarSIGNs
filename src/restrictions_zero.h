#ifndef _RESTRICTIONS_ZERO_H_
#define _RESTRICTIONS_ZERO_H_

#include <RcppArmadillo.h>

#include "utils.h"
#include "bsvarTOOLs.h"

using namespace Rcpp;
using namespace arma;


double weight_zero(
    const arma::field<arma::mat>& Z,
    const arma::mat&              B,
    const arma::mat&              h_inv,
    const arma::mat&              Q
);

arma::mat rzeroQ(
    const arma::field<arma::mat>& Z,
    const arma::mat&              irf_0
);

#endif  // _RESTRICTIONS_ZERO_H_