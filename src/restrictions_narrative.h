#ifndef _RESTRICTIONS_NARRATIVE_H_
#define _RESTRICTIONS_NARRATIVE_H_

#include <RcppArmadillo.h>

#include "utils.h"
#include "bsvarTOOLs.h"

using namespace Rcpp;
using namespace arma;


bool match_sign_narrative(
    const arma::mat&  Epsilon,
    const arma::mat&  sign_narrative,
    const arma::cube& irf
);


double weight_narrative(
    const int&                    T,
    arma::mat                     sign_narrative,
    const arma::cube&             irf
);

#endif  // _RESTRICTIONS_NARRATIVE_H_