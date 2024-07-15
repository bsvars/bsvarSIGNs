
#include <functional>
#include <iostream>
#include <RcppArmadillo.h>
#include "Rcpp/Rmath.h"
#include "progress.hpp"
#include "ramcmc.h"

using namespace arma;


// adaptive Metropolis algorithm for strictly positive parameters
// [[Rcpp:interfaces(cpp)]]
arma::mat metropolis(
    const int& T,
    const int& t0,
    arma::vec  x,
    arma::mat  Sigma,
    const std::function<double(const arma::vec&)>& log_target
) {
  
  int    n = x.n_elem;
  
  x        = log(x);
  vec xbar = x;
  double s = 2.38 / sqrt(n);
  double d = log_target(x) + sum(x);
  
  double new_d, a;
  vec    new_x, diff;
  
  mat X(n, T);
  X.col(0) = x;
  
  vec      progress = round(linspace(0, T, 50));
  Progress p(50, true);
  
  for (int t = 1; t < T; t++) {
    new_x = mvnrnd(x, s*s * Sigma);
    new_d = log_target(exp(new_x)) + sum(new_x);
    a     = std::min(1.0, exp(new_d - d));
    
    if (randu() < a) {
      x = new_x;
      d = new_d;
    }
    X.col(t) = x;
    
    if (t <= t0) {
      diff   = x - xbar;
      s     += pow(t, -0.6) * (a - 0.234);
      xbar  += diff / (t + 1);
      Sigma  = Sigma * t / (t + 1) + diff * diff.t() * t / (t + 1) / (t + 1);
    }
    
    if (any(progress == t)) {
      p.increment();
    }
  }
  
  X = exp(X);
  return X;
}







