
#include <functional>
#include <iostream>

#include <RcppArmadillo.h>

using namespace arma;


// QR decomposition, where the diagonal elements of R are positive
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
arma::mat qr_sign_cpp(const arma::mat& A) {
  int N = A.n_rows;
  
  arma::mat Q(N, N), R(N, N);
  arma::qr_econ(Q, R, A);

  // Check and modify the diagonal elements of R
  for(arma::uword i = 0; i < R.n_cols; ++i) {
    if(R(i, i) < 0) {
      // R.row(i) *= -1;  // Change sign of the column
      Q.col(i) *= -1;  // Change sign of the corresponding row in Q
    }
  }

  return Q;
}


// Sample uniformly from the space of NxN orthogonal matrices
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
arma::mat rortho_cpp(const int& N) {
  return qr_sign_cpp(arma::mat(N, N, fill::randn));
}


// If matches signs
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
bool match_sign(
    const arma::mat& A, 
    const arma::mat& sign
) {
  return accu(((A % sign) > 0)) == accu(sign != 0);
}


// Numerical derivative of f: R^n -> R^m, at x (Rcpp::export is not possible)
// [[Rcpp:interface(cpp)]]
arma::mat Df(
    const std::function<arma::vec(const arma::vec&)>& f,
    const arma::vec& x,
    const double     h = 1e-10
)
{
  colvec f_x = f(x);
  
  int n  = x.n_elem;
  int m  = f_x.n_elem;
  
  mat result(m, n);
  
  for (int i = 0; i < n; i++)
  {
    vec x_plus_h  = x;
    x_plus_h(i)  += h;
    
    vec f_plus_h = f(x_plus_h);
    
    for (int j = 0; j < m; j++)
    {
      result(j, i) = (f_plus_h(j) - f_x(j)) / h;
    }
  }
  
  return result;
}


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
  
  double s = 2.38 / sqrt(n);
  double d = log_target(x);
  
  double new_d, a;
  vec    new_x, xbar, diff;
  
  mat X(n, T);
  x        = log(x);
  X.col(0) = x;
  
  for (int t = 1; t < T; t++) {
    new_x = mvnrnd(x, s*s * Sigma);
    new_d = log_target(exp(new_x));
    a     = std::min(1.0, exp(new_d - d + sum(x - new_x)));
    
    if (randu() < a) {
      x = new_x;
      d = new_d;
    }
    X.col(t) = x;
    
    
    if (t == t0) {
      xbar   = mean(X.cols(0, t), 1);
      Sigma  = cov(X.cols(0, t).t());
    } else if (t > t0) {
      diff   = x - xbar;
      s     += pow(t, -0.6) * (a - 0.234);
      xbar  += diff / (t + 1);
      Sigma  = Sigma * t / (t + 1) + diff * diff.t() * t / (t + 1) / (t + 1);
    }
  }
  
  X = exp(X);
  return X;
}







