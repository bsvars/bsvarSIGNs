
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
    const arma::vec& x
)
{
  
  double h   = 1e-10;
  
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

