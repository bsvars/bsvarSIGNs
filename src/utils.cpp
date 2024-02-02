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
      R.col(i) *= -1;  // Change sign of the column
      Q.row(i) *= -1;  // Change sign of the corresponding row in Q
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
