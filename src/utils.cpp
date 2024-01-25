#include <RcppArmadillo.h>

using namespace arma;

// log determinant given Cholesky decomposition of some matrix
double log_det_lower_cpp(const arma::mat& L) {
  return 2 * sum(log(diagvec(L)));
}

// QR decomposition, where the diagonal elements of R are positive
// [[Rcpp:interface(cpp)]]
arma::mat qr_sign_cpp(const arma::mat& A) {
  arma::mat Q, R;
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
