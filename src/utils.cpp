#include <RcppArmadillo.h>

#include <bsvars.h>

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


// Compute impulse response functions with bsvars
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
arma::field<arma::cube> irf_cpp(
    arma::cube&         posterior_B,        // (N, N, S)
    arma::cube&         posterior_A,        // (N, K, S)
    const int           horizon,
    const int           p
) {
  
  // const int   S = posterior_Q.n_slices;
  // 
  // cube        QB(size(posterior_Q));
  // 
  // for(int s = 0; s < S; s++) {
  //   QB.slice(s) = posterior_Q.slice(s).t() * posterior_B.slice(s);
  // }
  
  return bsvars::bsvars_ir(posterior_B, posterior_A, horizon, p);
}









