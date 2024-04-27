#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;


// // [[Rcpp::export]]
// arma::mat matnrnd_cpp_byA(const arma::mat& M,
//                       const arma::mat& U,
//                       const arma::mat& V) {
//   
//   mat X = mat(size(M), fill::randn);
//   return M + chol(U).t() * X * chol(V);
// }

// [[Rcpp::export]]
arma::cube rIW (
    const int         n, // a positive int
    const arma::mat&  S, // an NxN positive definite matrix
    const double&     nu // a positive double
) {
  int N         = S.n_cols;
  cube sigma(N, N, n);
  mat S_inv_chol = chol(inv_sympd(S));
  for (int i=0; i<n; i++) {
    sigma.slice(i)   = iwishrnd(S, nu, S_inv_chol);
  } // END i loop
  return sigma;
}


/*** R
n = 1000000
N = 2
# K = 3

set.seed(12)
# M = matrix(rnorm(N * K), K, N)
# U = rWishart(1, K + 1, diag(K))[,,1]
# V = rWishart(1, N + 1, diag(N))[,,1]
S = rWishart(1, N + 1, diag(N))[,,1]
nu = N + 2

# matnrnd_cpp_byA(M, U, V)
# set.seed(12)
# matrwish_byT1(S, nu)
# set.seed(12)
# matrwish_byT2(S, nu)

set.seed(12)
ssT1 = rIW(n, S, nu)
apply(ssT1, 1:2, mean)
apply(ssT1, 1:2, sd)
*/



