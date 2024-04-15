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
arma::mat matrwish_byT1 (
    const arma::mat&  S, 
    const double&     nu
) {
  mat sigma   = iwishrnd(S, nu);
  return sigma;
}

// [[Rcpp::export]]
arma::mat matrwish_byT2 (
    const arma::mat&  S, 
    const double&     nu
) {
  // Based on algorithm B.4.4. from Appendinx B by Bauwens, Lubrano, Richard (1999) Bayesian Inference in Dynamic Econometric Models, Oxford Uni Press
  
  int N           = S.n_cols;

  mat s_chol      = chol(S, "lower");

  mat Q(N, N, fill::zeros);
  Q.diag()        = sqrt(pow(chi2rnd(nu - N + 1, N), -1));
  
  for (int i = 0; i < (N - 1); i++) {
    Q.submat(i + 1, i, N - 1, i) = randn(N - i - 1);
  }
  mat Q_inv       = inv(trimatu(Q));
  
  return s_chol * Q_inv.t() * Q_inv * s_chol.t();
}



/*** R
N = 20
# K = 3

set.seed(12)
# M = matrix(rnorm(N * K), K, N)
# U = rWishart(1, K + 1, diag(K))[,,1]
# V = rWishart(1, N + 1, diag(N))[,,1]
S = rWishart(1, N + 1, diag(N))[,,1]
nu = N + 1

# matnrnd_cpp_byA(M, U, V)
# set.seed(12)
# matrwish_byT1(S, nu)
# set.seed(12)
# matrwish_byT2(S, nu)

# microbenchmark::microbenchmark(
#   matrwish_byT1(S, nu),
#   matrwish_byT2(S, nu),
#   times = 100000
# )  

rbenchmark::benchmark(
  matrwish_byT1(S, nu),
  matrwish_byT2(S, nu),
  replications = 100000
)
  
*/



