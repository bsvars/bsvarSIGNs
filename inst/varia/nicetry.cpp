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
arma::cube rIW_byT1 (
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

// [[Rcpp::export]]
arma::cube rIW_byT1b (
    const int         n, // a positive int
    const arma::mat&  S, // an NxN positive definite matrix
    const double&     nu // a positive double
) {
  int N         = S.n_cols;
  cube sigma(N, N, n);
  for (int i=0; i<n; i++) {
    sigma.slice(i)   = iwishrnd(S, nu);
  } // END i loop
  return sigma;
}



// [[Rcpp::export]]
arma::cube rIW_byT2 (
    const int         n, // a positive int
    const arma::mat&  S, // an NxN positive definite matrix
    const double&     nu // a positive double
) {
  // Based on algorithm B.4.4. from Appendinx B by Bauwens, Lubrano, Richard (1999) Bayesian Inference in Dynamic Econometric Models, Oxford Uni Press
  
  int N           = S.n_cols;
  mat s_chol      = chol(S, "lower");
  cube sigma(N, N, n);

  for (int i=0; i<n; i++) {
    mat Q(N, N, fill::zeros);
    Q.diag()        = sqrt(pow(chi2rnd(nu - N + 1, N), -1));
    
    for (int i = 0; i < (N - 1); i++) {
      Q.submat(i + 1, i, N - 1, i) = randn(N - i - 1);
    }
    mat Q_inv       = inv(trimatu(Q));
    sigma.slice(i)  = s_chol * Q_inv.t() * Q_inv * s_chol.t();
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

# rbenchmark::benchmark(
#   rIW_byT1(n, S, nu),
#   rIW_byT2(n, S, nu),
#   replications = 1000
# )

set.seed(12)
ssT1 = rIW_byT1(n, S, nu)
# set.seed(12)
ssT1b = rIW_byT1b(n, S, nu)
# set.seed(12)
ssT2 = rIW_byT2(n, S, nu)

apply(ssT1, 1:2, mean)
apply(ssT1b, 1:2, mean)  
apply(ssT2, 1:2, mean)  

apply(ssT1, 1:2, sd)
apply(ssT1b, 1:2, sd)  
apply(ssT2, 1:2, sd)  
*/



