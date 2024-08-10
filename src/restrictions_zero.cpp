
#include <RcppArmadillo.h>

#include "utils.h"
#include "compute.h"

using namespace Rcpp;
using namespace arma;


// Z_j * irf_0
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::mat> ZIRF(
    const arma::field<arma::mat>& Z,
    const arma::mat&              irf_0
) {
  
  arma::field<arma::mat> ZIRF(Z.n_elem);
  
  for (int j=0; j<Z.n_elem; j++) {
    ZIRF(j) = Z(j) * irf_0;
  }
  
  return ZIRF;
}


// Zero restrictions
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::colvec zero_restrictions(
    const arma::field<arma::mat>& Z,
    const arma::vec               vec_structural
) {
  int N  = Z(0).n_cols;
  
  Rcout << "debug 1" << std::endl;
  
  mat A0 = vec_structural.rows(0, N * N - 1);
  A0     = reshape(A0, N, N);
  
  Rcout << "debug 2" << std::endl;
  
  // field<mat> ZF = ZIRF(Z, inv(A0.t()));
  
  Rcout << "debug 3" << std::endl;
  
  vec z = vec();
  
  Rcout << "debug 4" << std::endl;
  
  mat irf_0 = inv(A0.t());
  
  for (int j=0; j<Z.n_elem; j++) {
    mat ZF = Z(j) * irf_0;
    Rcout << "debug 5" << std::endl;
    vec ZF_j = ZF.col(j);
    Rcout << "debug 6" << std::endl;
    z = join_vert(z, ZF_j);
    Rcout << "debug 7" << std::endl;
  }
  
  Rcout << "debug 8" << std::endl;
  
  return z;
}


// g ∘ f_h
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::colvec g_fh(
    const arma::field<arma::mat>& Z,
    const arma::mat& A0,
    const arma::mat& Aplus
) {
  mat B     = Aplus * inv(A0);
  mat Sigma = inv_sympd(A0 * A0.t());
  mat irf_0 = chol(Sigma, "lower");
  mat Q     = irf_0.t() * A0;
  
  arma::field<arma::mat> ZF = ZIRF(Z, irf_0);
  
  colvec out = join_vert(vectorise(B), vectorise(Sigma));
  
  mat    M_j, K_j;
  colvec w_j;
  for (int j=0; j<ZF.n_elem; j++) {
    if (j == 0) {
      M_j = ZF(j);
    } else {
      M_j = join_horiz(Q.cols(0, j-1), ZF(j).t()).t();  
    }
    
    K_j = null(M_j);
    w_j = K_j.t() * Q.col(j);
    out = join_vert(out, w_j);
  }
  
  return out;
}


// g ∘ f_h with vectorized input
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::colvec g_fh_vec(
    const arma::field<arma::mat>& Z,
    const arma::colvec            vec_structural
) {
  int N = Z(0).n_cols;
  int M = vec_structural.n_rows;
  int K = (M - N*N)/N;
  
  mat A0    = reshape(vec_structural.rows(0, N*N-1), N, N);
  mat Aplus = reshape(vec_structural.rows(N*N, M-1), K, N);
  
  return g_fh(Z, A0, Aplus);
}

// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat D_z(
    const arma::field<arma::mat>& Z,
    const arma::vec& x,
    const double     h = 1e-10
)
{
  vec f_x = zero_restrictions(Z, x);
  
  int n  = x.n_elem;
  int m  = f_x.n_elem;
  
  mat result(m, n);
  
  for (int i = 0; i < n; i++)
  {
    vec x_plus_h  = x;
    x_plus_h(i)  += h;
    
    vec f_plus_h = zero_restrictions(Z, x_plus_h);
    
    for (int j = 0; j < m; j++)
    {
      result(j, i) = (f_plus_h(j) - f_x(j)) / h;
    }
  }
  
  return result;
}


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat D_gf(
    const arma::field<arma::mat>& Z,
    const arma::vec& x,
    const double     h = 1e-10
)
{
  vec f_x = g_fh_vec(Z, x);
  
  int n  = x.n_elem;
  int m  = f_x.n_elem;
  
  mat result(m, n);
  
  for (int i = 0; i < n; i++)
  {
    vec x_plus_h  = x;
    x_plus_h(i)  += h;
    
    vec f_plus_h = g_fh_vec(Z, x_plus_h);
    
    for (int j = 0; j < m; j++)
    {
      result(j, i) = (f_plus_h(j) - f_x(j)) / h;
    }
  }
  
  return result;
}


// log volume element
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
double log_volume_element(
    const arma::field<arma::mat>& Z,
    const arma::mat&              A0,
    const arma::mat&              Aplus
) {
  colvec vec_structural = join_vert(vectorise(A0), vectorise(Aplus));
  
  mat Dz  = D_z(Z, vec_structural);
  mat Dgf = D_gf(Z, vec_structural);
  
  mat DN  = Dgf * null(Dz);
  
  return 0.5 * log_det(DN.t() * DN).real();
}


// importance weight for zero restrictions
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
double weight_zero(
    const arma::field<arma::mat>& Z,
    const arma::mat&              B,
    const arma::mat&              h_inv,
    const arma::mat&              Q
) {
  int K     = B.n_rows;
  int N     = Q.n_cols;
  
  mat A0    = h_inv * Q;
  mat Aplus = B * h_inv * Q;
  
  double log_ve_f   = -(2*N+K+1) * log_det(A0).real();
  double log_ve_gfz = log_volume_element(Z, A0, Aplus);
  
  return exp(log_ve_f - log_ve_gfz);
}


// draw Q conditional on zero restrictions
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat rzeroQ(
    const arma::field<arma::mat>& Z,
    const arma::mat&              irf_0
) {
  // Algorithm 2 in ARRW (2018)
  int N = irf_0.n_rows;
  
  int    z_j;
  colvec x_j, w_j, q_j;
  mat    Q, Z_j, M_j, K_j;
  
  arma::field<arma::mat> ZF = ZIRF(Z, irf_0);
  
  for (int j=0; j<Z.n_elem; j++) {
    z_j = Z(j).n_rows;
    
    x_j = colvec(N-j-z_j, fill::randn);
    w_j = x_j / sqrt(sum(square(x_j)));
    
    M_j = join_horiz(Q, (ZF(j).t())).t();
    K_j = null(M_j);
    q_j = K_j * w_j;
    Q   = join_horiz(Q, q_j);
  }
  
  return Q;
}