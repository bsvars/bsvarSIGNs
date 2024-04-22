
#include <RcppArmadillo.h>

using namespace arma;


// compute historical decomposition from t to t+h of the i-th variable
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat hd1(
    const int&        var_i,  // i-th variable
    const int&        t,      // start at period t
    const int&        h,      // number of horizons
    const arma::mat&  U,      // structural shocks 
    const arma::cube& irf
) {
  
  int ii = var_i - 1;
  int tt = t - 1;
  int N  = U.n_rows;
  
  mat hd(N, h + 1);
  
  // for each shock jj
  for(int jj = 0; jj < N; jj++) {
    // for each horizon hh
    for(int hh = 0; hh <= h; hh++) {
      // for each lag ll
      for(int ll = 0; ll <= hh; ll++) {
        hd(jj, hh) += irf(ii, jj, ll) * U(jj, tt + hh - ll);
      }
    }
  }
  return hd;
}


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube ir11_cpp (
    arma::mat&    irf_0,              // (N, N)
    arma::mat&    aux_A,              // (N, K)
    const int     horizon,
    const int     p,
    const bool    standardise = false
) {
  
  const int       N = irf_0.n_rows;
  cube            aux_irfs(N, N, horizon + 1);  // + 0 horizons
  mat             A_bold_tmp(N * (p - 1), N * p, fill::eye);
  
  if ( standardise ) {
    irf_0             = irf_0 * diagmat(pow(diagvec(irf_0), -1));
  }
  mat   A_bold        = join_cols(aux_A.cols(0, N * p - 1), A_bold_tmp);
  mat   A_bold_power  = A_bold;
  
  aux_irfs.slice(0)   = irf_0;
  
  for (int h=1; h<horizon + 1; h++) {
    aux_irfs.slice(h) = A_bold_power.submat(0, 0, N-1, N-1) * irf_0;
    A_bold_power      = A_bold_power * A_bold;
  } // END h loop
  
  return aux_irfs;
} // END bsvars_ir1


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube ir1_cpp(
    const arma::mat& A, 
    const arma::mat& chol_SIGMA, 
    int              horizon,
    const int&       p
) {
  
  horizon++;
  int N = chol_SIGMA.n_cols;
  
  cube irf = zeros(N, N, horizon);
  irf.slice(0) = chol_SIGMA;
  
  for(int t = 2; t <= horizon; t++) {
    for(int j = 1; j <= std::min(p, t - 1); j++) {
      mat A_j = A.rows(1 + (j - 1) * N, j * N).t();
      irf.slice(t - 1) += A_j * irf.slice(t - j - 1);
    }
  }
  
  return irf;
}