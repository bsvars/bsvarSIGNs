
#include <RcppArmadillo.h>

using namespace arma;


// compute historical decomposition from t to t+h of the i-th variable
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat hd1_cpp(
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
} // END hd1_cpp


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::cube ir1_cpp(
    const arma::mat& B,           // KxN
    const arma::mat& Theta0,      // NxN
    int              horizon,
    const int&       p
) {
  
  horizon++;
  int N         = Theta0.n_cols;
  
  cube irf      = zeros(N, N, horizon);
  irf.slice(0)  = Theta0;
  
  for(int t = 2; t <= horizon; t++) {
    for(int j = 1; j <= std::min(p, t - 1); j++) {
      mat B_j           = B.rows((j - 1) * N, j * N - 1).t();
      irf.slice(t - 1) += B_j * irf.slice(t - j - 1);
    } // END j loop
  } // END t loop
  
  return irf;
} // END ir1_cpp