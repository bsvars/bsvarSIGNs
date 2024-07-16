
#include <RcppArmadillo.h>

using namespace arma;


// compute historical decomposition from t to t+h of the i-th variable
// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat hd1_cpp(
    const int&        var_i,   // i-th variable
    const int&        t,       // start at period t
    const int&        h,       // number of horizons
    const arma::mat&  Epsilon, // structural shocks, NxT
    const arma::cube& irf
) {
  
  int ii = var_i - 1;
  int tt = t - 1;
  int N  = Epsilon.n_rows;
  
  mat hd(N, h + 1);
  
  // for each shock jj
  for(int jj = 0; jj < N; jj++) {
    // for each horizon hh
    for(int hh = 0; hh <= h; hh++) {
      // for each lag ll
      for(int ll = 0; ll <= hh; ll++) {
        hd(jj, hh) += irf(ii, jj, ll) * Epsilon(jj, tt + hh - ll);
      }
    }
  }
  return hd;
} // END hd1_cpp

