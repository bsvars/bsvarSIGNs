#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
bool match_sign_cpp(const arma::mat& A, const arma::mat sign) {
  return accu(((A % sign) > 0)) == accu(abs(sign));
}

/*** R
A = matrix(rnorm(4),2,2); A
sign = matrix(c(0,0,-1,-1),2,2)
match_sign_cpp(A, sign)

*/
