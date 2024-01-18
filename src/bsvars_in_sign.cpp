
#include <RcppArmadillo.h>
#include "progress.hpp"
#include "Rcpp/Rmath.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::mat orthogonal_complement_matrix_TW (const arma::mat& x) {
  // # x is a mxn matrix and m>n
  // # the function returns a mx(m-n) matrix, out, that is an orthogonal complement of x, i.e.:
  // # t(x)%*%out = 0 and det(cbind(x,out))!=0
  int n_nrow     = x.n_rows;
  int n_ncol     = x.n_cols;
  mat Q;
  mat R;
  qr(Q, R, x);
  mat ocm = Q.tail_cols(n_nrow-n_ncol);
  return ocm;
} // END orthogonal_complement_matrix_TW


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
std::string ordinal(
    int n
) {
  std::string suffix;
  if (n % 10 == 1 && n % 100 != 11) {
    suffix = "st";
  } else if (n % 10 == 2 && n % 100 != 12) {
    suffix = "nd";
  } else if (n % 10 == 3 && n % 100 != 13) {
    suffix = "rd";
  } else {
    suffix = "th";
  }
  return std::to_string(n) + suffix;
} // END ordinal


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
void sample_hyperparameters (
    arma::mat&              aux_hyper,       // (2*N+1) x 2 :: col 0 for B, col 1 for A
    const arma::mat&        aux_B,            // NxN
    const arma::mat&        aux_A,
    const arma::field<arma::mat>& VB,
    const Rcpp::List&       prior
) {
  // the function returns aux_hyper by reference (filling it with a new draw)
  
  const int N = aux_B.n_rows;
  const int K = aux_A.n_cols;
  
  double prior_hyper_nu_B     = as<double>(prior["hyper_nu_B"]);
  double prior_hyper_a_B      = as<double>(prior["hyper_a_B"]);
  double prior_hyper_s_BB     = as<double>(prior["hyper_s_BB"]);
  double prior_hyper_nu_BB    = as<double>(prior["hyper_nu_BB"]);
  
  double prior_hyper_nu_A     = as<double>(prior["hyper_nu_A"]);
  double prior_hyper_a_A      = as<double>(prior["hyper_a_A"]);
  double prior_hyper_s_AA     = as<double>(prior["hyper_s_AA"]);
  double prior_hyper_nu_AA    = as<double>(prior["hyper_nu_AA"]);
  
  mat   prior_A               = as<mat>(prior["A"]);
  mat   prior_A_V_inv         = as<mat>(prior["A_V_inv"]);
  mat   prior_B_V_inv         = as<mat>(prior["B_V_inv"]);
  
  // aux_B - related hyper-parameters 
  vec     ss_tmp      = aux_hyper.submat(N, 0, 2 * N - 1, 0);
  double  scale_tmp   = prior_hyper_s_BB + 2 * sum(ss_tmp);
  double  shape_tmp   = prior_hyper_nu_BB + 2 * N * prior_hyper_a_B;
  aux_hyper(2 * N, 0) = scale_tmp / R::rchisq(shape_tmp);
  
  // aux_A - related hyper-parameters 
  ss_tmp              = aux_hyper.submat(N, 1, 2 * N - 1, 1);
  scale_tmp           = prior_hyper_s_AA + 2 * sum(ss_tmp);
  shape_tmp           = prior_hyper_nu_AA + 2 * N * prior_hyper_a_A;
  aux_hyper(2 * N, 1) = scale_tmp / R::rchisq(shape_tmp);
  
  for (int n=0; n<N; n++) {
    
    // count unrestricted elements of aux_B's row
    int rn            = VB(n).n_rows;
    
    // aux_B - related hyper-parameters 
    scale_tmp         = 1 / ((1 / (2 * aux_hyper(n, 0))) + (1 / aux_hyper(2 * N, 0)));
    shape_tmp         = prior_hyper_a_B + prior_hyper_nu_B / 2;
    aux_hyper(N + n, 0) = R::rgamma(shape_tmp, scale_tmp);
    
    scale_tmp         = aux_hyper(N + n, 0) + as_scalar(aux_B.row(n) * prior_B_V_inv * aux_B.row(n).t());
    shape_tmp         = prior_hyper_nu_B + rn;
    aux_hyper(n, 0)   = scale_tmp / R::rchisq(shape_tmp);
    
    // aux_A - related hyper-parameters 
    scale_tmp         = 1 / ((1 / (2 * aux_hyper(n, 1))) + (1 / aux_hyper(2 * N, 1)));
    shape_tmp         = prior_hyper_a_A + prior_hyper_nu_A / 2;
    aux_hyper(N + n, 1) = R::rgamma(shape_tmp, scale_tmp);
    
    scale_tmp         = aux_hyper(N + n, 1) + 
      as_scalar((aux_A.row(n) - prior_A.row(n)) * prior_A_V_inv * trans(aux_A.row(n) - prior_A.row(n)));
    shape_tmp         = prior_hyper_nu_A + K;
    aux_hyper(n, 1)   = scale_tmp / R::rchisq(shape_tmp);
  } // END n loop
} // END sample_hyperparameters


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
void sample_A_homosk1 (
    arma::mat&        aux_A,          // NxK
    const arma::mat&  aux_B,          // NxN
    const arma::mat&  aux_hyper,      // (2*N+1) x 2 :: col 0 for B, col 1 for A
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior           // a list of priors - original dimensions
) {
  // the function changes the value of aux_A by reference
  const int N         = aux_A.n_rows;
  const int K         = aux_A.n_cols;
  
  mat prior_A_mean    = as<mat>(prior["A"]);
  mat prior_A_Vinv    = as<mat>(prior["A_V_inv"]);
  rowvec    zerosA(K);
  
  for (int n=0; n<N; n++) {
    mat   A0          = aux_A;
    A0.row(n)         = zerosA;
    vec   zn          = vectorise( aux_B * (Y - A0 * X) );
    mat   Wn          = kron( trans(X), aux_B.col(n) );
    
    mat     precision = (pow(aux_hyper(n,1), -1) * prior_A_Vinv) + trans(Wn) * Wn;
    rowvec  location  = prior_A_mean.row(n) * (pow(aux_hyper(n,1), -1) * prior_A_Vinv) + trans(zn) * Wn;
    
    mat     precision_chol = trimatu(chol(precision));
    vec     xx(K, fill::randn);
    vec     draw      = solve(precision_chol, 
                              solve(trans(precision_chol), trans(location)) + xx);
    aux_A.row(n)      = trans(draw);
  } // END n loop
} // END sample_A_homosk1


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
void sample_B_homosk1 (
    arma::mat&        aux_B,          // NxN
    const arma::mat&  aux_A,          // NxK
    const arma::mat&  aux_hyper,      // (2*N+1) x 2 :: col 0 for B, col 1 for A
    const arma::mat&  Y,              // NxT dependent variables
    const arma::mat&  X,              // KxT dependent variables
    const Rcpp::List& prior,          // a list of priors - original dimensions
    const arma::field<arma::mat>& VB        // restrictions on B0
) {
  // the function changes the value of aux_B by reference
  const int N               = aux_B.n_rows;
  const int T               = Y.n_cols;
  
  const int posterior_nu    = T + as<int>(prior["B_nu"]);
  mat prior_SS_inv          = as<mat>(prior["B_V_inv"]);
  mat shocks                = Y - aux_A * X;
  mat posterior_SS_inv      = shocks * shocks.t();
  
  for (int n=0; n<N; n++) {
    mat posterior_S_inv     = VB(n) * ((pow(aux_hyper(n,0), -1) * prior_SS_inv) + posterior_SS_inv) * VB(n).t();
    posterior_S_inv         = 0.5*( posterior_S_inv + posterior_S_inv.t() );
    
    // sample B
    mat Un                  = chol(posterior_nu * inv_sympd(posterior_S_inv));
    mat B_tmp               = aux_B;
    B_tmp.shed_row(n);
    rowvec w                = trans(orthogonal_complement_matrix_TW(B_tmp.t()));
    vec w1_tmp              = trans(w * VB(n).t() * Un.t());
    double w1w1_tmp         = as_scalar(sum(pow(w1_tmp, 2)));
    mat w1                  = w1_tmp.t()/sqrt(w1w1_tmp);
    mat Wn;
    const int rn            = VB(n).n_rows;
    if (rn==1) {
      Wn                    = w1;
    } else {
      Wn                    = join_rows(w1.t(), orthogonal_complement_matrix_TW(w1.t()));
    }
    
    vec   alpha(rn);
    vec   u(posterior_nu + 1, fill::randn);
    u                      *= pow(posterior_nu, -0.5);
    alpha(0)                = sqrt(as_scalar(sum(square(u))));
    if (R::runif(0,1)<0.5) {
      alpha(0)       *= -1;
    }
    if (rn>1){
      vec nn(rn-1, fill::randn);
      nn                   *= pow(posterior_nu, -0.5);
      alpha.rows(1,rn-1)    = nn;
    }
    rowvec b0n              = alpha.t() * Wn * Un;
    aux_B.row(n)            = b0n * VB(n);
  } // END n loop
} // END sample_B_homosk1


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List bsvar_sign_cpp(
    const int&  S,                        // number of draws from the posterior
    const arma::mat&  Y,                  // NxT dependent variables
    const arma::mat&  X,                  // KxT dependent variables
    const arma::field<arma::mat>& VB,     // N-list
    const Rcpp::List& prior,              // a list of priors
    const Rcpp::List& starting_values,    // a list of starting values
    const int         thin = 100,         // introduce thinning
    const bool        show_progress = true
) {
  
  std::string oo = "";
  if ( thin != 1 ) {
    oo      = ordinal(thin) + " ";
  }
  
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, S, 50));
  if (show_progress) {
    Rcout << "**************************************************|" << endl;
    Rcout << "bsvars: Bayesian Structural Vector Autoregressions|" << endl;
    Rcout << "**************************************************|" << endl;
    Rcout << " Gibbs sampler for the SVAR model                 |" << endl;
    Rcout << "**************************************************|" << endl;
    Rcout << " Progress of the MCMC simulation for " << S << " draws" << endl;
    Rcout << "    Every " << oo << "draw is saved via MCMC thinning" << endl;
    Rcout << " Press Esc to interrupt the computations" << endl;
    Rcout << "**************************************************|" << endl;
  }
  Progress p(50, show_progress);
  
  const int N       = Y.n_rows;
  const int K       = X.n_rows;
  
  mat   aux_B       = as<mat>(starting_values["B"]);
  mat   aux_A       = as<mat>(starting_values["A"]);
  mat   aux_hyper   = as<mat>(starting_values["hyper"]);
  
  const int   SS    = floor(S / thin);
  
  cube  posterior_B(N, N, SS);
  cube  posterior_A(N, K, SS);
  cube  posterior_hyper(2 * N + 1, 2, SS);
  
  int   ss = 0;
  
  for (int s=0; s<S; s++) {
    
    // Increment progress bar
    if (any(prog_rep_points == s)) p.increment();
    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    sample_hyperparameters(aux_hyper, aux_B, aux_A, VB, prior);
    sample_A_homosk1(aux_A, aux_B, aux_hyper, Y, X, prior);
    sample_B_homosk1(aux_B, aux_A, aux_hyper, Y, X, prior, VB);
    
    if (s % thin == 0) {
      posterior_B.slice(ss)    = aux_B;
      posterior_A.slice(ss)    = aux_A;
      posterior_hyper.slice(ss)  = aux_hyper;
      ss++;
    }
  } // END s loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["B"]        = aux_B,
      _["A"]        = aux_A,
      _["hyper"]    = aux_hyper
    ),
    _["posterior"]  = List::create(
      _["B"]        = posterior_B,
      _["A"]        = posterior_A,
      _["hyper"]    = posterior_hyper
    )
  );
} // END bsvar_sign_cpp
