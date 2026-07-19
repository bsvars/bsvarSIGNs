#include "sample_SOE.h"
#include "compute.h"
#include "sample_NIW.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::field<arma::mat> sample_restricted_B_cpp(
    const arma::mat& post_B,
    const arma::mat& post_V,
    const arma::mat& Sigma,
    const int& p,
    const int& N,
    const int& Nf,
    const int& K
) {
    int d = K - p * N - 1;
    
    arma::uvec R2(p * (N - Nf));
    arma::uvec R1(p * Nf + 1 + d);
    
    int idx1 = 0;
    int idx2 = 0;
    for (int l = 0; l < p; l++) {
        for (int i = 0; i < Nf; i++) R1(idx1++) = l * N + i;
        for (int i = Nf; i < N; i++) R2(idx2++) = l * N + i;
    }
    R1(idx1++) = p * N;
    for (int i = 0; i < d; i++) R1(idx1++) = p * N + 1 + i;
    
    arma::mat V11 = post_V.submat(R1, R1);
    arma::mat V12 = post_V.submat(R1, R2);
    arma::mat V21 = post_V.submat(R2, R1);
    arma::mat V22 = post_V.submat(R2, R2);
    
    arma::mat V22_inv = arma::inv_sympd(V22);
    
    arma::mat post_B_f = post_B.cols(0, Nf - 1);
    arma::mat M1 = post_B_f.rows(R1);
    arma::mat M2 = post_B_f.rows(R2);
    
    arma::mat M1_cond = M1 - V12 * V22_inv * M2;
    arma::mat V11_cond = arma::symmatu(V11 - V12 * V22_inv * V21);
    
    arma::mat Sigma11 = arma::symmatu(Sigma(arma::span(0, Nf - 1), arma::span(0, Nf - 1)));
    
    arma::mat B_f_R1 = rmatnorm_cpp(M1_cond, V11_cond, Sigma11);
    
    arma::mat Bf = arma::zeros(K, Nf);
    Bf.rows(R1) = B_f_R1;
    // Bf.rows(R2) are already zero
    
    arma::mat Md = post_B.cols(Nf, N - 1);
    arma::mat Mf = post_B.cols(0, Nf - 1);
    
    arma::mat Sigma12 = Sigma(arma::span(0, Nf - 1), arma::span(Nf, N - 1));
    arma::mat Sigma21 = Sigma(arma::span(Nf, N - 1), arma::span(0, Nf - 1));
    arma::mat Sigma22 = Sigma(arma::span(Nf, N - 1), arma::span(Nf, N - 1));
    
    arma::mat Sigma11_inv = arma::inv_sympd(Sigma11);
    
    arma::mat Md_cond = Md + (Bf - Mf) * Sigma11_inv * Sigma12;
    arma::mat Sigma22_cond = arma::symmatu(Sigma22 - Sigma21 * Sigma11_inv * Sigma12);
    
    arma::mat Bd = rmatnorm_cpp(Md_cond, post_V, Sigma22_cond);
    
    arma::mat B = arma::join_rows(Bf, Bd);
    
    // Calculate log density of b2 = 0
    double val;
    double sign;
    arma::log_det(val, sign, V22);
    double log_det_V22 = val; // Assuming V22 is pos def, sign is positive
    
    arma::log_det(val, sign, Sigma11);
    double log_det_Sigma11 = val;
    
    double trace_term = arma::trace(Sigma11_inv * M2.t() * V22_inv * M2);
    
    double log_w_B = -0.5 * ( R2.n_elem * Nf * std::log(2 * M_PI) + 
                              Nf * log_det_V22 + 
                              R2.n_elem * log_det_Sigma11 + 
                              trace_term );
    
    arma::field<arma::mat> result(2);
    result(0) = B;
    result(1) = arma::mat(1, 1, arma::fill::zeros);
    result(1)(0, 0) = log_w_B;
    
    return result;
}
