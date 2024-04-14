// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/bsvarSIGNs.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bsvar_sign_cpp
Rcpp::List bsvar_sign_cpp(const int& S, const int& lags, const arma::mat& Y, const arma::mat& X, const arma::field<arma::mat>& VB, const arma::cube& sign_irf, const arma::mat& sign_narrative, const arma::mat& sign_B, const Rcpp::List& prior, const Rcpp::List& starting_values, const int thin, const bool show_progress, const int& max_tries);
static SEXP _bsvarSIGNs_bsvar_sign_cpp_try(SEXP SSEXP, SEXP lagsSEXP, SEXP YSEXP, SEXP XSEXP, SEXP VBSEXP, SEXP sign_irfSEXP, SEXP sign_narrativeSEXP, SEXP sign_BSEXP, SEXP priorSEXP, SEXP starting_valuesSEXP, SEXP thinSEXP, SEXP show_progressSEXP, SEXP max_triesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const int& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const int& >::type lags(lagsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::mat>& >::type VB(VBSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type sign_irf(sign_irfSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sign_narrative(sign_narrativeSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sign_B(sign_BSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type starting_values(starting_valuesSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const bool >::type show_progress(show_progressSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_tries(max_triesSEXP);
    rcpp_result_gen = Rcpp::wrap(bsvar_sign_cpp(S, lags, Y, X, VB, sign_irf, sign_narrative, sign_B, prior, starting_values, thin, show_progress, max_tries));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _bsvarSIGNs_bsvar_sign_cpp(SEXP SSEXP, SEXP lagsSEXP, SEXP YSEXP, SEXP XSEXP, SEXP VBSEXP, SEXP sign_irfSEXP, SEXP sign_narrativeSEXP, SEXP sign_BSEXP, SEXP priorSEXP, SEXP starting_valuesSEXP, SEXP thinSEXP, SEXP show_progressSEXP, SEXP max_triesSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_bsvarSIGNs_bsvar_sign_cpp_try(SSEXP, lagsSEXP, YSEXP, XSEXP, VBSEXP, sign_irfSEXP, sign_narrativeSEXP, sign_BSEXP, priorSEXP, starting_valuesSEXP, thinSEXP, show_progressSEXP, max_triesSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// bsvar_sign_gibbs_cpp
Rcpp::List bsvar_sign_gibbs_cpp(const int& S, const int& lags, const arma::mat& Y, const arma::mat& X, const arma::field<arma::mat>& VB, const arma::cube& sign_irf, const arma::mat& sign_narrative, const arma::mat& sign_B, const Rcpp::List& prior, const Rcpp::List& starting_values, const int thin, const bool show_progress, const int& max_tries);
static SEXP _bsvarSIGNs_bsvar_sign_gibbs_cpp_try(SEXP SSEXP, SEXP lagsSEXP, SEXP YSEXP, SEXP XSEXP, SEXP VBSEXP, SEXP sign_irfSEXP, SEXP sign_narrativeSEXP, SEXP sign_BSEXP, SEXP priorSEXP, SEXP starting_valuesSEXP, SEXP thinSEXP, SEXP show_progressSEXP, SEXP max_triesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const int& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const int& >::type lags(lagsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::mat>& >::type VB(VBSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type sign_irf(sign_irfSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sign_narrative(sign_narrativeSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sign_B(sign_BSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type starting_values(starting_valuesSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const bool >::type show_progress(show_progressSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_tries(max_triesSEXP);
    rcpp_result_gen = Rcpp::wrap(bsvar_sign_gibbs_cpp(S, lags, Y, X, VB, sign_irf, sign_narrative, sign_B, prior, starting_values, thin, show_progress, max_tries));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _bsvarSIGNs_bsvar_sign_gibbs_cpp(SEXP SSEXP, SEXP lagsSEXP, SEXP YSEXP, SEXP XSEXP, SEXP VBSEXP, SEXP sign_irfSEXP, SEXP sign_narrativeSEXP, SEXP sign_BSEXP, SEXP priorSEXP, SEXP starting_valuesSEXP, SEXP thinSEXP, SEXP show_progressSEXP, SEXP max_triesSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_bsvarSIGNs_bsvar_sign_gibbs_cpp_try(SSEXP, lagsSEXP, YSEXP, XSEXP, VBSEXP, sign_irfSEXP, sign_narrativeSEXP, sign_BSEXP, priorSEXP, starting_valuesSEXP, thinSEXP, show_progressSEXP, max_triesSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// niw_cpp
void niw_cpp(arma::mat& post_A, arma::mat& post_V, arma::mat& post_S, int& post_nu, const arma::mat& Y, const arma::mat& X, const Rcpp::List prior);
RcppExport SEXP _bsvarSIGNs_niw_cpp(SEXP post_ASEXP, SEXP post_VSEXP, SEXP post_SSEXP, SEXP post_nuSEXP, SEXP YSEXP, SEXP XSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type post_A(post_ASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type post_V(post_VSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type post_S(post_SSEXP);
    Rcpp::traits::input_parameter< int& >::type post_nu(post_nuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type prior(priorSEXP);
    niw_cpp(post_A, post_V, post_S, post_nu, Y, X, prior);
    return R_NilValue;
END_RCPP
}
// match_sign
bool match_sign(const arma::mat& A, const arma::mat& sign);
static SEXP _bsvarSIGNs_match_sign_try(SEXP ASEXP, SEXP signSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sign(signSEXP);
    rcpp_result_gen = Rcpp::wrap(match_sign(A, sign));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _bsvarSIGNs_match_sign(SEXP ASEXP, SEXP signSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_bsvarSIGNs_match_sign_try(ASEXP, signSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// match_sign_irf
bool match_sign_irf(const arma::mat& Q, const arma::cube& sign_irf, const arma::cube& irf);
static SEXP _bsvarSIGNs_match_sign_irf_try(SEXP QSEXP, SEXP sign_irfSEXP, SEXP irfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type sign_irf(sign_irfSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type irf(irfSEXP);
    rcpp_result_gen = Rcpp::wrap(match_sign_irf(Q, sign_irf, irf));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _bsvarSIGNs_match_sign_irf(SEXP QSEXP, SEXP sign_irfSEXP, SEXP irfSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_bsvarSIGNs_match_sign_irf_try(QSEXP, sign_irfSEXP, irfSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// hd1
arma::mat hd1(const int& var_i, const int& t, const int& h, const arma::mat& U, const arma::cube& irf);
static SEXP _bsvarSIGNs_hd1_try(SEXP var_iSEXP, SEXP tSEXP, SEXP hSEXP, SEXP USEXP, SEXP irfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const int& >::type var_i(var_iSEXP);
    Rcpp::traits::input_parameter< const int& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const int& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type irf(irfSEXP);
    rcpp_result_gen = Rcpp::wrap(hd1(var_i, t, h, U, irf));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _bsvarSIGNs_hd1(SEXP var_iSEXP, SEXP tSEXP, SEXP hSEXP, SEXP USEXP, SEXP irfSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_bsvarSIGNs_hd1_try(var_iSEXP, tSEXP, hSEXP, USEXP, irfSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// match_sign_narrative
bool match_sign_narrative(const arma::mat& U, const arma::mat& sign_narrative, const arma::cube& irf);
static SEXP _bsvarSIGNs_match_sign_narrative_try(SEXP USEXP, SEXP sign_narrativeSEXP, SEXP irfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sign_narrative(sign_narrativeSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type irf(irfSEXP);
    rcpp_result_gen = Rcpp::wrap(match_sign_narrative(U, sign_narrative, irf));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _bsvarSIGNs_match_sign_narrative(SEXP USEXP, SEXP sign_narrativeSEXP, SEXP irfSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_bsvarSIGNs_match_sign_narrative_try(USEXP, sign_narrativeSEXP, irfSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// approximate_w
double approximate_w(const int& T, arma::mat sign_narrative, const arma::cube& irf);
static SEXP _bsvarSIGNs_approximate_w_try(SEXP TSEXP, SEXP sign_narrativeSEXP, SEXP irfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const int& >::type T(TSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sign_narrative(sign_narrativeSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type irf(irfSEXP);
    rcpp_result_gen = Rcpp::wrap(approximate_w(T, sign_narrative, irf));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _bsvarSIGNs_approximate_w(SEXP TSEXP, SEXP sign_narrativeSEXP, SEXP irfSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_bsvarSIGNs_approximate_w_try(TSEXP, sign_narrativeSEXP, irfSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// sample_Q
arma::mat sample_Q(const int& lags, const arma::mat& Y, const arma::mat& X, double& aux_w, arma::mat& aux_A, arma::mat& aux_B, arma::mat& aux_hyper, const Rcpp::List& prior, const arma::field<arma::mat>& VB, const arma::cube& sign_irf, const arma::mat& sign_narrative, const arma::mat& sign_B, const int& max_tries, bool& success);
static SEXP _bsvarSIGNs_sample_Q_try(SEXP lagsSEXP, SEXP YSEXP, SEXP XSEXP, SEXP aux_wSEXP, SEXP aux_ASEXP, SEXP aux_BSEXP, SEXP aux_hyperSEXP, SEXP priorSEXP, SEXP VBSEXP, SEXP sign_irfSEXP, SEXP sign_narrativeSEXP, SEXP sign_BSEXP, SEXP max_triesSEXP, SEXP successSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const int& >::type lags(lagsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< double& >::type aux_w(aux_wSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type aux_A(aux_ASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type aux_B(aux_BSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type aux_hyper(aux_hyperSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::mat>& >::type VB(VBSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type sign_irf(sign_irfSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sign_narrative(sign_narrativeSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sign_B(sign_BSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_tries(max_triesSEXP);
    Rcpp::traits::input_parameter< bool& >::type success(successSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_Q(lags, Y, X, aux_w, aux_A, aux_B, aux_hyper, prior, VB, sign_irf, sign_narrative, sign_B, max_tries, success));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _bsvarSIGNs_sample_Q(SEXP lagsSEXP, SEXP YSEXP, SEXP XSEXP, SEXP aux_wSEXP, SEXP aux_ASEXP, SEXP aux_BSEXP, SEXP aux_hyperSEXP, SEXP priorSEXP, SEXP VBSEXP, SEXP sign_irfSEXP, SEXP sign_narrativeSEXP, SEXP sign_BSEXP, SEXP max_triesSEXP, SEXP successSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_bsvarSIGNs_sample_Q_try(lagsSEXP, YSEXP, XSEXP, aux_wSEXP, aux_ASEXP, aux_BSEXP, aux_hyperSEXP, priorSEXP, VBSEXP, sign_irfSEXP, sign_narrativeSEXP, sign_BSEXP, max_triesSEXP, successSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// matnrnd_cpp
arma::mat matnrnd_cpp(const arma::mat& M, const arma::mat& U, const arma::mat& V);
RcppExport SEXP _bsvarSIGNs_matnrnd_cpp(SEXP MSEXP, SEXP USEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(matnrnd_cpp(M, U, V));
    return rcpp_result_gen;
END_RCPP
}
// qr_sign_cpp
arma::mat qr_sign_cpp(const arma::mat& A);
RcppExport SEXP _bsvarSIGNs_qr_sign_cpp(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(qr_sign_cpp(A));
    return rcpp_result_gen;
END_RCPP
}
// rortho_cpp
arma::mat rortho_cpp(const int& N);
RcppExport SEXP _bsvarSIGNs_rortho_cpp(SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(rortho_cpp(N));
    return rcpp_result_gen;
END_RCPP
}
// irf_cpp
arma::field<arma::cube> irf_cpp(arma::cube& posterior_B, arma::cube& posterior_A, const int horizon, const int p);
RcppExport SEXP _bsvarSIGNs_irf_cpp(SEXP posterior_BSEXP, SEXP posterior_ASEXP, SEXP horizonSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type posterior_B(posterior_BSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type posterior_A(posterior_ASEXP);
    Rcpp::traits::input_parameter< const int >::type horizon(horizonSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(irf_cpp(posterior_B, posterior_A, horizon, p));
    return rcpp_result_gen;
END_RCPP
}

// validate (ensure exported C++ functions exist before calling them)
static int _bsvarSIGNs_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("Rcpp::List(*bsvar_sign_cpp)(const int&,const int&,const arma::mat&,const arma::mat&,const arma::field<arma::mat>&,const arma::cube&,const arma::mat&,const arma::mat&,const Rcpp::List&,const Rcpp::List&,const int,const bool,const int&)");
        signatures.insert("Rcpp::List(*bsvar_sign_gibbs_cpp)(const int&,const int&,const arma::mat&,const arma::mat&,const arma::field<arma::mat>&,const arma::cube&,const arma::mat&,const arma::mat&,const Rcpp::List&,const Rcpp::List&,const int,const bool,const int&)");
        signatures.insert("bool(*match_sign)(const arma::mat&,const arma::mat&)");
        signatures.insert("bool(*match_sign_irf)(const arma::mat&,const arma::cube&,const arma::cube&)");
        signatures.insert("arma::mat(*hd1)(const int&,const int&,const int&,const arma::mat&,const arma::cube&)");
        signatures.insert("bool(*match_sign_narrative)(const arma::mat&,const arma::mat&,const arma::cube&)");
        signatures.insert("double(*approximate_w)(const int&,arma::mat,const arma::cube&)");
        signatures.insert("arma::mat(*sample_Q)(const int&,const arma::mat&,const arma::mat&,double&,arma::mat&,arma::mat&,arma::mat&,const Rcpp::List&,const arma::field<arma::mat>&,const arma::cube&,const arma::mat&,const arma::mat&,const int&,bool&)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _bsvarSIGNs_RcppExport_registerCCallable() { 
    R_RegisterCCallable("bsvarSIGNs", "_bsvarSIGNs_bsvar_sign_cpp", (DL_FUNC)_bsvarSIGNs_bsvar_sign_cpp_try);
    R_RegisterCCallable("bsvarSIGNs", "_bsvarSIGNs_bsvar_sign_gibbs_cpp", (DL_FUNC)_bsvarSIGNs_bsvar_sign_gibbs_cpp_try);
    R_RegisterCCallable("bsvarSIGNs", "_bsvarSIGNs_match_sign", (DL_FUNC)_bsvarSIGNs_match_sign_try);
    R_RegisterCCallable("bsvarSIGNs", "_bsvarSIGNs_match_sign_irf", (DL_FUNC)_bsvarSIGNs_match_sign_irf_try);
    R_RegisterCCallable("bsvarSIGNs", "_bsvarSIGNs_hd1", (DL_FUNC)_bsvarSIGNs_hd1_try);
    R_RegisterCCallable("bsvarSIGNs", "_bsvarSIGNs_match_sign_narrative", (DL_FUNC)_bsvarSIGNs_match_sign_narrative_try);
    R_RegisterCCallable("bsvarSIGNs", "_bsvarSIGNs_approximate_w", (DL_FUNC)_bsvarSIGNs_approximate_w_try);
    R_RegisterCCallable("bsvarSIGNs", "_bsvarSIGNs_sample_Q", (DL_FUNC)_bsvarSIGNs_sample_Q_try);
    R_RegisterCCallable("bsvarSIGNs", "_bsvarSIGNs_RcppExport_validate", (DL_FUNC)_bsvarSIGNs_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_bsvarSIGNs_bsvar_sign_cpp", (DL_FUNC) &_bsvarSIGNs_bsvar_sign_cpp, 13},
    {"_bsvarSIGNs_bsvar_sign_gibbs_cpp", (DL_FUNC) &_bsvarSIGNs_bsvar_sign_gibbs_cpp, 13},
    {"_bsvarSIGNs_niw_cpp", (DL_FUNC) &_bsvarSIGNs_niw_cpp, 7},
    {"_bsvarSIGNs_match_sign", (DL_FUNC) &_bsvarSIGNs_match_sign, 2},
    {"_bsvarSIGNs_match_sign_irf", (DL_FUNC) &_bsvarSIGNs_match_sign_irf, 3},
    {"_bsvarSIGNs_hd1", (DL_FUNC) &_bsvarSIGNs_hd1, 5},
    {"_bsvarSIGNs_match_sign_narrative", (DL_FUNC) &_bsvarSIGNs_match_sign_narrative, 3},
    {"_bsvarSIGNs_approximate_w", (DL_FUNC) &_bsvarSIGNs_approximate_w, 3},
    {"_bsvarSIGNs_sample_Q", (DL_FUNC) &_bsvarSIGNs_sample_Q, 14},
    {"_bsvarSIGNs_matnrnd_cpp", (DL_FUNC) &_bsvarSIGNs_matnrnd_cpp, 3},
    {"_bsvarSIGNs_qr_sign_cpp", (DL_FUNC) &_bsvarSIGNs_qr_sign_cpp, 1},
    {"_bsvarSIGNs_rortho_cpp", (DL_FUNC) &_bsvarSIGNs_rortho_cpp, 1},
    {"_bsvarSIGNs_irf_cpp", (DL_FUNC) &_bsvarSIGNs_irf_cpp, 4},
    {"_bsvarSIGNs_RcppExport_registerCCallable", (DL_FUNC) &_bsvarSIGNs_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_bsvarSIGNs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
