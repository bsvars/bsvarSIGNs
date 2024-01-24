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
Rcpp::List bsvar_sign_cpp(const int& S, const arma::mat& Y, const arma::mat& X, const arma::field<arma::mat>& VB, const Rcpp::List& prior, const Rcpp::List& starting_values, const int thin, const bool show_progress);
static SEXP _bsvarSIGNs_bsvar_sign_cpp_try(SEXP SSEXP, SEXP YSEXP, SEXP XSEXP, SEXP VBSEXP, SEXP priorSEXP, SEXP starting_valuesSEXP, SEXP thinSEXP, SEXP show_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const int& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::mat>& >::type VB(VBSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type starting_values(starting_valuesSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const bool >::type show_progress(show_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(bsvar_sign_cpp(S, Y, X, VB, prior, starting_values, thin, show_progress));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _bsvarSIGNs_bsvar_sign_cpp(SEXP SSEXP, SEXP YSEXP, SEXP XSEXP, SEXP VBSEXP, SEXP priorSEXP, SEXP starting_valuesSEXP, SEXP thinSEXP, SEXP show_progressSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_bsvarSIGNs_bsvar_sign_cpp_try(SSEXP, YSEXP, XSEXP, VBSEXP, priorSEXP, starting_valuesSEXP, thinSEXP, show_progressSEXP));
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

// validate (ensure exported C++ functions exist before calling them)
static int _bsvarSIGNs_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("Rcpp::List(*bsvar_sign_cpp)(const int&,const arma::mat&,const arma::mat&,const arma::field<arma::mat>&,const Rcpp::List&,const Rcpp::List&,const int,const bool)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _bsvarSIGNs_RcppExport_registerCCallable() { 
    R_RegisterCCallable("bsvarSIGNs", "_bsvarSIGNs_bsvar_sign_cpp", (DL_FUNC)_bsvarSIGNs_bsvar_sign_cpp_try);
    R_RegisterCCallable("bsvarSIGNs", "_bsvarSIGNs_RcppExport_validate", (DL_FUNC)_bsvarSIGNs_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_bsvarSIGNs_bsvar_sign_cpp", (DL_FUNC) &_bsvarSIGNs_bsvar_sign_cpp, 8},
    {"_bsvarSIGNs_matnrnd_cpp", (DL_FUNC) &_bsvarSIGNs_matnrnd_cpp, 3},
    {"_bsvarSIGNs_RcppExport_registerCCallable", (DL_FUNC) &_bsvarSIGNs_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_bsvarSIGNs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
