// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_bsvarSIGNs_RCPPEXPORTS_H_GEN_
#define RCPP_bsvarSIGNs_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace bsvarSIGNs {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("bsvarSIGNs", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in bsvarSIGNs");
            }
        }
    }

    inline Rcpp::List bsvar_sign_cpp(const int& S, const int& p, const arma::mat& Y, const arma::mat& X, const arma::cube& sign_irf, const arma::mat& sign_narrative, const arma::mat& sign_B, const arma::field<arma::mat>& Z, const Rcpp::List& prior, const bool show_progress, const bool parallel, const int& max_tries) {
        typedef SEXP(*Ptr_bsvar_sign_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_bsvar_sign_cpp p_bsvar_sign_cpp = NULL;
        if (p_bsvar_sign_cpp == NULL) {
            validateSignature("Rcpp::List(*bsvar_sign_cpp)(const int&,const int&,const arma::mat&,const arma::mat&,const arma::cube&,const arma::mat&,const arma::mat&,const arma::field<arma::mat>&,const Rcpp::List&,const bool,const bool,const int&)");
            p_bsvar_sign_cpp = (Ptr_bsvar_sign_cpp)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_bsvar_sign_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bsvar_sign_cpp(Shield<SEXP>(Rcpp::wrap(S)), Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(sign_irf)), Shield<SEXP>(Rcpp::wrap(sign_narrative)), Shield<SEXP>(Rcpp::wrap(sign_B)), Shield<SEXP>(Rcpp::wrap(Z)), Shield<SEXP>(Rcpp::wrap(prior)), Shield<SEXP>(Rcpp::wrap(show_progress)), Shield<SEXP>(Rcpp::wrap(parallel)), Shield<SEXP>(Rcpp::wrap(max_tries)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline arma::cube bsvarSIGNs_structural_shocks(const arma::cube& posterior_B, const arma::cube& posterior_A, const arma::mat& Y, const arma::mat& X) {
        typedef SEXP(*Ptr_bsvarSIGNs_structural_shocks)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_bsvarSIGNs_structural_shocks p_bsvarSIGNs_structural_shocks = NULL;
        if (p_bsvarSIGNs_structural_shocks == NULL) {
            validateSignature("arma::cube(*bsvarSIGNs_structural_shocks)(const arma::cube&,const arma::cube&,const arma::mat&,const arma::mat&)");
            p_bsvarSIGNs_structural_shocks = (Ptr_bsvarSIGNs_structural_shocks)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_bsvarSIGNs_structural_shocks");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bsvarSIGNs_structural_shocks(Shield<SEXP>(Rcpp::wrap(posterior_B)), Shield<SEXP>(Rcpp::wrap(posterior_A)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::cube >(rcpp_result_gen);
    }

    inline arma::cube bsvarSIGNs_fitted_values(arma::cube& posterior_A, arma::cube& posterior_B, arma::cube& posterior_sigma, arma::mat& X) {
        typedef SEXP(*Ptr_bsvarSIGNs_fitted_values)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_bsvarSIGNs_fitted_values p_bsvarSIGNs_fitted_values = NULL;
        if (p_bsvarSIGNs_fitted_values == NULL) {
            validateSignature("arma::cube(*bsvarSIGNs_fitted_values)(arma::cube&,arma::cube&,arma::cube&,arma::mat&)");
            p_bsvarSIGNs_fitted_values = (Ptr_bsvarSIGNs_fitted_values)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_bsvarSIGNs_fitted_values");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bsvarSIGNs_fitted_values(Shield<SEXP>(Rcpp::wrap(posterior_A)), Shield<SEXP>(Rcpp::wrap(posterior_B)), Shield<SEXP>(Rcpp::wrap(posterior_sigma)), Shield<SEXP>(Rcpp::wrap(X)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::cube >(rcpp_result_gen);
    }

    inline arma::cube ir1_cpp(const arma::mat& B, const arma::mat& Theta0, int horizon, const int& p) {
        typedef SEXP(*Ptr_ir1_cpp)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_ir1_cpp p_ir1_cpp = NULL;
        if (p_ir1_cpp == NULL) {
            validateSignature("arma::cube(*ir1_cpp)(const arma::mat&,const arma::mat&,int,const int&)");
            p_ir1_cpp = (Ptr_ir1_cpp)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_ir1_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_ir1_cpp(Shield<SEXP>(Rcpp::wrap(B)), Shield<SEXP>(Rcpp::wrap(Theta0)), Shield<SEXP>(Rcpp::wrap(horizon)), Shield<SEXP>(Rcpp::wrap(p)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::cube >(rcpp_result_gen);
    }

    inline arma::field<arma::cube> bsvarSIGNs_ir(arma::cube& posterior_B, arma::cube& posterior_Theta0, const int horizon, const int p, const bool standardise = false) {
        typedef SEXP(*Ptr_bsvarSIGNs_ir)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_bsvarSIGNs_ir p_bsvarSIGNs_ir = NULL;
        if (p_bsvarSIGNs_ir == NULL) {
            validateSignature("arma::field<arma::cube>(*bsvarSIGNs_ir)(arma::cube&,arma::cube&,const int,const int,const bool)");
            p_bsvarSIGNs_ir = (Ptr_bsvarSIGNs_ir)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_bsvarSIGNs_ir");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bsvarSIGNs_ir(Shield<SEXP>(Rcpp::wrap(posterior_B)), Shield<SEXP>(Rcpp::wrap(posterior_Theta0)), Shield<SEXP>(Rcpp::wrap(horizon)), Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(standardise)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::field<arma::cube> >(rcpp_result_gen);
    }

    inline arma::field<arma::cube> bsvarSIGNs_hd(arma::field<arma::cube>& posterior_irf_T, arma::cube& structural_shocks, const bool show_progress = true) {
        typedef SEXP(*Ptr_bsvarSIGNs_hd)(SEXP,SEXP,SEXP);
        static Ptr_bsvarSIGNs_hd p_bsvarSIGNs_hd = NULL;
        if (p_bsvarSIGNs_hd == NULL) {
            validateSignature("arma::field<arma::cube>(*bsvarSIGNs_hd)(arma::field<arma::cube>&,arma::cube&,const bool)");
            p_bsvarSIGNs_hd = (Ptr_bsvarSIGNs_hd)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_bsvarSIGNs_hd");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bsvarSIGNs_hd(Shield<SEXP>(Rcpp::wrap(posterior_irf_T)), Shield<SEXP>(Rcpp::wrap(structural_shocks)), Shield<SEXP>(Rcpp::wrap(show_progress)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::field<arma::cube> >(rcpp_result_gen);
    }

    inline arma::mat hd1_cpp(const int& var_i, const int& t, const int& h, const arma::mat& Epsilon, const arma::cube& irf) {
        typedef SEXP(*Ptr_hd1_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_hd1_cpp p_hd1_cpp = NULL;
        if (p_hd1_cpp == NULL) {
            validateSignature("arma::mat(*hd1_cpp)(const int&,const int&,const int&,const arma::mat&,const arma::cube&)");
            p_hd1_cpp = (Ptr_hd1_cpp)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_hd1_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_hd1_cpp(Shield<SEXP>(Rcpp::wrap(var_i)), Shield<SEXP>(Rcpp::wrap(t)), Shield<SEXP>(Rcpp::wrap(h)), Shield<SEXP>(Rcpp::wrap(Epsilon)), Shield<SEXP>(Rcpp::wrap(irf)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::field<arma::cube> bsvarSIGNs_fevd(arma::field<arma::cube>& posterior_irf) {
        typedef SEXP(*Ptr_bsvarSIGNs_fevd)(SEXP);
        static Ptr_bsvarSIGNs_fevd p_bsvarSIGNs_fevd = NULL;
        if (p_bsvarSIGNs_fevd == NULL) {
            validateSignature("arma::field<arma::cube>(*bsvarSIGNs_fevd)(arma::field<arma::cube>&)");
            p_bsvarSIGNs_fevd = (Ptr_bsvarSIGNs_fevd)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_bsvarSIGNs_fevd");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bsvarSIGNs_fevd(Shield<SEXP>(Rcpp::wrap(posterior_irf)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::field<arma::cube> >(rcpp_result_gen);
    }

    inline arma::cube forecast_bsvarSIGNs(arma::cube& posterior_Sigma, arma::cube& posterior_A, arma::vec& X_T, arma::mat& exogenous_forecast, arma::mat& cond_forecast, const int& horizon) {
        typedef SEXP(*Ptr_forecast_bsvarSIGNs)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_forecast_bsvarSIGNs p_forecast_bsvarSIGNs = NULL;
        if (p_forecast_bsvarSIGNs == NULL) {
            validateSignature("arma::cube(*forecast_bsvarSIGNs)(arma::cube&,arma::cube&,arma::vec&,arma::mat&,arma::mat&,const int&)");
            p_forecast_bsvarSIGNs = (Ptr_forecast_bsvarSIGNs)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_forecast_bsvarSIGNs");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_forecast_bsvarSIGNs(Shield<SEXP>(Rcpp::wrap(posterior_Sigma)), Shield<SEXP>(Rcpp::wrap(posterior_A)), Shield<SEXP>(Rcpp::wrap(X_T)), Shield<SEXP>(Rcpp::wrap(exogenous_forecast)), Shield<SEXP>(Rcpp::wrap(cond_forecast)), Shield<SEXP>(Rcpp::wrap(horizon)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::cube >(rcpp_result_gen);
    }

    inline bool match_sign_narrative(const arma::mat& Epsilon, const arma::mat& sign_narrative, const arma::cube& irf) {
        typedef SEXP(*Ptr_match_sign_narrative)(SEXP,SEXP,SEXP);
        static Ptr_match_sign_narrative p_match_sign_narrative = NULL;
        if (p_match_sign_narrative == NULL) {
            validateSignature("bool(*match_sign_narrative)(const arma::mat&,const arma::mat&,const arma::cube&)");
            p_match_sign_narrative = (Ptr_match_sign_narrative)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_match_sign_narrative");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_match_sign_narrative(Shield<SEXP>(Rcpp::wrap(Epsilon)), Shield<SEXP>(Rcpp::wrap(sign_narrative)), Shield<SEXP>(Rcpp::wrap(irf)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<bool >(rcpp_result_gen);
    }

    inline double weight_narrative(const int& T, arma::mat sign_narrative, const arma::cube& irf) {
        typedef SEXP(*Ptr_weight_narrative)(SEXP,SEXP,SEXP);
        static Ptr_weight_narrative p_weight_narrative = NULL;
        if (p_weight_narrative == NULL) {
            validateSignature("double(*weight_narrative)(const int&,arma::mat,const arma::cube&)");
            p_weight_narrative = (Ptr_weight_narrative)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_weight_narrative");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_weight_narrative(Shield<SEXP>(Rcpp::wrap(T)), Shield<SEXP>(Rcpp::wrap(sign_narrative)), Shield<SEXP>(Rcpp::wrap(irf)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline arma::field<arma::mat> ZIRF(const arma::field<arma::mat>& Z, const arma::mat& irf_0) {
        typedef SEXP(*Ptr_ZIRF)(SEXP,SEXP);
        static Ptr_ZIRF p_ZIRF = NULL;
        if (p_ZIRF == NULL) {
            validateSignature("arma::field<arma::mat>(*ZIRF)(const arma::field<arma::mat>&,const arma::mat&)");
            p_ZIRF = (Ptr_ZIRF)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_ZIRF");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_ZIRF(Shield<SEXP>(Rcpp::wrap(Z)), Shield<SEXP>(Rcpp::wrap(irf_0)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::field<arma::mat> >(rcpp_result_gen);
    }

    inline arma::colvec zero_restrictions(const arma::field<arma::mat>& Z, const arma::colvec vec_structural) {
        typedef SEXP(*Ptr_zero_restrictions)(SEXP,SEXP);
        static Ptr_zero_restrictions p_zero_restrictions = NULL;
        if (p_zero_restrictions == NULL) {
            validateSignature("arma::colvec(*zero_restrictions)(const arma::field<arma::mat>&,const arma::colvec)");
            p_zero_restrictions = (Ptr_zero_restrictions)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_zero_restrictions");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_zero_restrictions(Shield<SEXP>(Rcpp::wrap(Z)), Shield<SEXP>(Rcpp::wrap(vec_structural)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::colvec >(rcpp_result_gen);
    }

    inline arma::colvec g_fh(const arma::field<arma::mat>& Z, const arma::mat& A0, const arma::mat& Aplus) {
        typedef SEXP(*Ptr_g_fh)(SEXP,SEXP,SEXP);
        static Ptr_g_fh p_g_fh = NULL;
        if (p_g_fh == NULL) {
            validateSignature("arma::colvec(*g_fh)(const arma::field<arma::mat>&,const arma::mat&,const arma::mat&)");
            p_g_fh = (Ptr_g_fh)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_g_fh");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_g_fh(Shield<SEXP>(Rcpp::wrap(Z)), Shield<SEXP>(Rcpp::wrap(A0)), Shield<SEXP>(Rcpp::wrap(Aplus)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::colvec >(rcpp_result_gen);
    }

    inline arma::colvec g_fh_vec(const arma::field<arma::mat>& Z, const arma::colvec vec_structural) {
        typedef SEXP(*Ptr_g_fh_vec)(SEXP,SEXP);
        static Ptr_g_fh_vec p_g_fh_vec = NULL;
        if (p_g_fh_vec == NULL) {
            validateSignature("arma::colvec(*g_fh_vec)(const arma::field<arma::mat>&,const arma::colvec)");
            p_g_fh_vec = (Ptr_g_fh_vec)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_g_fh_vec");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_g_fh_vec(Shield<SEXP>(Rcpp::wrap(Z)), Shield<SEXP>(Rcpp::wrap(vec_structural)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::colvec >(rcpp_result_gen);
    }

    inline double log_volume_element(const arma::field<arma::mat>& Z, const arma::mat& A0, const arma::mat& Aplus) {
        typedef SEXP(*Ptr_log_volume_element)(SEXP,SEXP,SEXP);
        static Ptr_log_volume_element p_log_volume_element = NULL;
        if (p_log_volume_element == NULL) {
            validateSignature("double(*log_volume_element)(const arma::field<arma::mat>&,const arma::mat&,const arma::mat&)");
            p_log_volume_element = (Ptr_log_volume_element)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_log_volume_element");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_log_volume_element(Shield<SEXP>(Rcpp::wrap(Z)), Shield<SEXP>(Rcpp::wrap(A0)), Shield<SEXP>(Rcpp::wrap(Aplus)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double weight_zero(const arma::field<arma::mat>& Z, const arma::mat& B, const arma::mat& h_inv, const arma::mat& Q) {
        typedef SEXP(*Ptr_weight_zero)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_weight_zero p_weight_zero = NULL;
        if (p_weight_zero == NULL) {
            validateSignature("double(*weight_zero)(const arma::field<arma::mat>&,const arma::mat&,const arma::mat&,const arma::mat&)");
            p_weight_zero = (Ptr_weight_zero)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_weight_zero");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_weight_zero(Shield<SEXP>(Rcpp::wrap(Z)), Shield<SEXP>(Rcpp::wrap(B)), Shield<SEXP>(Rcpp::wrap(h_inv)), Shield<SEXP>(Rcpp::wrap(Q)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline arma::mat rzeroQ(const arma::field<arma::mat>& Z, const arma::mat& irf_0) {
        typedef SEXP(*Ptr_rzeroQ)(SEXP,SEXP);
        static Ptr_rzeroQ p_rzeroQ = NULL;
        if (p_rzeroQ == NULL) {
            validateSignature("arma::mat(*rzeroQ)(const arma::field<arma::mat>&,const arma::mat&)");
            p_rzeroQ = (Ptr_rzeroQ)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_rzeroQ");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rzeroQ(Shield<SEXP>(Rcpp::wrap(Z)), Shield<SEXP>(Rcpp::wrap(irf_0)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline bool match_sign_irf(const arma::mat& Q, const arma::cube& sign_irf, const arma::cube& irf) {
        typedef SEXP(*Ptr_match_sign_irf)(SEXP,SEXP,SEXP);
        static Ptr_match_sign_irf p_match_sign_irf = NULL;
        if (p_match_sign_irf == NULL) {
            validateSignature("bool(*match_sign_irf)(const arma::mat&,const arma::cube&,const arma::cube&)");
            p_match_sign_irf = (Ptr_match_sign_irf)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_match_sign_irf");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_match_sign_irf(Shield<SEXP>(Rcpp::wrap(Q)), Shield<SEXP>(Rcpp::wrap(sign_irf)), Shield<SEXP>(Rcpp::wrap(irf)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<bool >(rcpp_result_gen);
    }

    inline arma::field<arma::mat> sample_Q(const int& p, const arma::mat& Y, const arma::mat& X, arma::mat& B, arma::mat& h_invp, arma::mat& chol_Sigma, const Rcpp::List& prior, const arma::cube& sign_irf, const arma::mat& sign_narrative, const arma::mat& sign_B, const arma::field<arma::mat>& Z, const int& max_tries) {
        typedef SEXP(*Ptr_sample_Q)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_sample_Q p_sample_Q = NULL;
        if (p_sample_Q == NULL) {
            validateSignature("arma::field<arma::mat>(*sample_Q)(const int&,const arma::mat&,const arma::mat&,arma::mat&,arma::mat&,arma::mat&,const Rcpp::List&,const arma::cube&,const arma::mat&,const arma::mat&,const arma::field<arma::mat>&,const int&)");
            p_sample_Q = (Ptr_sample_Q)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_sample_Q");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sample_Q(Shield<SEXP>(Rcpp::wrap(p)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(B)), Shield<SEXP>(Rcpp::wrap(h_invp)), Shield<SEXP>(Rcpp::wrap(chol_Sigma)), Shield<SEXP>(Rcpp::wrap(prior)), Shield<SEXP>(Rcpp::wrap(sign_irf)), Shield<SEXP>(Rcpp::wrap(sign_narrative)), Shield<SEXP>(Rcpp::wrap(sign_B)), Shield<SEXP>(Rcpp::wrap(Z)), Shield<SEXP>(Rcpp::wrap(max_tries)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::field<arma::mat> >(rcpp_result_gen);
    }

    inline arma::mat qr_sign_cpp(const arma::mat& A) {
        typedef SEXP(*Ptr_qr_sign_cpp)(SEXP);
        static Ptr_qr_sign_cpp p_qr_sign_cpp = NULL;
        if (p_qr_sign_cpp == NULL) {
            validateSignature("arma::mat(*qr_sign_cpp)(const arma::mat&)");
            p_qr_sign_cpp = (Ptr_qr_sign_cpp)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_qr_sign_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_qr_sign_cpp(Shield<SEXP>(Rcpp::wrap(A)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::mat rortho_cpp(const int& N) {
        typedef SEXP(*Ptr_rortho_cpp)(SEXP);
        static Ptr_rortho_cpp p_rortho_cpp = NULL;
        if (p_rortho_cpp == NULL) {
            validateSignature("arma::mat(*rortho_cpp)(const int&)");
            p_rortho_cpp = (Ptr_rortho_cpp)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_rortho_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rortho_cpp(Shield<SEXP>(Rcpp::wrap(N)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline bool match_sign(const arma::mat& A, const arma::mat& sign) {
        typedef SEXP(*Ptr_match_sign)(SEXP,SEXP);
        static Ptr_match_sign p_match_sign = NULL;
        if (p_match_sign == NULL) {
            validateSignature("bool(*match_sign)(const arma::mat&,const arma::mat&)");
            p_match_sign = (Ptr_match_sign)R_GetCCallable("bsvarSIGNs", "_bsvarSIGNs_match_sign");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_match_sign(Shield<SEXP>(Rcpp::wrap(A)), Shield<SEXP>(Rcpp::wrap(sign)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<bool >(rcpp_result_gen);
    }

}

#endif // RCPP_bsvarSIGNs_RCPPEXPORTS_H_GEN_
