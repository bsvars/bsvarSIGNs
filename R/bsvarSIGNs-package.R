#' @title Bayesian Estimation of Structural Vector Autoregressive Models Identified by Sign and Narrative Restrictions
#'
#' @description Implements efficient algorithms for the Bayesian estimation of
#' Stuructural Vector Autoregressive models identified by sign and narrative
#' restrictions following Rubio-Ramírez, Waggoner & Zha (2010)
#' <doi:10.1111/j.1467-937X.2009.00578.x> and Antolín-Díaz & Rubio-Ramírez (2018)
#' <doi:10.1257/aer.20161852>.
# #' @details
#' @name bsvarSIGNs-package
#' @aliases bsvarSIGNs-package bsvarSIGNs
#' "_PACKAGE"
#' @useDynLib bsvarSIGNs, .registration = TRUE
#' @importFrom R6 R6Class
#' @importFrom Rcpp sourceCpp
#' @import bsvars
#' @import RcppArmadillo
#' @import RcppProgress
#' @note This package is currently in active development. We give no
#' warranty that anything here works.
#' @author Xiaolei Wang \email{adamwang15@gmail.com} Tomasz Woźniak \email{wozniak.tom@pm.me}
#' @references
#'  Rubio-Ramírez, Waggoner & Zha (2010) Structural Vector Autoregressions: Theory of Identification and Algorithms for Inference, The Review of Economic Studies, 77(2), 665-696, <doi:10.1111/j.1467-937X.2009.00578.x>.
#'  Antolín-Díaz & Rubio-Ramírez (2018) Narrative Sign Restrictions for SVARs, American Economic Review, 108(10), 2802-29, <doi:10.1257/aer.20161852>.
NULL
