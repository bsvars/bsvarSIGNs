#  #####################################################################################
#  R package bsvarSIGNs by Xiaolei Wang Tomasz Woźniak Copyright (C) 2024
#
#  This file is part of the R package bsvarSIGNs: Bayesian Estimation of Structural 
#  Vector Autoregressions Identified by Sign, Zero, and Narrative Restrictions
#
#  The R package bsvarSIGNs is free software: you can redistribute it
#  and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation, either version 3 or
#  any later version of the License.
#
#  The R package bsvarSIGNs is distributed in the hope that it will be
#  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with the R package bsvarSIGNs. If that is not the case, please
#  refer to <http://www.gnu.org/licenses/>.
#  #####################################################################################
#
#' @title Bayesian Estimation of Structural Vector Autoregressions Identified 
#' by Sign, Zero, and Narrative Restrictions
#'
#' @description Implements state-of-the-art algorithms for the Bayesian analysis 
#' of Structural Vector Autoregressions identified by sign, zero, and narrative 
#' restrictions. The core model is based on a flexible Vector Autoregression with 
#' estimated hyper-parameters of the Minnesota prior and the dummy observation priors
#' as in Giannone, Lenza, Primiceri (2015) <doi:10.1162/REST_a_00483>. The sign 
#' restrictions are implemented employing the methods proposed by 
#' Rubio-Ramírez, Waggoner & Zha (2010) <doi:10.1111/j.1467-937X.2009.00578.x>, 
#' while identification through sign and zero restrictions follows the approach 
#' developed by Arias, Rubio-Ramírez, & Waggoner (2018) <doi:10.3982/ECTA14468>. 
#' Furthermore, our tool provides algorithms for identification via sign and 
#' narrative restrictions, in line with the methods introduced by 
#' Antolín-Díaz and Rubio-Ramírez (2018) <doi:10.1257/aer.20161852>. Users can 
#' also estimate a model with sign, zero, and narrative restrictions imposed at 
#' once. The package facilitates predictive and structural analyses using 
#' impulse responses, forecast error variance and historical decompositions, 
#' forecasting and conditional forecasting, as well as analyses of structural 
#' shocks and fitted values. All this is complemented by colourful plots, 
#' user-friendly summary functions, and comprehensive documentation. The 
#' 'bsvarSIGNs' package is aligned regarding objects, workflows, and code 
#' structure with the R package 'bsvars' by 
#' Woźniak (2024) <doi:10.32614/CRAN.package.bsvars>, and they constitute an 
#' integrated toolset.
#' 
#' @details
#' 
#' \strong{Models.} All the SVAR models in this package are specified by two 
#' equations, including the reduced form equation:
#' \deqn{y_t = Ax_t + \epsilon_t}
#' where \eqn{y_t} is an \code{N}-vector of dependent variables, 
#' \eqn{x_t} is a \code{K}-vector of explanatory variables, 
#' \eqn{\epsilon_t} is an \code{N}-vector of reduced form error terms, 
#' and \eqn{A} is an \code{NxK} matrix of autoregressive slope coefficients and 
#' parameters on deterministic terms in \eqn{x_t}.
#' 
#' The structural equation is given by:
#' \deqn{B\epsilon_t = u_t}
#' where \eqn{u_t} is an \code{N}-vector of structural shocks, and
#' \eqn{B} is an \code{NxN} matrix of contemporaneous relationships.
#' 
#' Finally, all of the models share the following assumptions regarding the 
#' structural shocks \code{u_t}, namely, joint conditional normality given the 
#' past observations collected in matrix \code{x_t}, and temporal and 
#' contemporaneous independence. The latter implies zero correlations and 
#' autocorrelations.
#' 
#' \strong{Identification.} The identification of the SVAR model is achieved by 
#' imposing:
#' \itemize{
#'   \item sign restrictions on the structural matrix \eqn{B},
#'   \item sign and zero restrictions on the zero-horizon impulse responses \eqn{\Theta_0 = B^{-1}},
#'   \item sign restrictions on the impulse responses at other horizons \eqn{\Theta_i} for \eqn{i = 1, 2, \ldots},
#'   \item sign restrictions on selected structural shocks \eqn{u_t},
#'   \item two types of sign restrictions on the historical decompositions.
#' }
#' These restrictions determine the sampling algorithms of the structural matrix 
#' \eqn{B} defined as
#' \deqn{B = Q'L}
#' where \eqn{Q} is an \code{NxN} orthogonal matrix and \eqn{L} is a lower-triangular 
#' matrix \eqn{L = chol(\Sigma)^{-1}}, and \eqn{\Sigma} is the \code{NxN} 
#' conditional covariance matrix of the reduced-form error term \eqn{\epsilon_t}.
#' Consult the original papers by Rubio-Ramírez, Waggoner & Zha (2010), 
#' Arias, Rubio-Ramírez, & Waggoner (2018) and Antolín-Díaz and Rubio-Ramírez (2018)
#' for more details.
#' 
#' \strong{Prior distributions.} All the models feature a hierarchical Minnesota 
#' prior following the specification proposed by Giannone, Lenza, Primiceri (2015)
#' and featuring:
#' \itemize{
#'   \item appropriate handling of unit-root non-stationary variables through 
#'   the prior mean of the autoregressive coefficients \eqn{A},
#'   \item normal prior shrinkage exhibiting exponential decay in the lag order 
#'   of the autoregressive matrices,
#'   \item sum-of-coefficients and dummy-initial-observation prior,
#'   \item estimated shrinkage hyper-parameters,
#'   \item inverse-Wishart prior for the reduced-form covariance matrix \eqn{\Sigma},
#'   \item estimated diagonal elements of the inverse-Wishart prior scale matrix.
#' }
#' 
#' @name bsvarSIGNs-package
#' @aliases bsvarSIGNs-package bsvarSIGNs
#' @useDynLib bsvarSIGNs, .registration = TRUE
#' @importFrom R6 R6Class
#' @importFrom Rcpp sourceCpp
#' @importFrom bsvars estimate forecast compute_impulse_responses compute_fitted_values compute_historical_decompositions compute_structural_shocks compute_variance_decompositions
#' @import bsvars
#' @import RcppArmadillo
#' @import RcppProgress
#' 
#' @note This package is currently in active development. Your comments,
#' suggestions and requests are warmly welcome!
#' 
#' @author Xiaolei Wang \email{adamwang15@gmail.com} & Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references
#'  Antolín-Díaz & Rubio-Ramírez (2018) Narrative Sign Restrictions for SVARs, American Economic Review, 108(10), 2802-29, <doi:10.1257/aer.20161852>.
#'  
#'  Arias, Rubio-Ramírez, & Waggoner (2018), Inference Based on Structural Vector Autoregressions Identified With Sign and Zero Restrictions: Theory and Applications, Econometrica, 86(2), 685-720, <doi:10.3982/ECTA14468>.
#'  
#'  Giannone, Lenza, Primiceri (2015) Prior Selection for Vector Autoregressions, Review of Economics and Statistics, 97(2), 436-451 <doi:10.1162/REST_a_00483>.
#'  
#'  Rubio-Ramírez, Waggoner & Zha (2010) Structural Vector Autoregressions: Theory of Identification and Algorithms for Inference, The Review of Economic Studies, 77(2), 665-696, <doi:10.1111/j.1467-937X.2009.00578.x>.
#'  
#'  Woźniak (2024) bsvars: Bayesian Estimation of Structural Vector Autoregressive Models. R package version 3.1, <doi:10.32614/CRAN.package.bsvars>.
#'  
#' @examples
#' # investigate the effects of the optimism shock
#' data(optimism)
#'
#' # specify identifying restrictions:
#' # + no effect on productivity (zero restriction)
#' # + positive effect on stock prices (positive sign restriction) 
#' sign_irf       = matrix(c(0, 1, rep(NA, 23)), 5, 5)
#' 
#' # specify the model and set seed
#' set.seed(123)
#' specification  = specify_bsvarSIGN$new(optimism, sign_irf = sign_irf)
#'                                        
#' # estimate the model
#' posterior      = estimate(specification, S = 10)
#' 
#' # compute and plot impulse responses
#' irf            = compute_impulse_responses(posterior, horizon = 40)
#' plot(irf, probability = 0.68)
#' 
NULL
