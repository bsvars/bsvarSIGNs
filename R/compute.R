

#' @method compute_structural_shocks PosteriorBSVARSIGN
#' 
#' @title Computes posterior draws of structural shocks
#'
#' @description Each of the draws from the posterior estimation of models from
#' packages \pkg{bsvars} or \pkg{bsvarSIGNs} is transformed into
#' a draw from the posterior distribution of the structural shocks. 
#' 
#' @param posterior posterior estimation outcome - an object of class 
#' \code{PosteriorBSVARSIGN} obtained by running the \code{estimate} function.
#' 
#' @return An object of class \code{PosteriorShocks}, that is, an \code{NxTxS} 
#' array with attribute \code{PosteriorShocks} containing \code{S} draws of the 
#' structural shocks.
#'
#' @seealso \code{\link{estimate.BSVARSIGN}}, \code{\link{summary}}, \code{\link{plot}}
#'
#' @author Xiaolei Wang \email{adamwang15@gmail.com} and Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' # upload data
#' data(optimism)
#' 
#' # specify the model and set seed
#' set.seed(123)
#' 
#' # + no effect on productivity (zero restriction)
#' # + positive effect on stock prices (positive sign restriction) 
#' sign_irf       = matrix(c(0, 1, rep(NA, 23)), 5, 5)
#' specification  = specify_bsvarSIGN$new(optimism, sign_irf = sign_irf)
#' 
#' # estimate the model
#' posterior      = estimate(specification, 10)
#' 
#' # compute structural shocks
#' shocks         = compute_structural_shocks(posterior)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' optimism |>
#'   specify_bsvarSIGN$new(sign_irf = sign_irf) |> 
#'   estimate(S = 20) |> 
#'   compute_structural_shocks() -> ss
#' 
#' @export
compute_structural_shocks.PosteriorBSVARSIGN <- function(posterior) {
  
  posterior_B     = posterior$posterior$B
  posterior_A     = posterior$posterior$A
  Y               = posterior$last_draw$data_matrices$Y
  X               = posterior$last_draw$data_matrices$X
  
  ss              = .Call(`_bsvarSIGNs_bsvarSIGNs_structural_shocks`, posterior_B, posterior_A, Y, X)
  class(ss)       = "PosteriorShocks"
  
  return(ss)
} # END compute_structural_shocks.PosteriorBSVARSIGN





#' @method compute_fitted_values PosteriorBSVARSIGN
#' 
#' @title Computes posterior draws from data predictive density
#'
#' @description Each of the draws from the posterior estimation of models from 
#' packages \pkg{bsvars} or \pkg{bsvarSIGNs} is transformed into
#' a draw from the data predictive density. 
#' 
#' @param posterior posterior estimation outcome - an object of class 
#' \code{PosteriorBSVARSIGN} obtained by running the \code{estimate} function.
#' 
#' @return An object of class \code{PosteriorFitted}, that is, an \code{NxTxS} 
#' array with attribute \code{PosteriorFitted} containing \code{S} draws from 
#' the data predictive density.
#'
#' @seealso \code{\link{estimate.BSVARSIGN}}, \code{\link{summary}}, \code{\link{plot}}
#'
#' @author Xiaolei Wang \email{adamwang15@gmail.com} and Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' # upload data
#' data(optimism)
#' 
#' # specify the model and set seed
#' set.seed(123)
#' 
#' # + no effect on productivity (zero restriction)
#' # + positive effect on stock prices (positive sign restriction) 
#' sign_irf       = matrix(c(0, 1, rep(NA, 23)), 5, 5)
#' specification  = specify_bsvarSIGN$new(optimism, sign_irf = sign_irf)
#' 
#' # estimate the model
#' posterior      = estimate(specification, 10)
#' 
#' # compute draws from in-sample predictive density
#' fitted         = compute_fitted_values(posterior)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' optimism |>
#'   specify_bsvarSIGN$new(sign_irf = sign_irf) |> 
#'   estimate(S = 20) |> 
#'   compute_fitted_values() -> fitted
#' 
#' @export
compute_fitted_values.PosteriorBSVARSIGN <- function(posterior) {
  
  posterior_A     = posterior$posterior$A
  posterior_B     = posterior$posterior$B
  
  N               = dim(posterior_A)[1]
  T               = dim(posterior$last_draw$data_matrices$X)[2]
  S               = dim(posterior_A)[3]
  posterior_sigma = array(1, c(N, T, S))
  X               = posterior$last_draw$data_matrices$X
  
  fv              = .Call(`_bsvarSIGNs_bsvarSIGNs_fitted_values`, posterior_A, posterior_B, posterior_sigma, X)
  class(fv)       = "PosteriorFitted"
  
  return(fv)
} # END compute_fitted_values.PosteriorBSVARSIGN





#' @method compute_impulse_responses PosteriorBSVARSIGN
#' 
#' @title Computes posterior draws of impulse responses 
#'
#' @description Each of the draws from the posterior estimation of models from 
#' packages \pkg{bsvars} or \pkg{bsvarSIGNs} is transformed into
#' a draw from the posterior distribution of the impulse responses. 
#' 
#' @param posterior posterior estimation outcome - an object of class 
#' \code{PosteriorBSVARSIGN} obtained by running the \code{estimate} function.
#' @param horizon a positive integer number denoting the forecast horizon for the impulse responses computations.
#' @param standardise a logical value. If \code{TRUE}, the impulse responses are standardised 
#' so that the variables' own shocks at horizon 0 are equal to 1. Otherwise, the parameter estimates 
#' determine this magnitude.
#' 
#' @return An object of class PosteriorIR, that is, an \code{NxNx(horizon+1)xS} array with attribute PosteriorIR 
#' containing \code{S} draws of the impulse responses.
#'
#' @seealso \code{\link{estimate.BSVARSIGN}}, \code{\link{summary}}, \code{\link{plot}}
#'
#' @author Xiaolei Wang \email{adamwang15@gmail.com} and Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Kilian, L., & Lütkepohl, H. (2017). Structural VAR Tools, Chapter 4, In: Structural vector autoregressive analysis. Cambridge University Press.
#' 
#' @examples
#' # upload data
#' data(optimism)
#' 
#' # specify the model and set seed
#' set.seed(123)
#' 
#' # + no effect on productivity (zero restriction)
#' # + positive effect on stock prices (positive sign restriction) 
#' sign_irf       = matrix(c(0, 1, rep(NA, 23)), 5, 5)
#' specification  = specify_bsvarSIGN$new(optimism, sign_irf = sign_irf)
#' 
#' # estimate the model
#' posterior      = estimate(specification, 10)
#' 
#' # compute impulse responses 2 years ahead
#' irf           = compute_impulse_responses(posterior, horizon = 8)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' optimism |>
#'   specify_bsvarSIGN$new(sign_irf = sign_irf) |> 
#'   estimate(S = 10) |> 
#'   compute_impulse_responses(horizon = 8) -> ir
#' 
#' 
#' @export
compute_impulse_responses.PosteriorBSVARSIGN <- function(posterior, horizon, standardise = FALSE) {
  
  if ( 
    any(diag(posterior$last_draw$identification$sign_irf[,,1]) == 0) &
    !is.na(any(diag(posterior$last_draw$identification$sign_irf[,,1]) == 0))
  ) {
    standardise = FALSE
    message("Argument standardise is forcibly set to FALSE due to zero restrictions imposed on the diagonal element(s) of the on-impact impulse response matrix.")
  }
  
  posterior_Theta0  = posterior$posterior$Theta0
  posterior_A       = posterior$posterior$A
  posterior_A       = aperm(posterior_A, c(2, 1, 3))
  N                 = dim(posterior_A)[2]
  p                 = posterior$last_draw$p
  S                 = dim(posterior_A)[3]
  
  qqq               = .Call(`_bsvarSIGNs_bsvarSIGNs_ir`, posterior_A, posterior_Theta0, horizon, p, standardise)
  
  irfs              = array(NA, c(N, N, horizon + 1, S))
  for (s in 1:S) irfs[,,,s] = qqq[s][[1]]
  class(irfs)       = "PosteriorIR"
  
  return(irfs)
} # END compute_impulse_responses.PosteriorBSVARSIGN




#' @method compute_historical_decompositions PosteriorBSVARSIGN
#' 
#' @title Computes posterior draws of historical decompositions
#'
#' @description Each of the draws from the posterior estimation of models from
#' packages \pkg{bsvars} or \pkg{bsvarSIGNs} is transformed into
#' a draw from the posterior distribution of the historical decompositions. 
#' IMPORTANT! The historical decompositions are interpreted correctly for 
#' covariance stationary data. Application to unit-root non-stationary data might
#' result in non-interpretable outcomes.
#' 
#' @param posterior posterior estimation outcome - an object of class 
#' \code{PosteriorBSVARSIGN} obtained by running the \code{estimate} function.
#' @param show_progress a logical value, if \code{TRUE} the estimation progress bar is visible
#' 
#' @return An object of class \code{PosteriorHD}, that is, an \code{NxNxTxS} array 
#' with attribute \code{PosteriorHD} containing \code{S} draws of the historical 
#' decompositions.
#'
#' @seealso \code{\link{estimate.BSVARSIGN}}, \code{\link{summary}}, \code{\link{plot}}
#'
#' @author Xiaolei Wang \email{adamwang15@gmail.com} and Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Kilian, L., & Lütkepohl, H. (2017). Structural VAR Tools, Chapter 4, In: Structural vector autoregressive analysis. Cambridge University Press.
#' 
#' @examples
#' # upload data
#' data(optimism)
#' 
#' # specify the model and set seed
#' set.seed(123)
#' 
#' # + no effect on productivity (zero restriction)
#' # + positive effect on stock prices (positive sign restriction) 
#' sign_irf       = matrix(c(0, 1, rep(NA, 23)), 5, 5)
#' specification  = specify_bsvarSIGN$new(optimism, sign_irf = sign_irf)
#' 
#' # estimate the model
#' posterior      = estimate(specification, 10)
#' 
#' # compute historical decompositions
#' hd            = compute_historical_decompositions(posterior)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' optimism |>
#'   specify_bsvarSIGN$new(sign_irf = sign_irf) |> 
#'   estimate(S = 10) |> 
#'   compute_historical_decompositions() -> hd
#'   
#' @export
compute_historical_decompositions.PosteriorBSVARSIGN <- function(posterior, show_progress = TRUE) {
  
  posterior_Theta0  = posterior$posterior$Theta0
  posterior_B       = posterior$posterior$B
  posterior_A       = posterior$posterior$A
  posterior_At      = aperm(posterior_A, c(2, 1, 3))
  
  Y                 = posterior$last_draw$data_matrices$Y
  X                 = posterior$last_draw$data_matrices$X
  
  N                 = nrow(Y)
  T                 = ncol(Y)
  p                 = posterior$last_draw$p
  S                 = dim(posterior_A)[3]
  
  standardise       = TRUE
  if ( 
    any(diag(posterior$last_draw$identification$sign_irf[,,1]) == 0) &
    !is.na(any(diag(posterior$last_draw$identification$sign_irf[,,1]) == 0))
  ) {
    standardise     = FALSE
  }
  
  ss                = .Call(`_bsvarSIGNs_bsvarSIGNs_structural_shocks`, posterior_B, posterior_A, Y, X)
  ir                = .Call(`_bsvarSIGNs_bsvarSIGNs_ir`, posterior_At, posterior_Theta0, T, p, standardise)
  qqq               = .Call(`_bsvarSIGNs_bsvarSIGNs_hd`, ir, ss, show_progress)
  
  hd                = array(NA, c(N, N, T, S))
  for (s in 1:S) hd[,,,s] = qqq[s][[1]]
  class(hd)         = "PosteriorHD"
  
  return(hd)
} # END compute_historical_decompositions.PosteriorBSVARSIGN




#' @method compute_variance_decompositions PosteriorBSVARSIGN
#' 
#' @title Computes posterior draws of the forecast error variance decomposition
#' 
#' @description Each of the draws from the posterior estimation of the model
#' is transformed into a draw from the posterior distribution of the forecast 
#' error variance decomposition. 
#' 
#' @param posterior posterior estimation outcome - an object of class 
#' \code{PosteriorBSVARSIGN} obtained by running the \code{estimate} function.
#' @param horizon a positive integer number denoting the forecast horizon for 
#' the impulse responses computations.
#' 
#' @return An object of class \code{PosteriorFEVD}, that is, an \code{NxNx(horizon+1)xS} 
#' array with attribute \code{PosteriorFEVD} containing \code{S} draws of the 
#' forecast error variance decomposition.
#'
#' @seealso \code{\link{compute_impulse_responses.PosteriorBSVARSIGN}}, \code{\link{estimate.BSVARSIGN}}, \code{\link{summary}}, \code{\link{plot}}
#'
#' @author Xiaolei Wang \email{adamwang15@gmail.com} and Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Kilian, L., & Lütkepohl, H. (2017). Structural VAR Tools, Chapter 4, In: Structural vector autoregressive analysis. Cambridge University Press.
#' 
#' @examples
#' # upload data
#' data(optimism)
#' 
#' # specify the model and set seed
#' set.seed(123)
#' 
#' # + no effect on productivity (zero restriction)
#' # + positive effect on stock prices (positive sign restriction) 
#' sign_irf       = matrix(c(0, 1, rep(NA, 23)), 5, 5)
#' specification  = specify_bsvarSIGN$new(optimism, sign_irf = sign_irf)
#' 
#' # estimate the model
#' posterior      = estimate(specification, 10)
#' 
#' # compute forecast error variance decomposition 2 years ahead
#' fevd           = compute_variance_decompositions(posterior, horizon = 8)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' optimism |>
#'   specify_bsvarSIGN$new(sign_irf = sign_irf) |> 
#'   estimate(S = 10) |> 
#'   compute_variance_decompositions(horizon = 8) -> fevd
#'   
#' @export
compute_variance_decompositions.PosteriorBSVARSIGN <- function(posterior, horizon) {
  
  posterior_Theta0  = posterior$posterior$Theta0
  posterior_A       = posterior$posterior$A
  posterior_A       = aperm(posterior_A, c(2, 1, 3))
  N                 = dim(posterior_A)[2]
  p                 = posterior$last_draw$p
  S                 = dim(posterior_A)[3]
  
  standardise       = TRUE
  if ( 
    any(diag(posterior$last_draw$identification$sign_irf[,,1]) == 0) &
    !is.na(any(diag(posterior$last_draw$identification$sign_irf[,,1]) == 0))
  ) {
    standardise     = FALSE
  }
  
  posterior_irf     = .Call(`_bsvarSIGNs_bsvarSIGNs_ir`, posterior_A, posterior_Theta0, horizon, p, standardise)
  qqq               = .Call(`_bsvarSIGNs_bsvarSIGNs_fevd`, posterior_irf)
  
  fevd              = array(NA, c(N, N, horizon + 1, S))
  for (s in 1:S) fevd[,,,s] = qqq[s][[1]]
  class(fevd)       = "PosteriorFEVD"
  
  return(fevd)
} # END compute_variance_decompositions.PosteriorBSVARSIGN




#' @method compute_conditional_sd PosteriorBSVARSIGN
#' 
#' @title Computes posterior draws of structural shock conditional standard deviations
#'
#' @description Each of the draws from the posterior estimation of models is 
#' transformed into a draw from the posterior distribution of the structural 
#' shock conditional standard deviations. 
#' 
#' @param posterior posterior estimation outcome - an object of class 
#' \code{PosteriorBSVARSIGN} obtained by running the \code{estimate} function.
#' 
#' @return An object of class \code{PosteriorSigma}, that is, an \code{NxTxS} 
#' array with attribute \code{PosteriorSigma} containing \code{S} draws of the 
#' structural shock conditional standard deviations.
#'
#' @seealso \code{\link{estimate.BSVARSIGN}}
#'
#' @author Xiaolei Wang \email{adamwang15@gmail.com} and Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' # upload data
#' data(optimism)
#' 
#' # specify the model and set seed
#' set.seed(123)
#' 
#' # + no effect on productivity (zero restriction)
#' # + positive effect on stock prices (positive sign restriction) 
#' sign_irf       = matrix(c(0, 1, rep(NA, 23)), 5, 5)
#' specification  = specify_bsvarSIGN$new(optimism, sign_irf = sign_irf)
#' 
#' # estimate the model
#' posterior      = estimate(specification, 10)
#' 
#' # compute structural shocks' conditional standard deviations
#' sigma          = compute_conditional_sd(posterior)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' optimism |>
#'   specify_bsvarSIGN$new(sign_irf = sign_irf) |> 
#'   estimate(S = 10) |> 
#'   compute_conditional_sd() -> csd
#' 
#' @export
compute_conditional_sd.PosteriorBSVARSIGN <- function(posterior) {
  
  Y     = posterior$last_draw$data_matrices$Y
  N     = nrow(Y)
  T     = ncol(Y)
  S     = dim(posterior$posterior$A)[3]
  
  posterior_sigma       = array(1, c(N, T, S))
  message("The model is homoskedastic. Returning an NxTxS matrix of conditional sd all equal to 1.")
  class(posterior_sigma)  = "PosteriorSigma"
  
  return(posterior_sigma)
} # END compute_conditional_sd.PosteriorBSVARSIGN
