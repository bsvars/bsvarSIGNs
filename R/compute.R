

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
#' @seealso \code{\link{estimate}}, \code{\link{summary}}, \code{\link{plot}}
#'
#' @author Xiaolei Wang \email{adamwang15@gmail.com} and Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' # upload data
#' data(oil)
#' 
#' # specify the model and set seed
#' set.seed(123)
#' sign_irf       = array(matrix(c(-1, -1, 1, rep(0, 6)), nrow = 3), dim = c(3, 3, 1))
#' specification  = specify_bsvarSIGN$new(oil, sign_irf = sign_irf)
#' 
#' # run the burn-in
#' posterior      = estimate(specification, 10)
#' 
#' # compute structural shocks
#' shocks         = compute_structural_shocks(posterior)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' oil |>
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
  
  ss              = .Call(`_bsvarSIGNs_structural_shocks`, posterior_B, posterior_A, Y, X)
  class(ss)       = "PosteriorShocks"
  
  return(ss)
} # END compute_structural_shocks.PosteriorBSVARSIGN

