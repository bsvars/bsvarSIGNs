#' @title Bayesian estimation of the model's hyper-parameters
#'
#' @description Estimates the hyper-parameters of the model including scalar
#' parameters \eqn{\mu}, \eqn{\delta}, and \eqn{\lambda}, as well as the 
#' N-vector \eqn{\psi} following the approach by Giannone, Lenza, & Primiceri (2015).
#' The estimation is performed marginally on other parameters of the model, namely,
#' the matrix parameters \eqn{A} and \eqn{\Sigma} using the Adaptive scaling 
#' within the Adaptive Random Walk Metropolis–Hastings by Andrieu, Moulines (2006) 
#' and Atchadé, Fort (2010). Our implementation of this algorithm closely follows
#' the exposition in Kalli, Griffin (2018). See more in \bold{Details}.
#' 
#' @details Let the \code{3+N} vector \eqn{\theta} collect the hyper-parameters  
#' \eqn{\mu}, \eqn{\delta}, \eqn{\lambda}, and \eqn{\psi}.
#' The hyper parameters are sampled from the marginal posterior 
#' distribution proportional to the marginalised likelihood \eqn{p(y|\theta)} 
#' and the prior distribution \eqn{p(\theta)}.
#' \deqn{p(\theta|y) \propto p(y|\theta)p(\theta)}
#' The likelihood \eqn{p(y|\theta)} is obtained by integrating out the model 
#' parameters \eqn{A} and \eqn{\Sigma}:
#' \deqn{p(y|\theta) = \int p(y|A,\Sigma)p(A,\Sigma|\theta)d(A,\Sigma)}
#' The analytical expression for this quantity is know and given by Giannone, 
#' Lenza, & Primiceri (2015). 
#' 
#' Let \eqn{\vartheta} denote a vector of logarithms of the corresponding elements
#' of vector \eqn{\theta}. The Adaptive Random Walk Metropolis–Hastings algorithm
#' provides draws sampled from \eqn{p(\theta|y)} using the candidate generating
#' density that in the \code{s}th iteration provides the candidate draw 
#' \eqn{\vartheta^*} sampled from the \code{3+N} variate normal distribution:
#' \deqn{\vartheta^* \sim\N_{3+N}( \vartheta^{(s-1)}, \sigma_{s}^2\Omega_{s} )}
#' where \eqn{\vartheta^{(s-1)}} is the current state of the Markov chain, 
#' \eqn{\Omega_{s}} is the covariance of the first \code{s-1} draws of \eqn{\vartheta^*},
#' and \eqn{\sigma_{s}^2} is the adaptive scaling parameter determined by the 
#' dynamic equation:
#' \deqn{\sigma_{s+1} = \sigma_{s} + s^{-0.6}(\alpha_s - 0.234)}
#' 
#' @return An object of class \code{BSVARSIGN} containing the Bayesian estimation 
#' output for the hyper-parameters.
#' 
#' @seealso \code{\link{specify_bsvarSIGN}}, \code{\link{estimate.BSVARSIGN}}
#'
#' @author Xiaolei Wang \email{adamwang15@gmail.com} and Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' 
#' Andrieu, Moulines (2006) On the Ergodicity Properties of Some Adaptive MCMC 
#' Algorithms." Annals of Applied Probability, 16(3), 1462-1505, 
#' <doi:10.1214/105051606000000286>.
#' 
#' Atchadé, Fort (2010) Limit Theorems for Some Adaptive MCMC Algorithms with 
#' Subgeometric Kernels, Bernoulli, 16(1) 116-154, <doi:10.3150/09-BEJ199>.
#' 
#' Giannone, Lenza, Primiceri (2015) Prior Selection for Vector Autoregressions, 
#' Review of Economics and Statistics, 97(2), 436-451 <doi:10.1162/REST_a_00483>.
#' 
#' Kalli, Griffin (2018) Bayesian Nonparametric Vector Autoregressive Models,
#' Journal of Econometrics, 203(2), 267-282, <doi:10.1016/j.jeconom.2017.11.009>.
#' 
#' @examples
#' # simple workflow
#' ############################################################
#' # investigate the effects of the optimism shock
#' data(optimism)
#'
#' # specify identifying restrictions:
#' sign_irf      = matrix(c(0, 1, rep(NA, 23)), 5, 5)
#' 
#' # specify the model and set seed
#' set.seed(123)
#' specification = specify_bsvarSIGN$new(
#'                   optimism * 100,
#'                   p        = 12,
#'                   sign_irf = sign_irf
#'                  )
#'                  
#' # estimate hyper-parameters
#' specification = estimate_hyper(specification, S = 10, S_burn = 5)
#' 
#' # estimate the model
#' posterior     = estimate(specification)
#' 
#' # workflow with a pipe
#' sign_irf       = matrix(c(0, 1, rep(NA, 23)), 5, 5)
#' 100 * optimism |> 
#'   specify_bsvarSIGN$new(p = 12, sign_irf = sign_irf) |> 
#'   estimate_hyper(S = 10, S_burn = 5) |> 
#'   estimate() -> posterior
#'   
#' @export
estimate_hyper <- function(specification, S, S_burn, thin = 1, show_progress = TRUE) {
  
  # check the inputs
  stopifnot("Argument S must be a positive integer number." = S > 1 & S %% 1 == 0)
  stopifnot("Argument S_burn must be a positive integer number." = S_burn > 1 & S_burn %% 1 == 0)
  stopifnot("Argument thin must be a positive integer number." = thin > 0 & thin %% 1 == 0)
  stopifnot("Argument S must be a positive integer multiplication of argument thin." = S %% thin == 0)
  stopifnot("Argument show_progress must be a logical value." = is.logical(show_progress))
  
  # call method
  UseMethod("estimate_hyper", specification)
}