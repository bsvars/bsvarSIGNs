#' @title Bayesian estimation of the model's hyper-parameters
#'
#' @description Estimates the hyper-parameters of the model including scalar
#' hyper-parameters:
#' \itemize{
#'   \item \eqn{\mu} - the hyper-parameter of the sum-of-coefficients prior 
#'   proposed by Doan, Litterman, Sims (1984). The value of this hyper-parameter 
#'   going to \code{0} implies the presence of a unit root in each equation and 
#'   rules out cointegration, whereas when it goes to infinity the prior becomes 
#'   uninformative.
#'   \item \eqn{\delta} - the hyper-parameter of the sum-of-coefficients prior
#'   proposed by Doan, Litterman, Sims (1984). As the value of this hyper-parameter
#'   goes to zero, all the variables in the system are set at their unconditional 
#'   mean, or the system includes an unspecified number of unit roots without 
#'   drift. This prior is consistent with cointegration.
#'   \item \eqn{\lambda} - the overall prior shrinkage of the autoregressive 
#'    parameters in matrix \eqn{\mathbf{A}}.
#' }
#' as well as: 
#' \itemize{
#'   \item \eqn{\psi} - the \code{N}-vector of the diagonal elements of the 
#'   \code{NxN} scale matrix of inverse Wishart prior distribution for error term
#'   covariance matrix \eqn{\Sigma}.
#' }  
#' following the approach by Giannone, Lenza, & Primiceri (2015).
#' 
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
#' and the prior distribution \eqn{p(\theta)}:
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
#' \eqn{\Omega_{s}} is the \code{(3+N)x(3+N)} covariance of the first \code{s-1} draws of \eqn{\vartheta^*},
#' and \eqn{\sigma_{s}^2} is the adaptive scaling parameter determined by the 
#' dynamic equation:
#' \deqn{\sigma_{s+1} = \sigma_{s} + s^{-0.6}(\alpha_s - 0.234)}
#' where \eqn{\alpha_s} is the Metropolis-Hastings algorithm acceptance rate
#' at the \code{s}th iteration.
#' 
#' @param specification  an object of class \code{BSVARSIGN} generated using the 
#' \code{specify_bsvarSIGN$new()} function.
#' @param S a positive integer specifying the number of MCMC iterations in the 
#' final run after the burn-in stage run to achieve convergence. The total number 
#' of iterations in the Adaptive Random Walk Metropolis–Hastings algorithm is 
#' equal to \code{S + S_burn}, whereas the number of draws returned by the function is 
#' equal to \code{(S + S_burn)/thin}.
#' @param S_burn a positive integer specifying the number of MCMC iterations in the 
#' burn-in run performed to achieve convergence. The total number 
#' of iterations in the Adaptive Random Walk Metropolis–Hastings algorithm is 
#' equal to \code{S + S_burn}, whereas the number of draws returned by the function is 
#' equal to \code{(S + S_burn)/thin}.
#' @param thin  a positive integer, specifying the frequency of MCMC output thinning.
#' @param hyper a logical vector of length \code{4} indicating whether the
#' particular hyper-parameter should be estimated. The ordering of the hyper-parameters
#' is as follows: \eqn{\mu}, \eqn{\delta}, \eqn{\lambda}, \eqn{\psi}. 
#' A value \code{TRUE} for the particular 
#' entry means that the corresponding hyper-parameter will be estimated.
#' An example of such a logical vector that would lead to all the hyper-parameters 
#' being estimated is:
#' \preformatted{
#'     mu  delta lambda   psi
#'   TRUE   TRUE   TRUE  TRUE
#' }
#' If no vector is provided, all the hyper-parameters are estimated.
#' @param show_progress a logical value indicating whether the progress bar for 
#' the procedure should be displayed.
#' 
#' @return An object of class \code{BSVARSIGN} containing the Bayesian estimation 
#' output for the hyper-parameters that can be found in \code{$prior$hyper}.
#' 
#' @seealso \code{\link{specify_bsvarSIGN}}, \code{\link{estimate.BSVARSIGN}}
#'
#' @author Xiaolei Wang \email{adamwang15@gmail.com} and Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' 
#' Andrieu, Moulines (2006) On the Ergodicity Properties of Some Adaptive MCMC 
#' Algorithms. Annals of Applied Probability, 16(3), 1462-1505, 
#' <doi:10.1214/105051606000000286>.
#' 
#' Atchadé, Fort (2010) Limit Theorems for Some Adaptive MCMC Algorithms with 
#' Subgeometric Kernels, Bernoulli, 16(1) 116-154, <doi:10.3150/09-BEJ199>.
#' 
#' Doan, Litterman, Sims (1984) Forecasting and Conditional Projection Using 
#' Realistic Prior Distributions, Econometric Reviews 3, 1–100, 
#' <doi:10.1080/07474938408800053>.
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
#' posterior     = estimate(specification, S = 5)
#' 
#' # workflow with a pipe
#' sign_irf       = matrix(c(0, 1, rep(NA, 23)), 5, 5)
#' (100 * optimism) |> 
#'   specify_bsvarSIGN$new(p = 12, sign_irf = sign_irf) |> 
#'   estimate_hyper(S = 10, S_burn = 5) |> 
#'   estimate(S = 10) -> posterior
#'   
#' @export
estimate_hyper <- function(
                    specification, 
                    S, 
                    S_burn, 
                    thin          = 1,
                    hyper,
                    show_progress = TRUE
                  ) {
  
  # check the inputs
  stopifnot("Argument S must be a positive integer number." = S > 1 & S %% 1 == 0)
  stopifnot("Argument S_burn must be a positive integer number." = S_burn > 1 & S_burn %% 1 == 0)
  stopifnot("Argument thin must be a positive integer number." = thin > 0 & thin %% 1 == 0)
  stopifnot("Argument S must be a positive integer multiplication of argument thin." = S %% thin == 0)
  stopifnot("Argument show_progress must be a logical value." = is.logical(show_progress))
  
  # call method
  UseMethod("estimate_hyper", specification)
}

#' @inherit estimate_hyper
#' 
#' @method estimate_hyper BSVARSIGN
#' 
#' @param specification  an object of class \code{BSVARSIGN} generated using the 
#' \code{specify_bsvarSIGN$new()} function.
#' 
#' @export
estimate_hyper.BSVARSIGN <- function(
    specification, 
    S, 
    S_burn, 
    thin          = 1,
    hyper,
    show_progress = TRUE
) {
  
  if ( missing(hyper) ) {
    hyper = c(TRUE, TRUE, TRUE, TRUE)
  } else {
    stopifnot("Argument hyper must be a logical vector of length 4." = length(hyper) == 4)
  }
  
  specification$prior$estimate_hyper(
    S = S, 
    burn_in = S_burn,
    mu = hyper[1], 
    delta = hyper[2], 
    lambda = hyper[3], 
    psi = hyper[4]
  )
  
  return(specification)
}