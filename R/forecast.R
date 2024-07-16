
#' @title Forecasting using Structural Vector Autoregression
#'
#' @description Samples from the joint predictive density of all of the dependent 
#' variables for models from packages \pkg{bsvars}, \pkg{bsvarSIGNs} or 
#' \pkg{bvarPANELs} at forecast horizons from 1 to \code{horizon} specified as 
#' an argument of the function. Also facilitates forecasting using models with 
#' exogenous variables and conditional forecasting given projected future 
#' trajcetories of (some of the) variables.
#' 
#' @method forecast PosteriorBSVARSIGN
#' @param posterior posterior estimation outcome - an object of class 
#' \code{PosteriorBSVARSIGN} obtained by running the \code{estimate} function.
#' @param horizon a positive integer, specifying the forecasting horizon.
#' @param exogenous_forecast a matrix of dimension \code{horizon x d} containing 
#' forecasted values of the exogenous variables. 
#' @param conditional_forecast a \code{horizon x N} matrix with forecasted values 
#' for selected variables. It should only contain \code{numeric} or \code{NA} 
#' values. The entries with \code{NA} values correspond to the values that are 
#' forecasted conditionally on the realisations provided as \code{numeric} values.
#' 
#' @return A list of class \code{Forecasts} containing the
#' draws from the predictive density and data. The output list includes element:
#' 
#' \describe{
#'  \item{forecasts}{an \code{NxhorizonxS} array with the draws from predictive density}
#'  \item{Y}{an \eqn{NxT} matrix with the data on dependent variables}
#' }
#' 
#' @seealso \code{\link{estimate.BSVARSIGN}}, \code{\link{summary}}, \code{\link{plot}}
#' 
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me} and Xiaolei Wang \email{adamwang15@gmail.com}
#' 
#' @examples
#' # upload data
#' data(oil)
#' 
#' # specify the model and set seed
#' set.seed(123)
#' sign_irf       = array(matrix(c(-1, -1, 1, rep(NA, 6)), nrow = 3), dim = c(3, 3, 1))
#' specification  = specify_bsvarSIGN$new(oil, sign_irf = sign_irf)
#' 
#' # estimate the model
#' posterior      = estimate(specification, 20)
#' 
#' # sample from predictive density 1 year ahead
#' predictive     = forecast(posterior, 4)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' oil |>
#'   specify_bsvarSIGN$new(sign_irf = sign_irf) |>
#'   estimate(S = 20) |> 
#'   forecast(horizon = 4) -> predictive
#' 
#' # conditional forecasting 2 quarters ahead conditioning on 
#' #  provided future values for the Gross Domestic Product 
#' ############################################################
#' cf        = matrix(NA , 2, 3)
#' cf[,3]    = tail(oil, 1)[3]   # conditional forecasts equal to the last gdp observation
#' predictive    = forecast(posterior, 2, conditional_forecast = cf)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' oil |>
#'   specify_bsvarSIGN$new(sign_irf = sign_irf) |>
#'   estimate(S = 10) |> 
#'   forecast(horizon = 2, conditional_forecast = cf) -> predictive
#' 
#' @export
forecast.PosteriorBSVARSIGN = function(
    posterior, 
    horizon = 1, 
    exogenous_forecast = NULL,
    conditional_forecast = NULL
) {
  
  posterior_Sigma = posterior$posterior$Sigma
  posterior_A     = posterior$posterior$A
  T               = ncol(posterior$last_draw$data_matrices$X)
  X_T             = posterior$last_draw$data_matrices$X[,T]
  Y               = posterior$last_draw$data_matrices$Y
  
  N               = nrow(posterior_Sigma)
  K               = length(X_T)
  d               = K - N * posterior$last_draw$p - 1
  S               = dim(posterior_Sigma)[3]
  
  # prepare forecasting with exogenous variables
  if (d == 0 ) {
    exogenous_forecast = matrix(NA, horizon, 1)
  } else {
    stopifnot("Forecasted values of exogenous variables are missing." = (d > 0) & !is.null(exogenous_forecast))
    stopifnot("The matrix of exogenous_forecast does not have a correct number of columns." = ncol(exogenous_forecast) == d)
    stopifnot("Provide exogenous_forecast for all forecast periods specified by argument horizon." = nrow(exogenous_forecast) == horizon)
    stopifnot("Argument exogenous has to be a matrix." = is.matrix(exogenous_forecast) & is.numeric(exogenous_forecast))
    stopifnot("Argument exogenous cannot include missing values." = sum(is.na(exogenous_forecast)) == 0 )
  }
  
  # prepare forecasting with conditional forecasts
  if ( is.null(conditional_forecast) ) {
    # this will not be used for forecasting, but needs to be provided
    conditional_forecast = matrix(NA, horizon, N)
  } else {
    stopifnot("Argument conditional_forecast must be a matrix with numeric values."
              = is.matrix(conditional_forecast) & is.numeric(conditional_forecast)
    )
    stopifnot("Argument conditional_forecast must have the number of rows equal to 
              the value of argument horizon."
              = nrow(conditional_forecast) == horizon
    )
    stopifnot("Argument conditional_forecast must have the number of columns 
              equal to the number of columns in the used data."
              = ncol(conditional_forecast) == N
    )
  }
  
  # perform forecasting
  for_y       = .Call(`_bsvarSIGNs_forecast_bsvarSIGNs`, 
                      posterior_Sigma,
                      posterior_A,
                      X_T,
                      exogenous_forecast,
                      conditional_forecast,
                      horizon
  ) # END .Call
  
  fore            = list()
  fore$forecasts  = for_y
  fore$Y          = Y
  class(fore)     = "Forecasts"
  
  return(fore)
} # END forecast.PosteriorBSVARSIGN
