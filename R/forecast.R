
#' @export
generics::forecast



#' @title R6 Class Representing Forecasts
#'
#' @description
#' R6 class representing draws from the predictive density of a Bayesian
#' Structural Vector Autoregression model.
#'
#' @details
#' The class contains the following objects:
#'
#' \describe{
#'   \item{\code{forecasts}}{An \code{N x horizon x S} array containing draws 
#'   from the predictive density.}
#'   \item{\code{forecast_mean}}{An \code{N x horizon x S} array containing the 
#'   conditional means of the predictive density.}
#'   \item{\code{forecast_covariance}}{An \code{N x N x horizon x S} array 
#'   containing the conditional covariance matrices of the predictive density.}
#'   \item{\code{Y}}{An \code{N x T} matrix containing the data on the dependent 
#'   variables used for estimation.}
#' }
#'
#' The method \code{as_list()} returns the contents of the \code{Forecasts}
#' object as a list.
#'
#' @param output A list containing the forecasting output, including
#' \code{forecasts}, \code{forecast_mean}, and \code{forecast_cov}.
#' @param Y An \code{N x T} matrix containing the data on the dependent variables.
#'
#' @return An object of class \code{Forecasts}.
#'
#' @examples
#' spec = specify_bsvar$new(us_fiscal_lsuw)
#' burn = estimate(spec, 5)
#' post = estimate(burn, 5)
#' fore = forecast(post, 4)
#' apply(fore$forecasts, 1:2, mean) # compute mean forecasts 
#'
#' @export
specify_forecasts = R6::R6Class(
  classname = "Forecasts",
  
  public = list(
    
    #' @field forecasts
    #' An \code{N x horizon x S} numeric array containing draws from the
    #' predictive density.
    forecasts = array(),
    
    #' @field forecast_mean
    #' An \code{N x horizon x S} numeric array containing the conditional
    #' means of the predictive density.
    forecast_mean = array(),
    
    #' @field forecast_covariance
    #' An \code{N x N x horizon x S} numeric array containing the conditional
    #' covariance matrices of the predictive density.
    forecast_covariance = array(),
    
    #' @field Y
    #' An \code{N x T} numeric matrix containing the data on the dependent
    #' variables used for estimation.
    Y = matrix(),
    
    #' @description
    #' Creates a new \code{Forecasts} object from the output of the forecasting
    #' procedure.
    #'
    #' @param output A list containing the forecasting output, including
    #' \code{forecasts}, \code{forecast_mean}, and \code{forecast_cov}.
    #' @param Y An \code{N x T} matrix containing the data on the dependent variables.
    #'
    #' @return An object of class \code{Forecasts}.
    initialize = function(output, Y) {
      
      N       = dim(output$forecasts)[1]
      horizon = dim(output$forecasts)[2]
      S       = dim(output$forecasts)[3]
      
      forecast_covariance = array(
        NA,
        c(N, N, horizon, S)
      )
      
      for (s in seq_len(S)) {
        forecast_covariance[, , , s] = output$forecast_cov[s, ][[1]]
      }
      
      self$forecasts           = output$forecasts
      self$forecast_mean       = output$forecast_mean
      self$forecast_covariance = forecast_covariance
      self$Y                   = Y
      
      invisible(self)
    },
    
    #' @description
    #' Converts the \code{Forecasts} object to a list.
    #'
    #' @return A list containing \code{forecasts}, \code{forecast_mean},
    #' \code{forecast_covariance}, and \code{Y}.
    get_forecasts = function() {
      
      list(
        forecasts           = self$forecasts,
        forecast_mean       = self$forecast_mean,
        forecast_covariance = self$forecast_covariance,
        Y                   = self$Y
      )
    }
  )
)


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
#' @param object posterior estimation outcome - an object of class 
#' \code{PosteriorBSVARSIGN} obtained by running the \code{estimate} function.
#' @param horizon a positive integer, specifying the forecasting horizon.
#' @param exogenous_forecast a matrix of dimension \code{horizon x d} containing 
#' forecasted values of the exogenous variables. 
#' @param conditional_forecast a \code{horizon x N} matrix with forecasted values 
#' for selected variables. It should only contain \code{numeric} or \code{NA} 
#' values. The entries with \code{NA} values correspond to the values that are 
#' forecasted conditionally on the realisations provided as \code{numeric} values.
#' @param ... Not used here.
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
#' # + no effect on productivity (zero restriction)
#' # + positive effect on stock prices (positive sign restriction) 
#' sign_irf       = matrix(c(0, 1, rep(NA, 23)), 5, 5)
#' specification  = specify_bsvarSIGN$new(optimism, sign_irf = sign_irf)
#' 
#' # estimate the model
#' posterior      = estimate(specification, 10)
#' 
#' # sample from predictive density 1 year ahead
#' predictive     = forecast(posterior, 4)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' optimism |>
#'   specify_bsvarSIGN$new(sign_irf = sign_irf) |>
#'   estimate(S = 20) |> 
#'   forecast(horizon = 4) -> predictive
#' 
#' # conditional forecasting 2 quarters ahead conditioning on 
#' #  provided future values for the Gross Domestic Product 
#' ############################################################
#' cf         = matrix(NA , 2, 5)
#' # # conditional forecasts equal to the last consumption observation
#' cf[,3]     = tail(optimism, 1)[3]
#' predictive = forecast(posterior, 2, conditional_forecast = cf)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' optimism |>
#'   specify_bsvarSIGN$new(sign_irf = sign_irf) |>
#'   estimate(S = 10) |> 
#'   forecast(horizon = 2, conditional_forecast = cf) -> predictive
#' 
#' @export
forecast.PosteriorBSVARSIGN = function(
    object, 
    horizon = 1, 
    exogenous_forecast = NULL,
    conditional_forecast = NULL,
    ...
) {
  stopifnot("forecast: specify horizon as integer." = horizon %% 1 == 0)
  
  posterior_Sigma = object$posterior$Sigma
  posterior_A     = object$posterior$A
  T               = ncol(object$last_draw$data_matrices$X)
  X_T             = object$last_draw$data_matrices$X[,T]
  Y               = object$last_draw$data_matrices$Y
  posterior_hyper = object$posterior$hyper
  covid           = object$last_draw$prior$covid
  covid           = ifelse(is.null(covid), -1, covid)
  
  N               = nrow(posterior_Sigma)
  K               = length(X_T)
  d               = K - N * object$last_draw$p - 1
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
  fore        = .Call(`_bsvarSIGNs_forecast_bsvarSIGNs`, 
                      posterior_Sigma,
                      posterior_A,
                      posterior_hyper,
                      X_T,
                      exogenous_forecast,
                      conditional_forecast,
                      covid,
                      T,
                      horizon
  ) # END .Call
  
  output = specify_forecasts$new(fore, Y)
  return(output)
} # END forecast.PosteriorBSVARSIGN
