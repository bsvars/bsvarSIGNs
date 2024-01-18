
#' R6 Class Representing PriorBSVAR
#'
#' @description
#' The class PriorBSVARSIGN presents a prior specification for the homoskedastic bsvar model.
#' 
#' 
#' @export
specify_prior_bsvarSIGN = R6::R6Class(
  "PriorBSVARSIGN",
  
  public = list(
    
    #' @field A an \code{NxK} matrix, the mean of the normal prior distribution for the parameter matrix \eqn{A}. 
    A          = matrix(),
    
    #' @field A_V_inv a \code{KxK} precision matrix of the normal prior distribution for each of 
    #' the row of the parameter matrix \eqn{A}. This precision matrix is equation invariant.
    A_V_inv    = matrix(),
    
    #' @field B_V_inv an \code{NxN} precision matrix of the generalised-normal prior distribution 
    #' for the structural matrix \eqn{B}. This precision matrix is equation invariant.
    B_V_inv    = matrix(),
    
    #' @field B_nu a positive integer greater of equal than \code{N}, a shape parameter of 
    #' the generalised-normal prior distribution for the structural matrix \eqn{B}.
    B_nu       = NA,
    
    #' @field hyper_nu_B a positive scalar, the shape parameter of the inverted-gamma 2 prior
    #' for the overall shrinkage parameter for matrix \eqn{B}.
    hyper_nu_B = NA,
    
    #' @field hyper_a_B a positive scalar, the shape parameter of the gamma prior
    #' for the second-level hierarchy for the overall shrinkage parameter for matrix \eqn{B}.
    hyper_a_B  = NA,
    
    #' @field hyper_s_BB a positive scalar, the scale parameter of the inverted-gamma 2 prior
    #' for the third-level of hierarchy for overall shrinkage parameter for matrix \eqn{B}.
    hyper_s_BB  = NA,
    
    #' @field hyper_nu_BB a positive scalar, the shape parameter of the inverted-gamma 2 prior
    #' for the third-level of hierarchy for overall shrinkage parameter for matrix \eqn{B}.
    hyper_nu_BB  = NA,
    
    #' @field hyper_nu_A a positive scalar, the shape parameter of the inverted-gamma 2 prior 
    #' for the overall shrinkage parameter for matrix \eqn{A}.
    hyper_nu_A  = NA,
    
    #' @field hyper_a_A a positive scalar, the shape parameter of the gamma prior
    #' for the second-level hierarchy for the overall shrinkage parameter for matrix \eqn{A}.
    hyper_a_A  = NA,
    
    #' @field hyper_s_AA a positive scalar, the scale parameter of the inverted-gamma 2 prior
    #' for the third-level of hierarchy for overall shrinkage parameter for matrix \eqn{A}.
    hyper_s_AA  = NA,
    
    #' @field hyper_nu_AA a positive scalar, the shape parameter of the inverted-gamma 2 prior
    #' for the third-level of hierarchy for overall shrinkage parameter for matrix \eqn{A}.
    hyper_nu_AA  = NA,
    
    #' @description
    #' Create a new prior specification PriorBSVARSIGN.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param p a positive integer - the autoregressive lag order of the SVAR model.
    #' @param d a positive integer - the number of \code{exogenous} variables in the model.
    #' @param stationary an \code{N} logical vector - its element set to \code{FALSE} sets 
    #' the prior mean for the autoregressive parameters of the \code{N}th equation to the white noise process, 
    #' otherwise to random walk.
    #' @return A new prior specification PriorBSVARSIGN.
    #' 
    initialize = function(N, p, d = 0, stationary = rep(FALSE, N)){
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      stopifnot("Argument d must be a non-negative integer number." = d >= 0 & d %% 1 == 0)
      stopifnot("Argument stationary must be a logical vector of length N." = length(stationary) == N & is.logical(stationary))
      
      K                 = N * p + 1 + d
      self$A            = cbind(diag(as.numeric(!stationary)), matrix(0, N, K - N))
      self$A_V_inv      = diag(c(kronecker((1:p)^2, rep(1, N) ), rep(1, d + 1)))
      self$B_V_inv      = diag(N)
      self$B_nu         = N
      self$hyper_nu_B   = 10
      self$hyper_a_B    = 10
      self$hyper_s_BB   = 100
      self$hyper_nu_BB  = 1
      self$hyper_nu_A   = 10
      self$hyper_a_A    = 10
      self$hyper_s_AA   = 10
      self$hyper_nu_AA  = 10
    }, # END initialize
    
    #' @description
    #' Returns the elements of the prior specification PriorBSVARSIGN as a \code{list}.
    #' 
    #' 
    get_prior = function(){
      list(
        A        = self$A,
        A_V_inv  = self$A_V_inv,
        B_V_inv  = self$B_V_inv,
        B_nu     = self$B_nu,
        hyper_nu_B  = self$hyper_nu_B,
        hyper_a_B   = self$hyper_a_B,
        hyper_s_BB  = self$hyper_s_BB,
        hyper_nu_BB = self$hyper_nu_BB,
        hyper_nu_A  = self$hyper_nu_A,
        hyper_a_A   = self$hyper_a_A,
        hyper_s_AA  = self$hyper_s_AA,
        hyper_nu_AA = self$hyper_nu_AA
      )
    } # END get_prior
    
  ) # END public
) # END specify_prior_bsvarSIGN


#' R6 Class Representing StartingValuesBSVARSIGN
#'
#' @description
#' The class StartingValuesBSVARSIGN presents starting values for the homoskedastic bsvar model.
#' 
#' 
#' @export
specify_starting_values_bsvarSIGN = R6::R6Class(
  "StartingValuesBSVARSIGN",
  
  public = list(
    
    #' @field A an \code{NxK} matrix of starting values for the parameter \eqn{A}. 
    A             = matrix(),
    
    #' @field B an \code{NxN} matrix of starting values for the parameter \eqn{B}. 
    B             = matrix(),
    
    #' @field hyper a \code{(2*N+1)x2} matrix of starting values for the shrinkage hyper-parameters of the 
    #' hierarchical prior distribution. 
    hyper         = matrix(),
    
    #' @description
    #' Create new starting values StartingValuesBSVARSIGN.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param p a positive integer - the autoregressive lag order of the SVAR model.
    #' @param d a positive integer - the number of \code{exogenous} variables in the model.
    #' @return Starting values StartingValuesBSVARSIGN.
    #' 
    initialize = function(N, p, d = 0){
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      stopifnot("Argument d must be a non-negative integer number." = d >= 0 & d %% 1 == 0)

      K                   = N * p + 1 + d
      self$B              = diag(N)
      self$A              = cbind(diag(runif(N)), matrix(0, N, K - N))
      self$hyper          = matrix(10, 2 * N + 1, 2)
    }, # END initialize
    
    #' @description
    #' Returns the elements of the starting values StartingValuesBSVARSIGN as a \code{list}.
    #' 
    #' 
    get_starting_values   = function(){
      list(
        B                 = self$B,
        A                 = self$A,
        hyper             = self$hyper
      )
    }, # END get_starting_values
    
    #' @description
    #' Returns the elements of the starting values StartingValuesBSVARSIGN as a \code{list}.
    #' @param last_draw a list containing the last draw of elements \code{B} - an \code{NxN} matrix, 
    #' \code{A} - an \code{NxK} matrix, and \code{hyper} - a vector of 5 positive real numbers.
    #' @return An object of class StartingValuesBSVARSIGN including the last draw of the current MCMC 
    #' as the starting value to be passed to the continuation of the MCMC estimation using \code{estimate()}.
    #' 
    #' 
    set_starting_values   = function(last_draw) {
        self$B            = last_draw$B
        self$A            = last_draw$A
        self$hyper        = last_draw$hyper
    } # END set_starting_values
  ) # END public
) # END specify_starting_values_bsvarSIGN



#' R6 Class Representing IdentificationBSVARSIGNs
#'
#' @description
#' The class IdentificationBSVARSIGNs presents the identifying restrictions for the bsvar models.
#' 
#' 
#' @export
specify_identification_bsvarSIGNs = R6::R6Class(
  "IdentificationBSVARSIGNs",
  
  public = list(
    
    #' @field VB a list of \code{N} matrices determining the unrestricted elements of matrix \eqn{B}. 
    VB    = list(),
    
    #' @description
    #' Create new identifying restrictions IdentificationBSVARSIGNs.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param B a logical \code{NxN} matrix containing value \code{TRUE} for the elements of 
    #' the structural matrix \eqn{B} to be estimated and value \code{FALSE} for exclusion restrictions 
    #' to be set to zero.
    #' @return Identifying restrictions IdentificationBSVARSIGNs.
    initialize = function(N, B) {
      if (missing(B)) {
          B     = matrix(FALSE, N, N)
          B[lower.tri(B, diag = TRUE)] = TRUE
      }

      stopifnot("Argument B must be an NxN matrix with logical values." = is.logical(B) & is.matrix(B) & prod(dim(B) == N))
      
      self$VB          <- vector("list", N)
      for (n in 1:N) {
        self$VB[[n]]   <- matrix(diag(N)[B[n,],], ncol = N)
      }
    }, # END initialize
    
    #' @description
    #' Returns the elements of the identification pattern IdentificationBSVARSIGNs as a \code{list}.
    #' 
    #' 
    get_identification = function() {
      as.list(self$VB)
    }, # END get_identification
    
    #' @description
    #' Set new starting values StartingValuesBSVARSIGN.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param B a logical \code{NxN} matrix containing value \code{TRUE} for the elements of 
    #' the structural matrix \eqn{B} to be estimated and value \code{FALSE} for exclusion restrictions 
    #' to be set to zero.
    #' 
    set_identification = function(N, B) {
      if (missing(B)) {
        B     = matrix(FALSE, N, N)
        B[lower.tri(B, diag = TRUE)] = TRUE
      }
      
      stopifnot("Argument B must be an NxN matrix with logical values." = is.logical(B) & is.matrix(B) & prod(dim(B) == N))
      
      self$VB          <- vector("list", N)
      for (n in 1:N) {
        self$VB[[n]]   <- matrix(diag(N)[B[n,],], ncol = N)
      }
    } # END set_identification
  ) # END public
) # END specify_identification_bsvarSIGNs



#' R6 Class Representing DataMatricesBSVARSIGN
#'
#' @description
#' The class DataMatricesBSVARSIGN presents the data matrices of dependent variables, \eqn{Y}, 
#' and regressors, \eqn{X}, for the homoskedastic bsvar model.
#' 
#' 
#' @export
specify_data_matrices = R6::R6Class(
  "DataMatricesBSVARSIGN",
  
  public = list(
    
    #' @field Y an \code{NxT} matrix of dependent variables, \eqn{Y}. 
    Y     = matrix(),
    
    #' @field X an \code{KxT} matrix of regressors, \eqn{X}. 
    X     = matrix(),
    
    #' @description
    #' Create new data matrices DataMatricesBSVARSIGN.
    #' @param data a \code{(T+p)xN} matrix with time series data.
    #' @param p a positive integer providing model's autoregressive lag order.
    #' @param exogenous a \code{(T+p)xd} matrix of exogenous variables. 
    #' This matrix should not include a constant term.
    #' @return New data matrices DataMatricesBSVARSIGN.
    initialize = function(data, p = 1L, exogenous = NULL) {
      if (missing(data)) {
        stop("Argument data has to be specified")
      } else {
        stopifnot("Argument data has to be a matrix." = is.matrix(data) & is.numeric(data))
        stopifnot("Argument data has to contain at least 2 columns and 3 rows." = (ncol(data) >= 2 & nrow(data) >= 3))
        stopifnot("Argument data cannot include missing values." = sum(is.na(data)) == 0 )
      }
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      
      if (is.null(exogenous)) {
        d = 0
      } else {
        stopifnot("Argument exogenous has to be a matrix." = is.matrix(exogenous) & is.numeric(exogenous))
        stopifnot("Argument exogenous has to contain at the same number of rows as argument data." = (ncol(exogenous) >= 1 & nrow(data) == nrow(exogenous)))
        stopifnot("Argument exogenous cannot include missing values." = sum(is.na(exogenous)) == 0 )
        d = ncol(exogenous)
        Td = nrow(exogenous)
        test_exogenous = 0
        for (i in 1:d) {
          test_exogenous = test_exogenous + prod(exogenous[,i] - mean(exogenous[,i]) == rep(0,Td))
        }
        stopifnot("Argument exogenous cannot include a constant term." = test_exogenous == 0 )
      }
      
      TT            = nrow(data)
      T             = TT - p
      
      self$Y        = t(data[(p + 1):TT,])
      X             = matrix(0, T, 0)
      for (i in 1:p) {
        X           = cbind(X, data[(p + 1):TT - i,])
      }
      X             = cbind(X, rep(1, T))
      if (!is.null(data)) {
        X           = cbind(X, exogenous[(p + 1):TT,])
      }
      self$X        = t(X)
    }, # END initialize
    
    #' @description
    #' Returns the data matrices DataMatricesBSVARSIGN as a \code{list}.
    #' 
    #' 
    get_data_matrices = function() {
      list(
        Y = self$Y,
        X = self$X
      )
    } # END get_data_matrices
  ) # END public
) # END specify_data_matrices



#' R6 Class representing the specification of the homoskedastic BSVARSIGN model
#'
#' @description
#' The class BSVARSIGN presents complete specification for the homoskedastic bsvar model.
#' 
#' 
#' 
#' @export
specify_bsvarSIGN = R6::R6Class(
  "BSVARSIGN",
  
  public = list(
    
    #' @field p a non-negative integer specifying the autoregressive lag order of the model. 
    p                      = numeric(),
    
    #' @field identification an object IdentificationBSVARSIGN with the identifying restrictions. 
    identification         = list(),
    
    #' @field prior an object PriorBSVARSIGN with the prior specification. 
    prior                  = list(),
    
    #' @field data_matrices an object DataMatricesBSVARSIGN with the data matrices.
    data_matrices          = list(),
    
    #' @field starting_values an object StartingValuesBSVARSIGN with the starting values.
    starting_values        = list(),
    
    #' @description
    #' Create a new specification of the homoskedastic bsvar model BSVARSIGN.
    #' @param data a \code{(T+p)xN} matrix with time series data.
    #' @param p a positive integer providing model's autoregressive lag order.
    #' @param B a logical \code{NxN} matrix containing value \code{TRUE} for the elements of 
    #' the structural matrix \eqn{B} to be estimated and value \code{FALSE} for exclusion restrictions 
    #' to be set to zero.
    #' @param exogenous a \code{(T+p)xd} matrix of exogenous variables. 
    #' @param stationary an \code{N} logical vector - its element set to \code{FALSE} sets 
    #' the prior mean for the autoregressive parameters of the \code{N}th equation to the white noise process, 
    #' otherwise to random walk.
    #' @return A new complete specification for the homoskedastic bsvar model BSVARSIGN.
    initialize = function(
      data,
      p = 1L,
      B,
      exogenous = NULL,
      stationary = rep(FALSE, ncol(data))
    ) {
      stopifnot("Argument p has to be a positive integer." = ((p %% 1) == 0 & p > 0))
      self$p     = p
      
      TT            = nrow(data)
      T             = TT - self$p
      N             = ncol(data)
      d             = 0
      if (!is.null(exogenous)) {
        d           = ncol(exogenous)
      }
      
      if (missing(B)) {
        message("The identification is set to the default option of lower-triangular structural matrix.")
        B     = matrix(FALSE, N, N)
        B[lower.tri(B, diag = TRUE)] = TRUE
      }
      stopifnot("Incorrectly specified argument B." = (is.matrix(B) & is.logical(B)) | (length(B) == 1 & is.na(B)))
      
      self$data_matrices   = specify_data_matrices$new(data, p, exogenous)
      self$identification  = specify_identification_bsvarSIGNs$new(N, B)
      self$prior           = specify_prior_bsvarSIGN$new(N, p, d, stationary)
      self$starting_values = specify_starting_values_bsvarSIGN$new(N, self$p, d)
    }, # END initialize
    
    #' @description
    #' Returns the data matrices as the DataMatricesBSVARSIGN object.
    #' 
    #' 
    get_data_matrices = function() {
      self$data_matrices$clone()
    }, # END get_data_matrices
    
    #' @description
    #' Returns the identifying restrictions as the IdentificationBSVARSIGNs object.
    #' 
    #' 
    get_identification = function() {
      self$identification$clone()
    }, # END get_identification
    
    #' @description
    #' Returns the prior specification as the PriorBSVARSIGN object.
    #' 
    #' 
    get_prior = function() {
      self$prior$clone()
    }, # END get_prior
    
    #' @description
    #' Returns the starting values as the StartingValuesBSVARSIGN object.
    #' 
    #' 
    get_starting_values = function() {
      self$starting_values$clone()
    } # END get_starting_values
  ) # END public
) # END specify_bsvarSIGN



#' R6 Class Representing PosteriorBSVARSIGN
#'
#' @description
#' The class PosteriorBSVARSIGN contains posterior output and the specification including 
#' the last MCMC draw for the homoskedastic bsvar model. 
#' Note that due to the thinning of the MCMC output the starting value in element \code{last_draw}
#' might not be equal to the last draw provided in element \code{posterior}.
#' 
#' 
#' 
#' @export
specify_posterior_bsvarSIGN = R6::R6Class(
  "PosteriorBSVARSIGN",
    
  private = list(
    normalised = FALSE
  ), # END private
  
  public = list(
    
    #' @field last_draw an object of class BSVARSIGN with the last draw of the current MCMC run as 
    #' the starting value to be passed to the continuation of the MCMC estimation using \code{estimate()}. 
    last_draw = list(),
    
    #' @field posterior a list containing Bayesian estimation output collected in elements 
    #' an \code{NxNxS} array \code{B}, an \code{NxKxS} array \code{A}, and a \code{5xS} matrix \code{hyper}.
    posterior = list(),
    
    #' @description
    #' Create a new posterior output PosteriorBSVARSIGN.
    #' @param specification_bsvarSIGN an object of class BSVARSIGN with the last draw of the current 
    #' MCMC run as the starting value.
    #' @param posterior_bsvarSIGN a list containing Bayesian estimation output collected in elements 
    #' an \code{NxNxS} array \code{B}, an \code{NxKxS} array \code{A}, and a \code{5xS} matrix \code{hyper}.
    #' @return A posterior output PosteriorBSVARSIGN.
    initialize = function(specification_bsvarSIGN, posterior_bsvarSIGN) {
      
      stopifnot("Argument specification_bsvarSIGN must be of class BSVARSIGN." = any(class(specification_bsvarSIGN) == "BSVARSIGN"))
      stopifnot("Argument posterior_bsvarSIGN must must contain MCMC output." = is.list(posterior_bsvarSIGN) & is.array(posterior_bsvarSIGN$B) & is.array(posterior_bsvarSIGN$A) & is.array(posterior_bsvarSIGN$hyper))
      
      self$last_draw    = specification_bsvarSIGN
      self$posterior    = posterior_bsvarSIGN
    }, # END initialize
    
    #' @description
    #' Returns a list containing Bayesian estimation output collected in elements 
    #' an \code{NxNxS} array \code{B}, an \code{NxKxS} array \code{A}, and a \code{5xS} matrix \code{hyper}.
    #' 
    #' 
    get_posterior       = function(){
      self$posterior
    }, # END get_posterior
    
    #' @description
    #' Returns an object of class BSVARSIGN with the last draw of the current MCMC run as 
    #' the starting value to be passed to the continuation of the MCMC estimation using \code{estimate()}.
    #' 
    #' 
    get_last_draw      = function(){
      self$last_draw$clone()
    }, # END get_last_draw
    
    #' @description
    #' Returns \code{TRUE} if the posterior has been normalised using \code{normalise_posterior()} 
    #' and \code{FALSE} otherwise.
    #' 
    #' 
    is_normalised      = function(){
      private$normalised
    }, # END is_normalised
    
    #' @description
    #' Sets the private indicator \code{normalised} to TRUE.
    #' @param value (optional) a logical value to be passed to indicator \code{normalised}.
    #' 
    #' 
    set_normalised     = function(value){
      if (missing(value)) {
        private$normalised <- TRUE
      } else {
        private$normalised <- value
      }
    } # END set_normalised

  ) # END public
) # END specify_posterior_bsvarSIGN
