

#' R6 Class Representing IdentificationBSVARSIGN
#'
#' @description
#' The class IdentificationBSVARSIGN presents the identifying restrictions for the bsvar models with sign and narrative restrictions.
#'
#' @examples 
#' specify_identification_bsvarSIGN$new(N = 3) # recursive specification for a 3-variable system
#' 
#' # an identification pattern with narrative sign restrictions
#' sign_narrative <- matrix(c(2, 0, 3, 2, 236, 0), ncol = 6)
#' specify_identification_bsvarSIGN$new(N = 3, sign_narrative = sign_narrative) 
#'
#' @export
specify_identification_bsvarSIGN = R6::R6Class(
  "IdentificationBSVARSIGN",
  
  public = list(
    
    #' @field VB a list of \code{N} matrices determining the unrestricted elements of matrix \eqn{B}.
    VB       = list(),
    #' @field sign_irf a \code{NxNxh} array of sign restrictions on the impulse response functions.
    sign_irf = array(),
    #' @field sign_narrative a \code{Kx6} matrix of narrative sign restrictions.
    sign_narrative  = matrix(),
    #' @field sign_B a \code{NxN} matrix of sign restrictions on contemporaneous relations.
    sign_B   = matrix(),
    #' @field max_tries a positive integer with the maximum number of iterations 
    #' for finding a rotation matrix \eqn{Q} that would satisfy sign restrictions.
    max_tries = 10000,
    
    #' @description
    #' Create new identifying restrictions IdentificationBSVARSIGN.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param sign_irf a \code{NxNxh} array of sign restrictions on the impulse response functions
    #' up to \code{h} horizons. Each \code{NxN} matrix only accept values in {-1 ,0, 1},
    #' -1: negative restriction, 0: no restriction, 1: positive restriction.
    #' @param sign_narrative a \code{Kx6} matrix of narrative sign restrictions,
    #' each row of the matrix corresponds to a different restriction,
    #' each column contains the following information: \cr
    #' | type | sign | var_i | shock_j | start | horizons | \cr
    #' type: 0 if no restriction; \cr
    #' 1 if restriction on structural shock; \cr
    #' 2 if type A restriction on historical decomposition
    #' i.e. hd of shock_j on var_i is greater (less) than 0; \cr
    #' 3 if type B restriction on historical decomposition
    #' i.e. hd of shock_j on var_i is the largest (smallest); \cr
    #' sign: depending on type, 1 if greater/largest, otherwise less/smallest; \cr
    #' var_i: \code{i}-th variable; \cr
    #' shock_j: \code{j}-th shock; \cr
    #' start: starting period of the restriction; \cr
    #' horizons: horizons of the restriction, input 0 if only one period \cr
    #' @param sign_B a \code{NxN} matrix of sign restrictions on contemporaneous relations.
    #' @param max_tries a positive integer with the maximum number of iterations
    #' for finding a rotation matrix \eqn{Q} that would satisfy sign restrictions.
    #' @return Identifying restrictions IdentificationBSVARSIGN.
    initialize = function(N, sign_irf, sign_narrative, sign_B, max_tries = 10000) {
        
      if (missing(sign_irf)) {
        sign_irf <- array(rep(0, N^2), dim = c(N, N, 1))
      }
      
      if (missing(sign_narrative)) {
        sign_narrative <- matrix(rep(0, 6), ncol = 6, nrow = 1)
      }
      
      if (missing(sign_B)) {
        sign_B <- matrix(rep(0, N^2), ncol = N, nrow = N)
      }
      
      B     = matrix(FALSE, N, N)
      B[lower.tri(B, diag = TRUE)] = TRUE
      
      self$VB          <- vector("list", N)
      for (n in 1:N) {
        self$VB[[n]]   <- matrix(diag(N)[B[n,],], ncol = N)
      }
      
      # TODO: verify sign_irf and sign_narrative
      self$sign_irf <- sign_irf
      self$sign_narrative  <- sign_narrative
      self$sign_B   <- sign_B
      self$max_tries <- max_tries
    }, # END initialize
    
    #' @description
    #' Returns the elements of the identification pattern IdentificationBSVARSIGN as a \code{list}.
    #'
    get_identification = function() {
      list(
        VB       = as.list(self$VB),
        sign_irf = as.array(self$sign_irf),
        sign_narrative  = as.matrix(self$sign_narrative),
        sign_B   = as.matrix(self$sign_B),
        max_tries       = self$max_tries
        )
    }, # END get_identification
    
    #' @description
    #' Set new starting values StartingValuesBSVARSIGN.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param sign_irf a \code{NxNxh} array of sign restrictions on the impulse response functions
    #' up to \code{h} horizons. Each \code{NxN} matrix only accept values in {-1 ,0, 1},
    #' -1: negative restriction, 0: no restriction, 1: positive restriction.
    #' @param sign_narrative a \code{Mx6} matrix of narrative sign restrictions,
    #' each row of the matrix corresponds to a different restriction,
    #' each column contains the following information: \cr
    #' | type | sign | var_i | shock_j | start | horizons | \cr
    #' type: 0 if no restriction; \cr
    #' 1 if restriction on structural shock; \cr
    #' 2 if type A restriction on historical decomposition
    #' i.e. hd of shock_j on var_i is greater (less) than 0; \cr
    #' 3 if type B restriction on historical decomposition
    #' i.e. hd of shock_j on var_i is the largest (smallest); \cr
    #' sign: depending on type, 1 if greater/largest, otherwise less/smallest; \cr
    #' var_i: \code{i}-th variable; \cr
    #' shock_j: \code{j}-th shock; \cr
    #' start: starting period of the restriction; \cr
    #' horizons: horizons of the restriction, input 0 if only one period \cr
    #' @param sign_B a \code{NxN} matrix of sign restrictions on contemporaneous relations.
    #' @param max_tries a positive integer with the maximum number of iterations
    #' for finding a rotation matrix \eqn{Q} that would satisfy sign restrictions.
    set_identification = function(N, sign_irf, sign_narrative, sign_B) {
      B     = matrix(FALSE, N, N)
      B[lower.tri(B, diag = TRUE)] = TRUE
      
      self$VB          <- vector("list", N)
      for (n in 1:N) {
        self$VB[[n]]   <- matrix(diag(N)[B[n,],], ncol = N)
      }
      
      # TODO: verify sign_irf and sign_narrative
      self$sign_irf <- sign_irf
      self$sign_narrative  <- sign_narrative
      self$sign_B   <- sign_B
    } # END set_identification
  ) # END public
) # END specify_identification_bsvarSIGN



#' R6 Class representing the specification of the BSVARSIGN model
#'
#' @description
#' The class BSVARSIGN presents complete specification for the bsvar model with sign and narrative restrictions.
#'
#' @seealso \code{\link{estimate}}, \code{\link{specify_posterior_bsvarSIGN}}
#' 
#' @examples 
#' data(oil)
#' specification = specify_bsvarSIGN$new(
#'    data = oil,
#'    p = 4
#' )
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
    #' Create a new specification of the bsvar model with sign and narrative restrictions BSVARSIGN.
    #' @param data a \code{(T+p)xN} matrix with time series data.
    #' @param p a positive integer providing model's autoregressive lag order.
    #' @param sign_irf a \code{NxNxh} array of sign restrictions on the impulse response functions
    #' up to \code{h} horizons. Each \code{NxN} matrix only accept values in {-1 ,0, 1},
    #' -1: negative restriction, 0: no restriction, 1: positive restriction.
    #' @param sign_narrative a \code{Mx6} matrix of narrative sign restrictions,
    #' each row of the matrix corresponds to a different restriction,
    #' each column contains the following information: \cr
    #' | type | sign | var_i | shock_j | start | horizons | \cr
    #' type: 0 if no restriction; \cr
    #' 1 if restriction on structural shock; \cr
    #' 2 if type A restriction on historical decomposition
    #' i.e. hd of shock_j on var_i is greater (less) than 0; \cr
    #' 3 if type B restriction on historical decomposition
    #' i.e. hd of shock_j on var_i is the largest (smallest); \cr
    #' sign: depending on type, 1 if greater/largest, otherwise less/smallest; \cr
    #' var_i: \code{i}-th variable; \cr
    #' shock_j: \code{j}-th shock; \cr
    #' start: starting period of the restriction; \cr
    #' horizons: horizons of the restriction, input 0 if only one period \cr
    #' @param sign_B a \code{NxN} matrix of sign restrictions on contemporaneous relations.
    #' @param max_tries a positive integer with the maximum number of iterations
    #' for finding a rotation matrix \eqn{Q} that would satisfy sign restrictions.
    #' @param exogenous a \code{(T+p)xd} matrix of exogenous variables.
    #' @param stationary an \code{N} logical vector - its element set to \code{FALSE} sets
    #' the prior mean for the autoregressive parameters of the \code{N}th equation to the white noise process,
    #' otherwise to random walk.
    #' @return A new complete specification for the homoskedastic bsvar model BSVARSIGN.
    initialize = function(
    data,
    p = 1L,
    sign_irf,
    sign_narrative,
    sign_B,
    max_tries = 10000L,
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
      
      if (missing(sign_irf)) {
        sign_irf <- array(rep(0, N^2), dim = c(N, N, 1))
      }
      
      if (missing(sign_narrative)) {
        sign_narrative <- matrix(rep(0, 6), ncol = 6, nrow = 1)
      }
      
      if (missing(sign_B)) {
        sign_B <- matrix(rep(0, N^2), ncol = N, nrow = N)
      }

      # stopifnot("Incorrectly specified argument B." = (is.matrix(B) & is.logical(B)) | (length(B) == 1 & is.na(B)))
      
      B     = matrix(FALSE, N, N)
      B[lower.tri(B, diag = TRUE)] = TRUE
      
      self$data_matrices   = bsvars::specify_data_matrices$new(data, p, exogenous)
      self$identification  = specify_identification_bsvarSIGN$new(N, sign_irf, sign_narrative, sign_B, max_tries)
      self$prior           = bsvars::specify_prior_bsvar$new(N, p, d, stationary)
      self$starting_values = bsvars::specify_starting_values_bsvar$new(N, self$p, d)
    }, # END initialize
    
    #' @description
    #' Returns the data matrices as the DataMatricesBSVAR object.
    #'
    #' @examples 
    #' data(oil)
    #' spec = specify_bsvarSIGN$new(
    #'    data = oil,
    #'    p = 4
    #' )
    #' spec$get_data_matrices()
    #'
    get_data_matrices = function() {
      self$data_matrices$clone()
    }, # END get_data_matrices
    
    #' @description
    #' Returns the identifying restrictions as the IdentificationBSVARSIGN object.
    #'
    #' @examples 
    #' data(oil)
    #' spec = specify_bsvarSIGN$new(
    #'    data = oil,
    #'    p = 4
    #' )
    #' spec$get_identification()
    #'
    get_identification = function() {
      self$identification$clone()
    }, # END get_identification
    
    #' @description
    #' Returns the prior specification as the PriorBSVAR object.
    #'
    #' @examples 
    #' data(oil)
    #' spec = specify_bsvarSIGN$new(
    #'    data = oil,
    #'    p = 4
    #' )
    #' spec$get_prior()
    #'
    get_prior = function() {
      self$prior$clone()
    }, # END get_prior
    
    #' @description
    #' Returns the starting values as the StartingValuesBSVAR object.
    #'
    #' @examples 
    #' data(oil)
    #' spec = specify_bsvarSIGN$new(
    #'    data = oil,
    #'    p = 4
    #' )
    #' spec$get_starting_values()
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
#' the last MCMC draw for the bsvar model with sign and narrative restrictions.
#' Note that due to the thinning of the MCMC output the starting value in element \code{last_draw}
#' might not be equal to the last draw provided in element \code{posterior}.
#'
#' @seealso \code{\link{estimate}}, \code{\link{specify_bsvarSIGN}}
#' 
#' @examples 
#' # This is a function that is used within estimate()
#' data(oil)
#' specification  = specify_bsvarSIGN$new(oil, p = 4)
#' set.seed(123)
#' posterior      = estimate(specification, 50)
#' class(posterior)
#'
#' @export
specify_posterior_bsvarSIGN = R6::R6Class(
  "PosteriorBSVARSIGN",
  
  private = list(
    normalised = TRUE
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
    #' @examples 
    #' data(oil)
    #' specification  = specify_bsvarSIGN$new(oil)
    #' set.seed(123)
    #' estimate       = estimate(specification, 50)
    #' estimate$get_posterior()
    #'
    get_posterior       = function(){
      self$posterior
    }, # END get_posterior
    
    #' @description
    #' Returns an object of class BSVARSIGN with the last draw of the current MCMC run as
    #' the starting value to be passed to the continuation of the MCMC estimation using \code{estimate()}.
    #'
    #' @examples
    #' data(oil)
    #' 
    #' # specify the model and set seed
    #' specification  = specify_bsvarSIGN$new(oil)
    #' set.seed(123)
    #' 
    #' # run the burn-in
    #' burn_in        = estimate(specification, 10)
    #' 
    #' # estimate the model
    #' posterior      = estimate(burn_in, 10)
    #'
    get_last_draw      = function(){
      self$last_draw$clone()
    }, # END get_last_draw
    
    #' @description
    #' Returns \code{TRUE} if the posterior has been normalised using \code{normalise_posterior()}
    #' and \code{FALSE} otherwise.
    #'
    #' @examples
    #' data(oil)
    #' specification  = specify_bsvarSIGN$new(oil)
    #' set.seed(123)
    #' estimate       = estimate(specification, 20)
    #' 
    #' # check normalisation status afterwards
    #' posterior$is_normalised()
    #'
    is_normalised      = function(){
      private$normalised
    } # END is_normalised
    
  ) # END public
) # END specify_posterior_bsvarSIGN
