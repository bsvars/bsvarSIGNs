
# construction Z_j matrices
get_Z = function(zero_irf) {
  if (sum(zero_irf) == 0) {
    return(NULL)
  }
  
  if (!(all(zero_irf %in% c(0, 1)))) {
    stop("Zero restriction matrix has entries that are not in {0, 1}.")
  }
  
  N = dim(zero_irf)[2]
  
  Z = list()
  for (j in 1:N) {
    Z_j          = diag(zero_irf[, j])  
    nonzero_rows = rowSums(Z_j) > 0
    Z_j          = as.matrix(Z_j[nonzero_rows, ])
    
    if (dim(Z_j)[2] == 1) {
      Z_j = as.matrix(t(Z_j))
    }
    
    if (dim(Z_j)[1] > N-j) {
      stop("Too many zero restrictions for shock ", j)
    }
    
    Z[[j]]       = Z_j
  }
  
  Z
}

# verify if the matrix is NxN and has entries in {-1, 0, 1}
verify_traditional = function(N, A) {
  if (!(is.matrix(A) && all(dim(A) == c(N, N)))) {
    stop("Sign restriction matrix is not NxN.")
  }
  if (!(all(A %in% c(-1, 0, 1)))) {
    stop("Sign restriction matrix has entries that are not in {-1, 0, 1}.")
  }
}

# verify if narrative sign restriciton matrix is valid, i.e.
# 1. has 6 columns
# 2. each column satisfies its definition
verify_narrative = function(N, A) {
  if (!(is.matrix(A) && dim(A)[2] == 6)) {
    stop("Narrative sign restriction matrix does not have exactly 6 columns.")
  }
  if (!(all(A[,1] %in% c(0, 1, 2, 3)))) {
    stop("Narrative sign restriction matrix column 1 (type) has entries that are not in {0, 1, 2, 3}.")
  }
  if (!(all(A[,2] %in% c(-1, 1)))) {
    stop("Narrative sign restriction matrix column 2 (sign) has entries that are not in {-1, 1}.")
  }
  if (!(all(A[,3] %in% 1:N))) {
    stop("Narrative sign restriction matrix column 3 (var_i) has entries that are not in 1:N.")
  }
  if (!(all(A[,4] %in% 1:N))) {
    stop("Narrative sign restriction matrix column 4 (shock_j) has entries that are not in 1:N.")
  }
  if (!(all(A[,5] == floor(A[,5])))) {
    stop("Narrative sign restriction matrix column 5 (start_t) has entries that are not in 1:T.")
  }
  if (!(all(A[,6] == floor(A[,6])))) {
    stop("Narrative sign restriction matrix column 6 (horizons_h) has entries that are not in 1:(T-start).")
  }
}

# verify all restrictions
verify_all = function(N, sign_irf, sign_narrative, sign_B) {
  verify_traditional(N, sign_B)
  
  verify_narrative(N, sign_narrative)
  
  dim_irf = dim(sign_irf)
  if (length(dim_irf) != 3) {
    stop("Sign restriction array is not 3-dimensional.")
  }
  for (h in 1:dim(sign_irf)[3]) {
    verify_traditional(N, sign_irf[,,h])
  }
}


# Minnesota prior of Normal-Inverse-Wishart form
niw_prior <- function(Y,
                      p,
                      non_stationary,
                      lambda_1 = 1 / .Machine$double.eps,
                      lambda_2 = 0.2,
                      lambda_3 = 1) {
  T <- nrow(Y)
  N <- ncol(Y)
  K <- 1 + N * p
  
  ar_s2 <- vector(length = N)
  for (n in 1:N) {
    resid <- stats::ar(Y[, n], order.max = p)$resid |> stats::na.omit()
    ar_s2[n] <- t(resid) %*% resid / (T - p - 1)
  }
  
  B <- matrix(0, K, N)
  B[1:N, 1:N] <- diag(non_stationary)
  
  V <- matrix(0, K, K)
  V[K, K] <- lambda_1
  V[1:(K-1), 1:(K-1)] <- diag(lambda_2^2 * rep((1:p)^-(2 * lambda_3), each = N) * rep(ar_s2^-1, p))
  
  S <- diag(ar_s2)
  nu <- N + 2
  
  list(B = B, V = V, S = S, nu = nu)
}


#' R6 Class Representing PriorBSVAR
#'
#' @description
#' The class PriorBSVARSIGN presents a prior specification for the homoskedastic bsvar model.
#' 
#' @examples
#' # a prior for 3-variable example with one lag 
#' data(oil)
#' prior = specify_prior_bsvarSIGN$new(oil, N = 3, p = 1)
#' prior$B                                        # show autoregressive prior mean
#' 
#' @export
specify_prior_bsvarSIGN = R6::R6Class(
  "PriorBSVARSIGN",
  
  public = list(
    
    #' @field B an \code{KxN} matrix, the mean of the normal prior distribution for the parameter matrix \eqn{B}. 
    B          = matrix(),
    
    #' @field V a \code{KxK} covariance matrix of the normal prior distribution for each of 
    #' the column of the parameter matrix \eqn{B}. This covariance matrix is equation invariant.
    V          = matrix(),
    
    #' @field S an \code{NxN} scale matrix of the inverse-Wishart prior distribution 
    #' for the covariance matrix \eqn{\Sigma}. This scale matrix is equation invariant.
    S          = matrix(),
    
    #' @field nu a positive real number greater of equal than \code{N}, a degree of freedom parameter
    #' of the inverse-Wishart prior distribution for the covariance matrix \eqn{\Sigma}.
    nu         = NA,
    
    #' @description
    #' Create a new prior specification PriorBSVAR.
    #' @param data the \code{TxN} data matrix.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param p a positive integer - the autoregressive lag order of the SVAR model.
    #' @param d a positive integer - the number of \code{exogenous} variables in the model.
    #' @param stationary an \code{N} logical vector - its element set to \code{FALSE} sets 
    #' the prior mean for the autoregressive parameters of the \code{N}th equation to the white noise process, 
    #' otherwise to random walk.
    #' @return A new prior specification PriorBSVARSIGN.
    #' @examples 
    #' # a prior for 3-variable example with one lag and stationary data
    #' data(oil)
    #' prior = specify_prior_bsvarSIGN$new(oil, N = 3, p = 1, stationary = rep(TRUE, 3))
    #' prior$B # show autoregressive prior mean
    #' 
    initialize = function(data, N, p, d = 0, stationary = rep(FALSE, N)){
      stopifnot("Argument N must be a positive integer number." = N > 0 & N %% 1 == 0)
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      stopifnot("Argument d must be a non-negative integer number." = d >= 0 & d %% 1 == 0)
      stopifnot("Argument stationary must be a logical vector of length N." = length(stationary) == N & is.logical(stationary))
      
      prior   = niw_prior(data, p, !stationary)
      
      self$B  = prior$B
      self$V  = prior$V
      self$S  = prior$S
      self$nu = prior$nu
      
    }, # END initialize
    
    #' @description
    #' Returns the elements of the prior specification PriorBSVAR as a \code{list}.
    #' 
    #' @examples 
    #' # a prior for 3-variable example with four lags
    #' prior = specify_prior_bsvar$new(N = 3, p = 4)
    #' prior$get_prior() # show the prior as list
    #' 
    get_prior = function(){
      list(
        B  = self$B,
        V  = self$V,
        S  = self$S,
        nu = self$nu
      )
    } # END get_prior
    
  ) # END public
) # END specify_prior_bsvarSIGN


#' R6 Class Representing IdentificationBSVARSIGN
#'
#' @description
#' The class IdentificationBSVARSIGN presents the identifying restrictions for the Bayesian Structural VAR models with sign and narrative restrictions.
#'
#' @examples 
#' specify_identification_bsvarSIGN$new(N = 3) # recursive specification for a 3-variable system
#' 
#' # an identification pattern with narrative sign restrictions
#' sign_narrative <- matrix(c(2, -1, 3, 2, 236, 0), ncol = 6)
#' specify_identification_bsvarSIGN$new(N = 3, sign_narrative = sign_narrative) 
#'
#' @export
specify_identification_bsvarSIGN = R6::R6Class(
  "IdentificationBSVARSIGN",
  
  public = list(
    
    #' @field VB a list of \code{N} matrices determining the unrestricted elements of matrix \eqn{B}.
    VB       = list(),
    #' @field sign_irf a \code{NxNxH} array of sign restrictions on the impulse response functions.
    sign_irf = array(),
    #' @field sign_narrative a \code{Mx6} matrix of narrative sign restrictions.
    sign_narrative  = matrix(),
    #' @field sign_B a \code{NxN} matrix of sign restrictions on contemporaneous relations.
    sign_B   = matrix(),
    #' @field zero_irf a \code{NxNxH} array of zero restrictions on the impulse response functions.
    zero_irf = array(),
    #' @field max_tries a positive integer with the maximum number of iterations 
    #' for finding a rotation matrix \eqn{Q} that would satisfy sign restrictions.
    max_tries = 1,
    
    #' @description
    #' Create new identifying restrictions IdentificationBSVARSIGN.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param sign_irf a \code{NxNxH} array with entries in (-1 ,0, 1), sign restrictions on the
    #' impulse response functions, the \code{h}-th slice \code{NxN} matrix contains the
    #' sign restrictions on the \code{h-1} horizon, e.g. \code{sign_irf[,,0]} contains restrictions
    #' on the contemporaneous impulse response function.
    #' @param sign_narrative a \code{Mx6} matrix of narrative sign restrictions,
    #' each row of the matrix corresponds to a different restriction,
    #' columns are (type, sign, var_i, shock_j, start_t, horizons_h) with detailed definitions: \cr\cr
    #' Column 1 (type):
    #' 0 if no restriction;
    #' 1 if restriction on structural shock;
    #' 2 if type A restriction on historical decomposition
    #' i.e. historical decomposition of shock_j on var_i is greater (less) than 0;
    #' 3 if type B restriction on historical decomposition
    #' i.e. historical decomposition of shock_j on var_i is the largest (smallest); \cr
    #' Column 2 (sign): depending on type, 1 if greater/largest, -1 if less/smallest. \cr
    #' Column 3 (var_i): an integer in 1:N, index of the restricted variable. \cr
    #' Column 4 (shock_j): an integer in 1:N, index of the restricted shock. \cr
    #' Column 5 (start_t): an integer in 1:T, starting period of the restriction; \cr
    #' Column 6 (horizons_h): an integer in 1:(T-start_t), number horizons of the restriction,
    #' if start=t and horizons=h the restriction in on periods t to t+h,
    #' e.g. when h=0 the restriction in only placed on period t.
    #' @param sign_B a \code{NxN} matrix with entries in (-1 ,0, 1), sign restrictions on the
    #' contemporaneous relations \code{B} between reduced-form errors \code{E} and
    #' structural shocks \code{U}. Recall the structural equation \code{BE=U}, the inverse
    #' of \code{B} is the contemporaneous impulse response function.
    #' @param zero_irf a \code{NxN} matrix with entries in (0, 1), zero restrictions on the
    #' contemporaneous impulse response functions.
    #' @param max_tries a positive integer with the maximum number of iterations
    #' for finding a rotation matrix \eqn{Q} that would satisfy sign restrictions.
    #' @return Identifying restrictions IdentificationBSVARSIGN.
    initialize = function(N, sign_irf, sign_narrative, sign_B, zero_irf, max_tries = 1) {
        
      missing_all   = TRUE
      if (missing(sign_irf)) {
        sign_irf = array(rep(0, N^2), dim = c(N, N, 1))
      } else {
        missing_all = FALSE
      }
      if (missing(sign_narrative)) {
        sign_narrative = matrix(c(0, 1, 1, 1, 1, 0), ncol = 6, nrow = 1)
      } else {
        missing_all = FALSE
      }
      if (missing(sign_B)) {
        if (missing_all) {
          sign_B = diag(N)
        } else {
          sign_B = matrix(rep(0, N^2), ncol = N, nrow = N)  
        }
      }
      if (missing(zero_irf)) {
        zero_irf = matrix(rep(0, N^2), ncol = N, nrow = N)
      }
      verify_all(N, sign_irf, sign_narrative, sign_B)
      
      B     = matrix(FALSE, N, N)
      B[lower.tri(B, diag = TRUE)] = TRUE
      
      self$VB = vector("list", N)
      for (n in 1:N) {
        self$VB[[n]] = matrix(diag(N)[B[n,],], ncol = N)
      }
      
      self$sign_irf       = sign_irf
      self$sign_narrative = sign_narrative
      self$sign_B         = sign_B
      self$zero_irf       = get_Z(zero_irf)
      self$max_tries      = max_tries
    }, # END initialize
    
    #' @description
    #' Returns the elements of the identification pattern IdentificationBSVARSIGN as a \code{list}.
    #'
    get_identification = function() {
      list(
        VB             = as.list(self$VB),
        sign_irf       = as.array(self$sign_irf),
        sign_narrative = as.matrix(self$sign_narrative),
        sign_B         = as.matrix(self$sign_B),
        zero_irf       = as.list(self$zero_irf),
        max_tries      = self$max_tries
        )
    }, # END get_identification
    
    #' @description
    #' Set new starting values StartingValuesBSVARSIGN.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param sign_irf a \code{NxNxH} array with entries in (-1 ,0, 1), sign restrictions on the
    #' impulse response functions, the \code{h}-th slice \code{NxN} matrix contains the
    #' sign restrictions on the \code{h-1} horizon, e.g. \code{sign_irf[,,0]} contains restrictions
    #' on the contemporaneous impulse response function.
    #' @param sign_narrative a \code{Mx6} matrix of narrative sign restrictions,
    #' each row of the matrix corresponds to a different restriction,
    #' columns are (type, sign, var_i, shock_j, start_t, horizons_h) with detailed definitions: \cr\cr
    #' Column 1 (type):
    #' 0 if no restriction;
    #' 1 if restriction on structural shock;
    #' 2 if type A restriction on historical decomposition
    #' i.e. historical decomposition of shock_j on var_i is greater (less) than 0;
    #' 3 if type B restriction on historical decomposition
    #' i.e. historical decomposition of shock_j on var_i is the largest (smallest); \cr
    #' Column 2 (sign): depending on type, 1 if greater/largest, -1 if less/smallest. \cr
    #' Column 3 (var_i): an integer in 1:N, index of the restricted variable. \cr
    #' Column 4 (shock_j): an integer in 1:N, index of the restricted shock. \cr
    #' Column 5 (start_t): an integer in 1:T, starting period of the restriction; \cr
    #' Column 6 (horizons_h): an integer in 1:(T-start_t), number horizons of the restriction,
    #' if start=t and horizons=h the restriction in on periods t to t+h,
    #' e.g. when h=0 the restriction in only placed on period t.
    #' @param sign_B a \code{NxN} matrix with entries in (-1 ,0, 1), sign restrictions on the
    #' contemporaneous relations \code{B} between reduced-form errors \code{E} and
    #' structural shocks \code{U}. Recall the structural equation \code{BE=U}, the inverse
    #' of \code{B} is the contemporaneous impulse response function.
    #' @param zero_irf a \code{NxN} matrix with entries in (0, 1), zero restrictions on the
    #' contemporaneous impulse response functions.
    #' @param max_tries a positive integer with the maximum number of iterations
    #' for finding a rotation matrix \eqn{Q} that would satisfy sign restrictions.
    set_identification = function(N, sign_irf, sign_narrative, sign_B, zero_irf) {
      B     = matrix(FALSE, N, N)
      B[lower.tri(B, diag = TRUE)] = TRUE
      
      self$VB          <- vector("list", N)
      for (n in 1:N) {
        self$VB[[n]]   <- matrix(diag(N)[B[n,],], ncol = N)
      }
      
      missing_all   = TRUE
      if (missing(sign_irf)) {
        sign_irf = array(rep(0, N^2), dim = c(N, N, 1))
      } else {
        missing_all = FALSE
      }
      if (missing(sign_narrative)) {
        sign_narrative = matrix(c(0, 1, 1, 1, 1, 0), ncol = 6, nrow = 1)
      } else {
        missing_all = FALSE
      }
      if (missing(sign_B)) {
        if (missing_all) {
          sign_B = diag(N)
        } else {
          sign_B = matrix(rep(0, N^2), ncol = N, nrow = N)  
        }
      }
      if (missing(zero_irf)) {
        zero_irf = matrix(rep(0, N^2), ncol = N, nrow = N)
      }
      verify_all(N, sign_irf, sign_narrative, sign_B)
      
      self$sign_irf       = sign_irf
      self$sign_narrative = sign_narrative
      self$sign_B         = sign_B
      self$zero_irf       = get_Z(zero_irf)
    } # END set_identification
  ) # END public
) # END specify_identification_bsvarSIGN



#' R6 Class representing the specification of the BSVARSIGN model
#'
#' @description
#' The class BSVARSIGN presents complete specification for the Bayesian Structural VAR model with sign and narrative restrictions.
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
    #' Create a new specification of the Bayesian Structural VAR model with sign and narrative restrictions BSVARSIGN.
    #' @param data a \code{(T+p)xN} matrix with time series data.
    #' @param p a positive integer providing model's autoregressive lag order.
    #' @param sign_irf a \code{NxNxH} array with entries in (-1 ,0, 1), sign restrictions on the
    #' impulse response functions, the \code{h}-th slice \code{NxN} matrix contains the
    #' sign restrictions on the \code{h-1} horizon, e.g. \code{sign_irf[,,0]} contains restrictions
    #' on the contemporaneous impulse response function.
    #' @param sign_narrative a \code{Mx6} matrix of narrative sign restrictions,
    #' each row of the matrix corresponds to a different restriction,
    #' columns are (type, sign, var_i, shock_j, start_t, horizons_h) with detailed definitions: \cr\cr
    #' Column 1 (type):
    #' 0 if no restriction;
    #' 1 if restriction on structural shock;
    #' 2 if type A restriction on historical decomposition
    #' i.e. historical decomposition of shock_j on var_i is greater (less) than 0;
    #' 3 if type B restriction on historical decomposition
    #' i.e. historical decomposition of shock_j on var_i is the largest (smallest); \cr
    #' Column 2 (sign): depending on type, 1 if greater/largest, -1 if less/smallest. \cr
    #' Column 3 (var_i): an integer in 1:N, index of the restricted variable. \cr
    #' Column 4 (shock_j): an integer in 1:N, index of the restricted shock. \cr
    #' Column 5 (start_t): an integer in 1:T, starting period of the restriction; \cr
    #' Column 6 (horizons_h): an integer in 1:(T-start_t), number horizons of the restriction,
    #' if start=t and horizons=h the restriction in on periods t to t+h,
    #' e.g. when h=0 the restriction in only placed on period t.
    #' @param sign_B a \code{NxN} matrix with entries in (-1 ,0, 1), sign restrictions on the
    #' contemporaneous relations \code{B} between reduced-form errors \code{E} and
    #' structural shocks \code{U}. Recall the structural equation \code{BE=U}, the inverse
    #' of \code{B} is the contemporaneous impulse response function.
    #' @param zero_irf a \code{NxN} matrix with entries in (0, 1), zero restrictions on the
    #' contemporaneous impulse response functions.
    #' @param max_tries a positive integer with the maximum number of iterations
    #' for finding a rotation matrix \eqn{Q} that would satisfy sign restrictions.
    #' @param exogenous a \code{(T+p)xd} matrix of exogenous variables.
    #' @param stationary an \code{N} logical vector - its element set to \code{FALSE} sets
    #' the prior mean for the autoregressive parameters of the \code{N}th equation to the white noise process,
    #' otherwise to random walk.
    #' @return A new complete specification for the Bayesian Structural VAR model BSVARSIGN.
    initialize = function(
    data,
    p = 1L,
    sign_irf,
    sign_narrative,
    sign_B,
    zero_irf,
    max_tries = 1,
    exogenous = NULL,
    stationary = rep(FALSE, ncol(data))
    ) {
      stopifnot("Argument p has to be a positive integer." = ((p %% 1) == 0 & p > 0))
      self$p        = p
      
      TT            = nrow(data)
      T             = TT - self$p
      N             = ncol(data)
      d             = 0
      if (!is.null(exogenous)) {
        d           = ncol(exogenous)
      }
      
      missing_all   = TRUE
      if (missing(sign_irf)) {
        sign_irf = array(rep(0, N^2), dim = c(N, N, 1))
      } else {
        missing_all = FALSE
      }
      if (missing(sign_narrative)) {
        sign_narrative = matrix(c(0, 1, 1, 1, 1, 0), ncol = 6, nrow = 1)
      } else {
        missing_all = FALSE
      }
      if (missing(sign_B)) {
        if (missing_all) {
          sign_B = diag(N)
        } else {
          sign_B = matrix(rep(0, N^2), ncol = N, nrow = N)  
        }
      }
      if (missing(zero_irf)) {
        zero_irf = matrix(rep(0, N^2), ncol = N, nrow = N)
      }
      verify_all(N, sign_irf, sign_narrative, sign_B)
      
      B                            = matrix(FALSE, N, N)
      B[lower.tri(B, diag = TRUE)] = TRUE
      
      self$data_matrices           = bsvars::specify_data_matrices$new(data, p, exogenous)
      self$identification          = specify_identification_bsvarSIGN$new(N,
                                                                          sign_irf,
                                                                          sign_narrative,
                                                                          sign_B,
                                                                          zero_irf,
                                                                          max_tries)
      self$prior                   = specify_prior_bsvarSIGN$new(data, N, p, stationary)
      self$starting_values         = bsvars::specify_starting_values_bsvar$new(N, self$p, d)
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
#' the last MCMC draw for the Bayesian Structural VAR model with sign and narrative restrictions.
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
    
    #' @field posterior a list containing Bayesian estimation output including:
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

