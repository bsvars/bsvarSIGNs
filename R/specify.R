
#' vector specifying one narrative restriction
#'
#' @description
#' The class narrative specifies a single narrative restriction.
#' 
#' @examples
#' # a prior for 5-variable example with one lag 
#' specify_narrative(start = 166, periods = 1, type = "S", sign = 1, shock = 1, var = 6)
#' 
#' @export
specify_narrative = function(start, periods = 1, type = "S", sign = 1, shock = 1, var = NA) {
  
  if (start %% 1 != 0 || start <= 0) {
    stop("start must be a positive integer")
  }
  if (periods %% 1 != 0 || periods <= 0) {
    stop("periods must be a positive integer")
  }
  if (!(type %in% c("S", "A", "B"))) {
    stop("type must be one of 'S', 'A', 'B'")
  }
  if (!(sign %in% c(-1, 1))) {
    stop("sign must be one of -1, 1")
  }
  if (shock %% 1 != 0 || shock <= 0) {
    stop("shock must be a positive integer")
  }
  if (!is.na(var)) {
    if (var %% 1 != 0 || var <= 0) {
      stop("var must be a positive integer")
    }
  }
  
  narrative = list(
    start   = start,
    periods = periods,
    type    = type,
    sign    = sign,
    shock   = shock,
    var     = var
  )
  class(narrative) = "narrative"
  narrative
}

# construct Z_j matrices
get_Z = function(sign_irf) {
  h = dim(sign_irf)[3]
  if (h >= 2) {
    test = sign_irf[, , 1:h]
    if (any(test[!is.na(test)] == 0)) {
      stop("Zero restrictions are not allowed for horizons >= 1")
    }
  }
  
  zero_irf = sign_irf[, , 1] == 0
  zero_irf[is.na(zero_irf)] = 0
  
  if (sum(zero_irf) == 0) {
    return(NULL)
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
  if (!(all(A %in% c(-1, 0, 1, NA)))) {
    stop("Sign restriction matrix has entries that are not in {-1, 0, 1, NA}.")
  }
}

# verify all restrictions
verify_all = function(N, sign_irf, sign_narrative, sign_relation) {
  verify_traditional(N, sign_relation)
  if (any(sign_relation[!is.na(sign_relation)] == 0)) {
    stop("Zero restrictions are not allowed for sign_relation")
  }
  
  if (!is.list(sign_narrative)) {
    stop("sign_narrative must be a list")
  } else if (length(sign_narrative) > 0) {
    for (i in 1:length(sign_narrative)) {
      if (!inherits(sign_narrative[[i]], "narrative")) {
        stop("Each element of sign_narrative must be of class 'narrative'")
      }
    }
  }
  
  for (h in 1:dim(sign_irf)[3]) {
    verify_traditional(N, sign_irf[,,h])
  }
}


# Minnesota prior of Normal-Inverse-Wishart form
# niw_prior = function(Y,
#                      p,
#                      non_stationary,
#                      lambda = 0.2) {
#   T = nrow(Y)
#   N = ncol(Y)
#   K = 1 + N * p
#   
#   B           = matrix(0, K, N)
#   B[1:N, 1:N] = diag(non_stationary)
#   
#   sigma2                  = sapply(1:N, \(i) summary(lm(Y[2:T, i] ~ Y[1:(T - 1), i]))$sigma^2)
#   V                       = matrix(0, K, K)
#   V[K, K]                 = 1e+6
#   V[1:(K - 1), 1:(K - 1)] = diag(lambda^2 * kronecker((1:p)^-2, sigma2^-1))
#   
#   S  = diag(sigma2)
#   nu = N + 2
#   
#   list(B = B, V = V, S = S, nu = nu)
# }

gamma_scale = function(mode, variance) {
  (2 * mode * variance) / (mode^2 + sqrt(mode^4 + 4 * variance * mode^2))
}

gamma_shape = function(mode, variance) {
  (mode^2 + sqrt(mode^4 + 4 * (mode^2) * variance) + 2 * variance) / (2 * variance)
}

igamma_scale = function(mode, variance) {
  (0.5 * mode * (sqrt(variance * (variance + 24 * mode^2)) - 5 * variance)) / (mode^2 - variance)
}

igamma_shape = function(mode, variance) {
  0.5 * (sqrt(variance * (variance + 24 * mode^2)) - 2 * mode^2 - 3 * variance ) / (mode^2 - variance)
}

#' R6 Class Representing PriorBSVAR
#'
#' @description
#' The class PriorBSVARSIGN presents a prior specification for the homoskedastic bsvar model.
#' 
#' @examples
#' # a prior for 5-variable example with one lag 
#' data(optimism)
#' prior = specify_prior_bsvarSIGN$new(optimism, p = 1)
#' prior$A  # show autoregressive prior mean
#' 
#' @export
specify_prior_bsvarSIGN = R6::R6Class(
  "PriorBSVARSIGN",
  
  public = list(
    #' @field p a positive integer - the number of lags.
    p           = 1,
    
    #' @field hyper a \code{(N+3)xS} matrix of hyper-parameters \eqn{\mu, \delta, \lambda, \psi}.
    hyper      = matrix(),
    
    #' @field A a \code{NxK} normal prior mean matrix for the autoregressive 
    #' parameters.
    A          = matrix(),
    
    #' @field V a \code{KxK} matrix determining  the normal prior column-specific 
    #' covariance for the autoregressive parameters.
    V          = matrix(),
    
    #' @field S an \code{NxN} matrix determining the inverted-Wishart prior scale 
    #' of error terms covariance matrix.
    S          = matrix(),
    
    #' @field nu a positive scalar greater than \code{N+1} - the shape of the 
    #' inverted-Wishart prior for error terms covariance matrix.
    nu          = NA,
    
    #' @field data an \code{TxN} matrix of observations.
    data        = matrix(),
    
    #' @field Y an \code{NxT} matrix of dependent variables.
    Y           = matrix(),
    
    #' @field X an \code{KxT} matrix of independent variables.
    X           = matrix(),
    
    #' @field Ysoc an \code{NxN} matrix with the sum-of-coefficients dummy observations.
    Ysoc        = matrix(),
    
    #' @field Xsoc an \code{KxN} matrix with the sum-of-coefficients dummy observations.
    Xsoc        = matrix(),
    
    #' @field Ysur an \code{NxN} matrix with the single-unit-root dummy observations.
    Ysur        = matrix(),
    
    #' @field Xsur an \code{KxN} matrix with the single-unit-root dummy observations.
    Xsur        = matrix(),
    
    #' @field mu.scale a positive scalar - the shape of the gamma prior for \eqn{\mu}.
    mu.scale    = NA,
    
    #' @field mu.shape a positive scalar - the shape of the gamma prior for \eqn{\mu}.
    mu.shape    = NA,
    
    #' @field delta.scale a positive scalar - the shape of the gamma prior for \eqn{\delta}.
    delta.scale = NA,

    #' @field delta.shape a positive scalar - the shape of the gamma prior for \eqn{\delta}.
    delta.shape = NA,
    
    #' @field lambda.scale a positive scalar - the shape of the gamma prior for \eqn{\lambda}.
    lambda.scale = NA,
    
    #' @field lambda.shape a positive scalar - the shape of the gamma prior for \eqn{\lambda}.
    lambda.shape = NA,
    
    #' @field psi.scale a positive scalar - the shape of the inverted gamma prior for \eqn{\psi}.
    psi.scale   = NA,
    
    #' @field psi.shape a positive scalar - the shape of the inverted gamma prior for \eqn{\psi}.
    psi.shape   = NA,
    
    #' @description
    #' Create a new prior specification PriorBSVAR.
    #' @param data the \code{TxN} data matrix of observations.
    #' @param p a positive integer - the autoregressive lag order of the SVAR model.
    #' @param exogenous a \code{Txd} matrix of exogenous variables.
    #' @param stationary an \code{N} logical vector - its element set to \code{FALSE} sets 
    #' the prior mean for the autoregressive parameters of the \code{N}th equation to the white noise process, 
    #' otherwise to random walk.
    #' @return A new prior specification PriorBSVARSIGN.
    #' @examples 
    #' # a prior for 5-variable example with one lag and stationary data
    #' data(optimism)
    #' prior = specify_prior_bsvarSIGN$new(optimism, p = 1)
    #' prior$B # show autoregressive prior mean
    #' 
    initialize = function(data, p, exogenous = NULL, stationary = rep(FALSE,  ncol(data))) {
      
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      
      data_m  = bsvars::specify_data_matrices$new(data, p, exogenous)
      Y       = t(data_m$Y)
      X       = t(data_m$X)
      N       = ncol(Y)
      
      stopifnot("Argument stationary must be a logical vector of length equal to the number of columns in data." = length(stationary) == N & is.logical(stationary))
      
      d       = 0
      if (!is.null(exogenous)) {
        d     = ncol(exogenous)
      }
      T       = nrow(Y)
      K       = N * p + 1 + d
      
      B       = matrix(0, K, N)
      B[1:N,] = diag(!stationary)
      
      V       = diag(c(kronecker((1:p)^-2, rep(1, N)), rep(100, 1 + d)))
      
      s2.ols  = rep(NA, N)
      for (n in 1:N) {
        y     = as.matrix(Y[(p + 5 + 1):T, n])
        x     = matrix(1, T - p - 5, 1)
        for (i in 1:(p + 5)) {
          x   = cbind(x, Y[(p + 5 + 1):T - i, n])
        }
        s2.ols[n] = sum(((diag(T - p - 5) - x %*% solve(t(x) %*% x) %*% t(x)) %*% y)^2) / (T - p - 5)
      }
      
      hyper              = matrix(NA, N + 3, 1)
      hyper[1:3]         = c(1, 1, 0.2)
      hyper[4:(N + 3),]  = s2.ols
      
      scale   = gamma_scale(1, 1)
      shape   = gamma_shape(1, 1)
      
      ybar    = colMeans(matrix(Y[1:p,], ncol = N))
      Ysoc    = diag(ybar)
      Ysur    = t(ybar)
      Xsoc    = cbind(kronecker(t(rep(1, p)), Ysoc), matrix(0, N, d + 1))
      Xsur    = cbind(kronecker(t(rep(1, p)), Ysur), 1, matrix(0, 1, d))
      
      # Ystar   = rbind(diag(ybar), ybar)
      # Xstar   = Ystar
      # if (p > 1) {
      #   for (i in 2:p) {
      #     Xstar = cbind(Xstar, Ystar)
      #   }
      # }
      # Xstar   = cbind(Xstar, c(rep(0, N), 1), matrix(0, N + 1, d))
      
      self$p             = p
      self$hyper         = hyper
      self$A             = t(B)
      self$V             = V
      self$S             = diag(N)
      self$nu            = N + 2
      self$Y             = t(Y)
      self$X             = t(X)
      self$Ysoc          = t(Ysoc)
      self$Xsoc          = t(Xsoc)
      self$Ysur          = t(Ysur)
      self$Xsur          = t(Xsur)
      self$mu.scale      = scale
      self$mu.shape      = shape
      self$delta.scale   = scale
      self$delta.shape   = shape
      self$lambda.scale  = gamma_scale(0.2, 0.4)
      self$lambda.shape  = gamma_shape(0.2, 0.4)
      self$psi.scale     = igamma_scale(0.02^2, 0.02^2)
      self$psi.shape     = igamma_shape(0.02^2, 0.02^2)
    }, # END initialize
    
    #' @description
    #' Returns the elements of the prior specification PriorBSVAR as a \code{list}.
    #' 
    #' @examples 
    #' # a prior for 5-variable example with four lags
    #' prior = specify_prior_bsvar$new(N = 5, p = 4)
    #' prior$get_prior() # show the prior as list
    #' 
    get_prior = function(){
      list(
        p            = self$p,
        hyper        = self$hyper,
        A            = self$A,
        V            = self$V,
        S            = self$S,
        nu           = self$nu,
        Ysoc         = self$Ysoc,
        Xsoc         = self$Xsoc,
        Ysur         = self$Ysur,
        Xsur         = self$Xsur,
        mu.scale     = self$mu.scale,
        mu.shape     = self$mu.shape,
        delta.scale  = self$delta.scale,
        delta.shape  = self$delta.shape,
        lambda.scale = self$lambda.scale,
        lambda.shape = self$lambda.shape,
        psi.scale    = self$psi.scale,
        psi.shape    = self$psi.shape
      )
    }, # END get_prior
    
    #' @description
    #' Estimates hyper-parameters with adaptive Metropolis algorithm.
    #' 
    #' @param mu whether to estimate the hyper-parameter in the 
    #' sum-of-coefficients dummy prior.
    #' @param delta whether to estimate the hyper-parameter in the 
    #' single-unit-root dummy prior.
    #' @param lambda whether to estimate the hyper-parameter of the 
    #' shrinkage in the Minnesota prior.
    #' @param psi whether to estimate the hyper-parameter of the 
    #' variances in the Minnesota prior.
    #' @param S number of MCMC draws.
    #' @param burn_in number of burn-in draws.
    #' 
    #' @examples 
    #' # a prior for 5-variable example with four lags
    #' data(optimism)
    #' prior = specify_prior_bsvarSIGN$new(optimism, p = 1)
    #' prior$estimate_hyper(S = 5)
    #' 
    estimate_hyper = function(
      S = 10000, burn_in = S / 2,
      mu = FALSE, delta = FALSE, lambda = TRUE, psi = FALSE
      ) {
      
      model = c(mu, delta, lambda, psi)
      
      if (all(!model)) {
        stop("At least one of the hyper-parameters must be estimated.")
      }
      
      hyper  = matrix(self$hyper[, ncol(self$hyper)])
      init   = narrow_hyper(model, hyper)
      prior  = self$get_prior()
      
      prior$B    = t(prior$A)
      prior$Ysoc = t(prior$Ysoc)
      prior$Xsoc = t(prior$Xsoc)
      prior$Ysur = t(prior$Ysur)
      prior$Xsur = t(prior$Xsur)
      
      result = stats::optim(
        init,
        \(x) -log_posterior_hyper(extend_hyper(hyper, model, matrix(x)), 
                                  model, t(self$Y), t(self$X), prior),
        method  = 'L-BFGS-B',
        lower   = rep(0, length(init)),
        upper   = init * 100,
        hessian = TRUE
        )

      mode       = extend_hyper(hyper, model, matrix(result$par))
      variance   = result$hessian

      if (length(init) == 1){
        variance = 1 / variance
      } else {
        e        = eigen(variance)
        variance = e$vectors %*% diag(as.vector(1 / abs(e$values))) %*% t(e$vectors)
      }
      
      self$hyper = sample_hyper(S, burn_in, mode, model, 
                                t(self$Y), t(self$X), variance, prior)
      self$hyper = self$hyper[, -(1:burn_in)]
    } # END estimate_hyper
    
  ) # END public
) # END specify_prior_bsvarSIGN


#' R6 Class Representing IdentificationBSVARSIGN
#'
#' @description
#' The class IdentificationBSVARSIGN presents the identifying restrictions for the Bayesian Structural VAR models with sign and narrative restrictions.
#'
#' @examples 
#' specify_identification_bsvarSIGN$new(N = 5) # recursive specification for a 5-variable system
#' 
#' # an identification pattern with narrative sign restrictions
#' sign_irf = matrix(c(0, 1, rep(NA, 23)), 5, 5)
#' specify_identification_bsvarSIGN$new(N = 5, sign_irf = sign_irf) 
#'
#' @export
specify_identification_bsvarSIGN = R6::R6Class(
  "IdentificationBSVARSIGN",
  
  public = list(
    
    #' @field VB a list of \code{N} matrices determining the unrestricted elements of matrix \eqn{B}.
    VB       = list(),
    #' @field sign_irf a \code{NxNxH} array of sign restrictions on the impulse response functions.
    sign_irf = array(),
    #' @field sign_narrative a \code{ANYx6} matrix of narrative sign restrictions.
    sign_narrative  = matrix(),
    #' @field sign_relation a \code{NxN} matrix of sign restrictions on contemporaneous relations.
    sign_relation   = matrix(),
    #' @field max_tries a positive integer with the maximum number of iterations 
    #' for finding a rotation matrix \eqn{Q} that would satisfy sign restrictions.
    max_tries = 1,
    
    #' @description
    #' Create new identifying restrictions IdentificationBSVARSIGN.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param sign_irf a \code{NxNxH} array - sign and zero restrictions 
    #' on the impulse response functions, ±1 for positive/negative sign restriction
    #' 0 for zero restrictions and NA for no restrictions,
    #' the \code{h}-th slice \code{NxN} matrix contains the
    #' restrictions on the \code{h-1} horizon.
    #' @param sign_narrative a list of objects of class "narrative" - narrative sign restrictions.
    #' @param sign_relation a \code{NxN} matrix with entries ±1 or NA - sign restrictions on the
    #' contemporaneous relations \code{B} between reduced-form errors \code{E} and
    #' structural shocks \code{U} where \code{BE=U}.
    #' @param max_tries a positive integer with the maximum number of iterations
    #' for finding a rotation matrix \eqn{Q} that would satisfy sign restrictions.
    #' @return Identifying restrictions IdentificationBSVARSIGN.
    initialize = function(N, sign_irf, sign_narrative, sign_relation, max_tries = 1) {
        
      missing_all   = TRUE
      if (missing(sign_irf)) {
        sign_irf = array(rep(NA, N^2), dim = c(N, N, 1))
      } else {
        missing_all = FALSE
      }
      if (missing(sign_narrative)) {
        sign_narrative = list()
      } else {
        missing_all = FALSE
      }
      if (missing(sign_relation)) {
        sign_relation = matrix(rep(NA, N^2), ncol = N, nrow = N)
        if (missing_all) {
          diag(sign_relation) = 1
        }
      }
      
      if (is.matrix(sign_irf)) {
        sign_irf = array(sign_irf, dim = c(dim(sign_irf), 1))
      }
      verify_all(N, sign_irf, sign_narrative, sign_relation)
      
      B     = matrix(FALSE, N, N)
      B[lower.tri(B, diag = TRUE)] = TRUE
      
      self$VB = vector("list", N)
      for (n in 1:N) {
        self$VB[[n]] = matrix(diag(N)[B[n,],], ncol = N)
      }
      
      self$sign_irf       = sign_irf
      self$sign_narrative = sign_narrative
      self$sign_relation  = sign_relation
      self$max_tries      = max_tries
    }, # END initialize
    
    #' @description
    #' Returns the elements of the identification pattern IdentificationBSVARSIGN as a \code{list}.
    #'
    get_identification = function() {
      list(
        VB             = as.list(self$VB),
        sign_irf       = as.array(self$sign_irf),
        sign_narrative = self$sign_narrative,
        sign_relation  = as.matrix(self$sign_relation),
        max_tries      = self$max_tries
        )
    }, # END get_identification
    
    #' @description
    #' Set new starting values StartingValuesBSVARSIGN.
    #' @param N a positive integer - the number of dependent variables in the model.
    #' @param sign_irf a \code{NxNxH} array - sign and zero restrictions 
    #' on the impulse response functions, ±1 for positive/negative sign restriction
    #' 0 for zero restrictions and NA for no restrictions,
    #' the \code{h}-th slice \code{NxN} matrix contains the
    #' restrictions on the \code{h-1} horizon.
    #' @param sign_narrative a list of objects of class "narrative" - narrative sign restrictions.
    #' @param sign_relation a \code{NxN} matrix with entries ±1 or NA - sign restrictions on the
    #' contemporaneous relations \code{B} between reduced-form errors \code{E} and
    #' structural shocks \code{U} where \code{BE=U}.
    #' @param max_tries a positive integer with the maximum number of iterations
    #' for finding a rotation matrix \eqn{Q} that would satisfy sign restrictions.
    set_identification = function(N, sign_irf, sign_narrative, sign_relation) {
      B     = matrix(FALSE, N, N)
      B[lower.tri(B, diag = TRUE)] = TRUE
      
      self$VB          <- vector("list", N)
      for (n in 1:N) {
        self$VB[[n]]   <- matrix(diag(N)[B[n,],], ncol = N)
      }
      
      missing_all   = TRUE
      if (missing(sign_irf)) {
        sign_irf = array(rep(NA, N^2), dim = c(N, N, 1))
      } else {
        missing_all = FALSE
      }
      if (missing(sign_narrative)) {
        sign_narrative = list()
      } else {
        missing_all = FALSE
      }
      if (missing(sign_relation)) {
        sign_relation = matrix(rep(NA, N^2), ncol = N, nrow = N)
        if (missing_all) {
          diag(sign_relation) = 1
        }
      }
      
      if (is.matrix(sign_irf)) {
        sign_irf = array(sign_irf, dim = c(dim(sign_irf), 1))
      }
      verify_all(N, sign_irf, sign_narrative, sign_relation)
      
      self$sign_irf       = sign_irf
      self$sign_narrative = sign_narrative
      self$sign_relation  = sign_relation
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
#' data(optimism)
#' specification = specify_bsvarSIGN$new(
#'    data = optimism,
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
    #' @param sign_irf a \code{NxNxH} array - sign and zero restrictions 
    #' on the impulse response functions, ±1 for positive/negative sign restriction
    #' 0 for zero restrictions and NA for no restrictions,
    #' the \code{h}-th slice \code{NxN} matrix contains the
    #' restrictions on the \code{h-1} horizon.
    #' @param sign_narrative a list of objects of class "narrative" - narrative sign restrictions.
    #' @param sign_relation a \code{NxN} matrix with entries ±1 or NA - sign restrictions on the
    #' contemporaneous relations \code{B} between reduced-form errors \code{E} and
    #' structural shocks \code{U} where \code{BE=U}.
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
    sign_relation,
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
        sign_irf = array(rep(NA, N^2), dim = c(N, N, 1))
      } else {
        missing_all = FALSE
      }
      if (missing(sign_narrative)) {
        sign_narrative = list()
      } else {
        missing_all = FALSE
      }
      if (missing(sign_relation)) {
        sign_relation = matrix(rep(NA, N^2), ncol = N, nrow = N)
        if (missing_all) {
          diag(sign_relation) = 1
        }
      }
      
      if (is.matrix(sign_irf)) {
        sign_irf = array(sign_irf, dim = c(dim(sign_irf), 1))
      }
      verify_all(N, sign_irf, sign_narrative, sign_relation)
      
      B                            = matrix(FALSE, N, N)
      B[lower.tri(B, diag = TRUE)] = TRUE
      
      self$data_matrices           = bsvars::specify_data_matrices$new(data, p, exogenous)
      self$identification          = specify_identification_bsvarSIGN$new(N,
                                                                          sign_irf,
                                                                          sign_narrative,
                                                                          sign_relation,
                                                                          max_tries)
      self$prior                   = specify_prior_bsvarSIGN$new(data, p, exogenous,
                                                                 stationary)
      self$starting_values         = bsvars::specify_starting_values_bsvar$new(N, self$p, d)
    }, # END initialize
    
    #' @description
    #' Returns the data matrices as the DataMatricesBSVAR object.
    #'
    #' @examples 
    #' data(optimism)
    #' spec = specify_bsvarSIGN$new(
    #'    data = optimism,
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
    #' data(optimism)
    #' spec = specify_bsvarSIGN$new(
    #'    data = optimism,
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
    #' data(optimism)
    #' spec = specify_bsvarSIGN$new(
    #'    data = optimism,
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
    #' data(optimism)
    #' spec = specify_bsvarSIGN$new(
    #'    data = optimism,
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
#' data(optimism)
#' specification  = specify_bsvarSIGN$new(optimism, p = 4)
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
    #' data(optimism)
    #' specification  = specify_bsvarSIGN$new(optimism)
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
    #' data(optimism)
    #' specification  = specify_bsvarSIGN$new(optimism)
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

