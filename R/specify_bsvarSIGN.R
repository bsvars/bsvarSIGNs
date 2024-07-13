
# construct Z_j matrices
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
  if (!(all(A[,3] %in% c(1:N, NA)))) {
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
#' # a prior for 3-variable example with one lag 
#' data(oil)
#' prior = specify_prior_bsvarSIGN$new(oil, p = 1, stationary = rep(TRUE, 3))
#' prior$B                                        # show autoregressive prior mean
#' 
#' @export
specify_prior_bsvarSIGN = R6::R6Class(
  "PriorBSVARSIGN",
  
  public = list(
    #' @field hyper a \code{(N+3)xS} matrix of hyper-parameters \eqn{\mu, \delta, \lambda, \psi}.
    hyper      = matrix(),
    
    #' @field B a \code{KxN} normal prior mean matrix for the autoregressive 
    #' parameters.
    B          = matrix(),
    
    #' @field V a \code{KxK} matrix determining  the normal prior column-specific 
    #' covariance for the autoregressive parameters.
    V          = matrix(),
    
    #' @field Vp a \code{px1} vector of lag shrinkage level.
    Vp         = matrix(),
    
    #' @field Vd a \code{(d+1)x1} vector of (large) variances for constants.
    Vd         = matrix(),
    
    #' @field S an \code{NxN} matrix determining the inverted-Wishart prior scale 
    #' of error terms covariance matrix.
    S          = matrix(),
    
    #' @field nu a positive scalar greater than \code{N+1} - the shape of the 
    #' inverted-Wishart prior for error terms covariance matrix.
    nu          = NA,
    
    #' @field data an \code{TxN} matrix of observations.
    data        = matrix(),
    
    #' @field Y an \code{TxN} matrix of dependent variables.
    Y           = matrix(),
    
    #' @field X an \code{TxK} matrix of independent variables.
    X           = matrix(),
    
    #' @field Ysoc an \code{NxN} matrix with the sum-of-coefficients dummy observations.
    Ysoc        = matrix(),
    
    #' @field Xsoc an \code{NxK} matrix with the sum-of-coefficients dummy observations.
    Xsoc        = matrix(),
    
    #' @field Ysur an \code{NxN} matrix with the single-unit-root dummy observations.
    Ysur        = matrix(),
    
    #' @field Xsur an \code{NxK} matrix with the single-unit-root dummy observations.
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
    #' # a prior for 3-variable example with one lag and stationary data
    #' data(oil)
    #' prior = specify_prior_bsvarSIGN$new(oil, p = 1, stationary = rep(TRUE, 3))
    #' prior$B # show autoregressive prior mean
    #' 
    initialize = function(data, p, exogenous = NULL, stationary = rep(FALSE, ncol(data))) {
      
      stopifnot("Argument p must be a positive integer number." = p > 0 & p %% 1 == 0)
      
      data    = bsvars::specify_data_matrices$new(data, p, exogenous)
      Y       = t(data$Y)
      X       = t(data$X)
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
      
      Vp      = (1:p)^-2
      Vd      = rep(100, 1 + d)
      V       = diag(c(kronecker(Vp, rep(1, N)), Vd))
      
      s2.ols  = rep(NA, N)
      for (n in 1:N) {
        y     = as.matrix(Y[(p + 5 + 1):T, n])
        x     = matrix(1, T - p - 5, 1)
        for (i in 1:(p + 5)) {
          x   = cbind(x, Y[(p + 5 + 1):T - i, n])
        }
        s2.ols[n] = sum(((diag(T - p - 5) - x %*% solve(t(x) %*% x) %*% t(x)) %*% y)^2) / 
          (T - p - 5)
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
      
      self$Y             = Y
      self$X             = X
      self$Vp            = Vp
      self$Vd            = Vd
      
      self$hyper         = hyper
      self$B             = B
      self$V             = V
      self$S             = diag(N)
      self$nu            = N + 2
      self$Ysoc          = Ysoc
      self$Xsoc          = Xsoc
      self$Ysur          = Ysur
      self$Xsur          = Xsur
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
    #' # a prior for 3-variable example with four lags
    #' prior = specify_prior_bsvar$new(N = 3, p = 4)
    #' prior$get_prior() # show the prior as list
    #' 
    get_prior = function(){
      list(
        hyper        = self$hyper,
        B            = self$B,
        V            = self$V,
        Vp           = self$Vp,
        Vd           = self$Vd,
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
    #' @param mu whether to estimate the hyper-parameter in the sum-of-coefficients dummy prior.
    #' @param delta whether to estimate the hyper-parameter in the single-unit-root dummy prior.
    #' @param lambda whether to estimate the hyper-parameter of shrinkage in the Minnesota prior.
    #' @param psi whether to estimate the hyper-parameter of variances in the Minnesota prior.
    #' @param S number of draws.
    #' @param start starting point of the adaptive Metropolis algorithm.
    #' @param burn number of burn-in draws.
    #' 
    #' @examples 
    #' # a prior for 3-variable example with four lags
    #' data(oil)
    #' prior = specify_prior_bsvarSIGN$new(oil, p = 1, stationary = rep(TRUE, 3))
    #' prior$estimate_hyper(S = 5)
    #' 
    estimate_hyper = function(mu = FALSE, delta = FALSE, lambda = TRUE, psi = FALSE,
                              S = 10000, start = 1000, burn = 5000){
      
      model = c(mu, delta, lambda, psi)
      
      if (all(!model)) {
        stop("At least one of the hyper-parameters must be estimated.")
      }
      
      hyper  = matrix(self$hyper[, ncol(self$hyper)])
      init   = narrow_hyper(model, hyper)
      prior  = self$get_prior()
      
      if (!mu) {
        prior$Ysoc = matrix(0, 0, ncol(self$Ysoc))
        prior$Xsoc = matrix(0, 0, ncol(self$Xsoc))
      }
      if (!delta) {
        prior$Ysur = matrix(0, 0, ncol(self$Ysur))
        prior$Xsur = matrix(0, 0, ncol(self$Xsur))
      }
      
      result = stats::optim(
        init,
        \(x) -log_posterior_hyper(extend_hyper(hyper, model, matrix(x)), 
                                  model, self$Y, self$X, prior),
        method  = 'L-BFGS-B',
        lower   = rep(0, length(init)),
        upper   = init * 100,
        hessian = TRUE
        )

      mode   = extend_hyper(hyper, model, matrix(result$par))
      W      = result$hessian

      if (length(init) == 1){
        W = 1 / W
      } else {
        e = eigen(W)
        W = e$vectors %*% diag(as.vector(1 / abs(e$values))) %*% t(e$vectors)
      }
      
      self$hyper = sample_hyper(S, start, mode, model, self$Y, self$X, W, prior)
      self$hyper = self$hyper[, -(1:burn)]
    } # END estimate_hyper
    
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
    #' Column 3 (var_i): an integer in 1:N (or NA when type = 0), index of the restricted variable. \cr
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
    #' Column 3 (var_i): an integer in 1:N (or NA when type = 0), index of the restricted variable. \cr
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
    #' Column 3 (var_i): an integer in 1:N (or NA when type = 0), index of the restricted variable. \cr
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
      self$data_matrices$Y         = t(self$data_matrices$Y)
      self$data_matrices$X         = t(self$data_matrices$X)
      self$identification          = specify_identification_bsvarSIGN$new(N,
                                                                          sign_irf,
                                                                          sign_narrative,
                                                                          sign_B,
                                                                          zero_irf,
                                                                          max_tries)
      self$prior                   = specify_prior_bsvarSIGN$new(data, p, exogenous, stationary)
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

