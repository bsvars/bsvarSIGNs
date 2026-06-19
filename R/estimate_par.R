#' @title Parallel Bayesian estimation of a Structural Vector Autoregression
#' with traditional and narrative sign restrictions
#'
#' @description Parallelised version of \code{estimate} using the \code{parallel} package.
#' Estimates Bayesian Structural Vector Autoregression model
#' using independent sampling. 
#' 
#' @param specification an object of class BSVARSIGN generated using the \code{specify_bsvarSIGN$new()} function.
#' @param S a positive integer, the number of posterior draws to be generated
#' @param thin a positive integer, specifying the frequency of MCMC output thinning
#' @param show_progress a logical value, if \code{TRUE} the estimation progress bar is visible
#' @param ncores a positive integer, the number of cores to be used for parallel execution
#' 
#' @return An object of class \code{PosteriorBSVARSIGN}
#' 
#' @export
estimate_par = function(specification, S, thin = 1, show_progress = TRUE, ncores = parallel::detectCores() - 1) {
  UseMethod("estimate_par", specification)
}

#' @export
estimate_par.BSVARSIGN = function(specification, S, thin = 1, show_progress = TRUE, ncores = parallel::detectCores() - 1) {
  
  # get the inputs to estimation
  prior               = specification$prior$get_prior()
  identification      = specification$identification$get_identification()
  max_tries           = identification$max_tries
  max_tries           = ifelse(max_tries == Inf, 0, max_tries)
  data_matrices       = specification$data_matrices$get_data_matrices()
  p                   = specification$p
  
  prior$B             = t(prior$A)
  prior$Ysoc          = t(prior$Ysoc)
  prior$Xsoc          = t(prior$Xsoc)
  prior$Ysur          = t(prior$Ysur)
  prior$Xsur          = t(prior$Xsur)
  Y                   = t(data_matrices$Y)
  X                   = t(data_matrices$X)
  
  Z                   = get_Z(identification$sign_irf)
  sign                = identification$sign_irf
  sign[is.na(sign)]   = 0
  
  n_narratives        = length(identification$sign_narrative)
  get_type            = list("S" = 1, "A" = 2, "B" = 3)
  if (n_narratives > 0) {
    narrative         = matrix(NA, n_narratives, 6)
    for (i in 1:n_narratives) {
      narrative_list  = identification$sign_narrative[[i]]
      narrative[i, 1] = get_type[[narrative_list$type]]
      narrative[i, 2] = narrative_list$sign
      narrative[i, 3] = narrative_list$var
      narrative[i, 4] = narrative_list$shock
      narrative[i, 5] = narrative_list$start - p
      narrative[i, 6] = narrative_list$periods - 1
    }
  } else {
    narrative         = t(c(0, 1, 1, 1, 1, 1))
  }
  struc               = identification$sign_structural
  struc[is.na(struc)] = 0
  Nf                  = specification$num_foreign_vars
  
  T_obs = nrow(Y)
  N     = ncol(Y)
  K     = ncol(X)
  
  if (show_progress) {
    message("**************************************************|")
    message(" bsvarSIGNs: Bayesian Structural VAR with sign,   |")
    message("             zero and narrative restrictions      |")
    message("**************************************************|")
    message(" Parallel execution on ", ncores, " cores")
    message(" Progress of simulation for ", S, " independent draws")
    message(" Press Esc to interrupt the computations")
    message("**************************************************|")
  }
  
  # Ensure integer S
  S = as.integer(S)
  
  # Determine OS for parallel method
  is_windows = .Platform$OS.type == "windows"
  
  if (is_windows) {
    cl = parallel::makeCluster(ncores, outfile = "")
    on.exit(parallel::stopCluster(cl))
    
    # Export necessary variables
    parallel::clusterExport(cl, varlist = c("p", "Y", "X", "sign", "narrative", "struc", "Z", "Nf", "prior", "max_tries"), envir = environment())
    
    # Set RNG for reproducibility
    parallel::clusterSetRNGStream(cl)
    
    results = parallel::parLapply(cl, 1:S, function(s) {
      if (show_progress && (s %% max(1, round(S/50)) == 0)) cat(".")
      .Call(`_bsvarSIGNs_bsvar_sign_par_cpp`, p, Y, X, sign, narrative, struc, Z, Nf, prior, max_tries)
    })
    
  } else {
    # Set RNG for reproducibility across forks
    if (show_progress) message("Running with mclapply...")
    
    results = parallel::mclapply(1:S, function(s) {
      if (show_progress && (s %% max(1, round(S/50)) == 0)) cat(".")
      .Call(`_bsvarSIGNs_bsvar_sign_par_cpp`, p, Y, X, sign, narrative, struc, Z, Nf, prior, max_tries)
    }, mc.cores = ncores, mc.set.seed = TRUE)
  }
  
  if (show_progress) cat("\n")
  
  # Check if results have errors
  if (inherits(results[[1]], "try-error")) {
    stop("Error in parallel execution: ", results[[1]])
  }
  
  # Pre-allocate arrays
  posterior_w      = matrix(NA, nrow = S, ncol = 1)
  posterior_hyper  = matrix(NA, nrow = nrow(prior$hyper), ncol = S)
  posterior_A      = array(NA, dim = c(N, K, S))
  posterior_B      = array(NA, dim = c(N, N, S))
  posterior_Q      = array(NA, dim = c(N, N, S))
  posterior_Sigma  = array(NA, dim = c(N, N, S))
  posterior_Theta0 = array(NA, dim = c(N, N, S))
  posterior_shocks = array(NA, dim = c(N, T_obs, S))
  
  for (s in 1:S) {
    res = results[[s]]
    posterior_w[s, 1]     = res$w
    posterior_hyper[, s]  = res$hyper
    posterior_A[, , s]    = res$A
    posterior_B[, , s]    = res$B
    posterior_Q[, , s]    = res$Q
    posterior_Sigma[, , s] = res$Sigma
    posterior_Theta0[, , s]= res$Theta0
    posterior_shocks[, , s]= res$shocks
  }
  
  qqq = list(
    posterior = list(
      w      = posterior_w,
      hyper  = posterior_hyper,
      A      = posterior_A,
      B      = posterior_B,
      Q      = posterior_Q,
      Sigma  = posterior_Sigma,
      Theta0 = posterior_Theta0,
      shocks = posterior_shocks
    ),
    last_draw = results[[S]]
  )
  
  # specification$starting_values$set_starting_values(qqq$last_draw)
  output = specify_posterior_bsvarSIGN$new(specification, qqq$posterior)
  output = importance_sampling(output)
  
  return(output)
}
