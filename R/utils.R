
#' Importance sampling
#'
#' @description
#' Given posterior draws with importance weights, sample with replacement
#'
#' @param posterior posterior estimation outcome
#'
#' @export
importance_sampling <- function(posterior) {
  
  w <- posterior$posterior$w
  posterior$posterior <- posterior$posterior[-4] # remove weights
  
  if (posterior$last_draw$identification$sign_narrative[1, 1] == 0) {
    return(posterior)
  }
  
  indices <- sample(1:nrow(w), nrow(w), replace = TRUE, prob = w)
  
  posterior$posterior$A     <- posterior$posterior$A[, , indices]
  posterior$posterior$B     <- posterior$posterior$B[, , indices]
  posterior$posterior$hyper <- posterior$posterior$hyper[, , indices]
  posterior$posterior$Q     <- posterior$posterior$Q[, , indices]
  
  return(posterior)
}


#' Compute impulse response functions
#'
#' @description
#' Given posterior draws compute impulse response functions with bsvars
#'
#' @param posterior posterior estimation outcome
#' @param horizon a positive integer number denoting the forecast horizon for the impulse responses computations.
#'
#' @export
compute_irf <- function(posterior, horizon = 24) {
  
  irf_list <-
    irf_cpp(
      posterior$posterior$B,
      posterior$posterior$A,
      horizon,
      posterior$last_draw$p
    )
  
  S <- length(irf_list)
  irf <- array(dim = c(dim(irf_list[1][[1]]), S))
  for (s in 1:S) {
    irf[, , , s] <- irf_list[s][[1]]
  }
  
  return(irf)
}
