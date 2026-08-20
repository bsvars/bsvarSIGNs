
# Resample posterior draws with unequal importance weights
importance_sampling <- function(posterior) {
  
  w <- posterior$posterior$w
  posterior$posterior$w <- NULL
  posterior$posterior$ess <- sum(w)^2/sum(w^2)
  
  if (all(w == 1)) {
    return(posterior)
  }
  
  S <- length(w)
  indices <- sample(seq_len(S), S, replace = TRUE, prob = w)
  
  posterior$posterior$A      = posterior$posterior$A[, , indices]
  posterior$posterior$B      = posterior$posterior$B[, , indices]
  posterior$posterior$hyper  = posterior$posterior$hyper[, indices]
  posterior$posterior$Q      = posterior$posterior$Q[, , indices]
  posterior$posterior$Sigma  = posterior$posterior$Sigma[, , indices]
  posterior$posterior$Theta0 = posterior$posterior$Theta0[, , indices]
  return(posterior)
}
