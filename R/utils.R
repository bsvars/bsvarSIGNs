
# Given posterior draws with importance weights, sample with replacement
importance_sampling <- function(posterior) {
  
  w <- posterior$posterior$w
  posterior$posterior <- posterior$posterior[-1] # remove weights
  
  # if (posterior$last_draw$identification$sign_narrative[1, 1] == 0) {
  #   return(posterior)
  # }
  
  indices <- sample(1:nrow(w), nrow(w), replace = TRUE, prob = w)
  
  posterior$posterior$A      = posterior$posterior$A[, , indices]
  posterior$posterior$B      = posterior$posterior$B[, , indices]
  posterior$posterior$hyper  = posterior$posterior$hyper[, , indices]
  posterior$posterior$Q      = posterior$posterior$Q[, , indices]
  posterior$posterior$Sigma  = posterior$posterior$Sigma[, , indices]
  posterior$posterior$Theta0 = posterior$posterior$Theta0[, , indices]
  
  posterior$posterior$ess    = sum(w)^2/sum(w^2)
  cat("The effective sample size is:", round(posterior$posterior$ess, 0), "\n")
  
  return(posterior)
}
