
log_det = function(A) {
  determinant(A, logarithm = TRUE)$modulus
}

dlambda = function(lambda, p, Y, X) {
  
  N     = ncol(Y)
  T     = nrow(Y)
  prior = niw_prior(Y, p, rep(1, N), lambda)
  
  b     = prior$B
  Omega = prior$V
  Psi   = prior$S
  d     = prior$nu
  
  inv_Omega = solve(Omega)
  Bhat      = solve(t(X) %*% X + inv_Omega) %*% (t(X) %*% Y + inv_Omega %*% b)
  ehat      = Y - X %*% Bhat
  
  lprior    = dgamma(
    lambda,
    shape = (9 + sqrt(17)) / 8,
    scale = 8 / (5 * (1 + sqrt(17))),
    log = TRUE
  )
  
  llike     = -N / 2 * log_det(Omega)
  llike     = llike - N / 2 * log_det(t(X) %*% X + inv_Omega)
  A         = Psi + t(ehat) %*% ehat + t(Bhat - b) %*% inv_Omega %*% (Bhat - b)
  llike     = llike - (T + d) / 2 * log_det(A)
  
  return(lprior + llike)
}


sample_lambda = function(S, p, Y, X) {
  init   = 0.2
  result = optim(
    init,
    \(lambda) {
      nlogp = -dlambda(lambda, p, Y, X)
      if (!is.finite(nlogp)) {
        nlogp = 1e+6
      }
      return(nlogp)
    },
    method  = "L-BFGS-B",
    lower = init / 100,
    upper = init * 100,
    hessian = TRUE,
    control = list(trace = 1, maxit = 1e4)
  )
  
  lambda    = numeric(S)
  lambda[1] = result$par
  lambda[1] = 0.2
  dlambda   = result$value
  
  c         = 0.5
  W         = c * 1 / result$hessian
  
  for (s in 2:S) {
    new_lambda  = rnorm(1, lambda[s - 1], sqrt(W))
    new_dlambda = dlambda(new_lambda, p, Y, X)
    
    if (runif(1) < exp(new_dlambda - dlambda)) {
      lambda[s] = new_lambda
      dlambda   = new_dlambda
    } else {
      lambda[s] = lambda[s - 1]
    }
  }
  
  return(lambda)
}


library(R.matlab)
data = readMat('../data/DataSW.mat')
name = sapply(1:7, \(i) data$ShortDescr[[i]][[1]])
data = data$y
colnames(data) = name
data = data[, c(1, 2, 7)]

p = 4
d = bsvars::specify_data_matrices$new(data, p, NULL)
Y = t(d$Y)
X = t(d$X)

lambda = sample_lambda(1000, p, Y, X)

par(mfrow = c(2, 1))
plot(lambda, type = 'l')
hist(lambda, breaks = 100)





