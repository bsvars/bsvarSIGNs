rm(list = ls())

niw = function(Y,
               p,
               non_stationary,
               lambda,
               psi) {
  T = nrow(Y)
  N = ncol(Y)
  K = 1 + N * p
  
  B           = matrix(0, K, N)
  B[1:N, 1:N] = diag(non_stationary)
  
  V                       = matrix(0, K, K)
  V[K, K]                 = 1e+6
  V[1:(K - 1), 1:(K - 1)] = diag(lambda^2 * kronecker((1:p)^-2, psi^-1))
  
  S  = diag(psi)
  nu = N + 2
  
  list(B = B, V = V, S = S, nu = nu)
}

log_det = function(A) {
  determinant(A, logarithm = TRUE)$modulus
}

dhyper = function(p, hyper, model, Y, X) {
  N     = ncol(Y)
  
  if (model[3]) {
    lambda = hyper[3]
  } else {
    lambda = 0.2
  }
  
  if (model[4]) {
    psi = hyper[4:length(hyper)]
  } else {
    psi = sapply(1:N, \(i) summary(lm(Y[2:T, i] ~ Y[1:(T - 1), i]))$sigma^2)
  }
  
  prior = niw(Y, p, rep(1, N), lambda, psi)
  
  extended = extend_dummy(p, hyper, model, Y, X)
  Ystar    = extended$Ystar
  Xstar    = extended$Xstar
  Yplus    = extended$Yplus
  Xplus    = extended$Xplus
  
  logml_plus = logml(prior, p, Yplus, Xplus)
  if (dim(Ystar)[1] > 0) {
    logml_star = logml(prior, p, Ystar, Xstar)
  } else {
    logml_star = 0
  }
  
  as.numeric(logml_plus - logml_star)
}

logml = function(prior, p, Y, X) {
  
  N     = ncol(Y)
  T     = nrow(Y)
  
  b     = prior$B
  Omega = prior$V
  Psi   = prior$S
  d     = prior$nu
  
  llike = -1e10
  tryCatch({
    inv_Omega = solve(Omega)
    Bhat      = solve(t(X) %*% X + inv_Omega, t(X) %*% Y + inv_Omega %*% b)
    ehat      = Y - X %*% Bhat
    
    llike     = - N * T / 2 * log(pi) 
    llike     = llike + log_mvgamma(N, (T + d) / 2) - log_mvgamma(N, d / 2)
    llike     = llike - N / 2 * log_det(Omega)
    llike     = llike + d / 2 * log_det(Psi)
    llike     = llike - N / 2 * log_det(t(X) %*% X + inv_Omega)
    A         = Psi + t(ehat) %*% ehat + t(Bhat - b) %*% inv_Omega %*% (Bhat - b)
    llike     = llike - (T + d) / 2 * log_det(A)
  }, error = function(e) {
  })
  
  if (!is.finite(llike)) {
    llike = -1e+10
  }
  
  return(llike)
}

set.seed(123)

dhyper = log_posterior_hyper

library(R.matlab)
data = readMat('../data/DataSW.mat')
name = sapply(1:7, \(i) data$ShortDescr[[i]][[1]])
data = data$y
colnames(data) = name
# data = data[, c(1, 2, 7)]

p = 4
d = bsvars::specify_data_matrices$new(data, p, NULL)
Y = t(d$Y)
X = t(d$X)
T = nrow(Y)
N = ncol(Y)

init   = c(1, 1, 0.2)
init   = c(init, rep(0.02^2, N))
model  = c(FALSE, FALSE, TRUE, TRUE)
# model[4] = TRUE

result = optim(narrow_hyper(model, matrix(init)),
               \(x) -dhyper(p, extend_hyper(init, model, matrix(x)), model, Y, X),
               method  = 'L-BFGS-B',
               control = list(trace = 1, maxit = 1e5),
               lower   = rep(0, length(init)),
               upper   = init * 100,
               hessian = TRUE
               )

h = result$par
d = dhyper(p, extend_hyper(init, model, matrix(h)), model, Y, X)
h = log(h)


c = 1
W = result$hessian
if (length(h) == 1){
  W = 1 / W
} else {
  e = eigen(W)
  W = e$vectors %*% diag(as.vector(1 / abs(e$values))) %*% t(e$vectors)
}


S          = 10000
start      = S / 10
start_time = Sys.time()


hyper = sample_hyper(S, start, p, 
                     extend_hyper(init, model, matrix(result$par)), 
                     model, Y, X, W)


hyper = t(hyper)
plot.ts(hyper)


Sys.time() - start_time

