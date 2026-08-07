S     <- 4L
draws <- array(seq_len(S), dim = c(1, 1, S))

posterior <- list(
  posterior = list(
    w      = matrix(1, nrow = S, ncol = 1),
    A      = draws,
    B      = draws + S,
    hyper  = matrix(seq_len(S) + 2L * S, nrow = 1),
    Q      = draws + 3L * S,
    Sigma  = draws + 4L * S,
    Theta0 = draws + 5L * S
  )
)

expected_draws <- posterior$posterior[-1]

set.seed(1)
rng_before <- .Random.seed
actual     <- bsvarSIGNs:::importance_sampling(posterior)
rng_after  <- .Random.seed

expect_identical(
  actual$posterior[names(expected_draws)],
  expected_draws,
  info = "importance_sampling: unit weights preserve posterior draws."
)

expect_false(
  "w" %in% names(actual$posterior),
  info = "importance_sampling: importance weights are removed."
)

expect_equal(
  actual$posterior$ess,
  S,
  info = "importance_sampling: unit weights imply full effective sample size."
)

expect_identical(
  rng_after,
  rng_before,
  info = "importance_sampling: unit weights do not consume random numbers."
)
