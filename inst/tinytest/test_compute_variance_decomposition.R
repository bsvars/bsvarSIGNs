
data(optimism)

set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvarSIGN$new(optimism)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
fevd                <- compute_variance_decompositions(run_no1, horizon = 2)

set.seed(1)
suppressMessages(
  fevd2               <- optimism |>
    specify_bsvarSIGN$new() |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    compute_variance_decompositions(horizon = 2)
)

expect_error(
  compute_variance_decompositions(run_no1),
  info = "compute_variance_decompositions: specify horizon."
)

expect_equal(
  sum(fevd[1,,1,1]), 100,
  info = "compute_variance_decompositions: sum to 100%."
)

expect_identical(
  fevd[3,3,3,3], fevd2[3,3,3,3],
  info = "compute_variance_decompositions: identical for normal and pipe workflow."
)


