
data(optimism)

set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvarSIGN$new(optimism)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
hd                  <- compute_historical_decompositions(run_no1)

set.seed(1)
suppressMessages(
  hd2               <- optimism |>
    specify_bsvarSIGN$new() |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    compute_historical_decompositions()
)



expect_equal(
  length(dim(hd)), length(dim(hd2)),
  info = "compute_historical_decompositions: same output dimensions for normal and pipe workflow."
)

expect_identical(
  hd[3,3,3,3], hd2[3,3,3,3],
  info = "compute_historical_decompositions: identical for normal and pipe workflow."
)
