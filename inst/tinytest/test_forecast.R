
data(optimism)

# for bsvar
set.seed(1)
suppressMessages(
  specification_no1 <- specify_bsvarSIGN$new(optimism, p = 1)
)
run_no1             <- estimate(specification_no1, 3, 1, show_progress = FALSE)
ff                  <- forecast(run_no1, horizon = 2)

set.seed(1)
suppressMessages(
  ff2              <- optimism |>
    specify_bsvarSIGN$new(p = 1) |>
    estimate(S = 3, thin = 1, show_progress = FALSE) |>
    forecast(horizon = 2)
)


expect_identical(
  ff$forecasts[1,1,1], ff2$forecasts[1,1,1],
  info = "forecast: forecast identical for normal and pipe workflow."
)

expect_true(
  is.numeric(ff$forecasts) & is.array(ff$forecasts),
  info = "forecast: returns numeric array."
)


expect_error(
  specify_bsvar$new(us_fiscal_lsuw) |> forecast(horizon = 3),
  info = "forecast: wrong input provided."
)

expect_error(
  forecast(run_no1, horizon = 1.5),
  info = "forecast: specify horizon as integer."
)

