
# minimum example of estimate.BSVAR and estimate.PosteriorBSVARSIGN works

# run 1
set.seed(123)
data(oil)
spec1 = specify_bsvarSIGN$new(oil)
post1 = estimate(spec1,
                 S = 3,
                 thin = 1,
                 show_progress = FALSE)

# run 2
set.seed(123)
data(oil)
spec2 = specify_bsvarSIGN$new(oil)
post2 = estimate(spec2,
                 S = 3,
                 thin = 1,
                 show_progress = FALSE)

# run 3 (pipe workflow)
set.seed(123)
data(oil)
post3 = oil |>
  specify_bsvarSIGN$new() |>
  estimate(S = 3,
           thin = 1,
           show_progress = FALSE)


# tests
expect_identical(class(post1)[1],
                 class(post2)[1],)

expect_identical(class(post2)[1],
                 class(post3)[1],)

expect_identical(class(post1)[1],
                 "PosteriorBSVARSIGN")

expect_identical(post1$posterior$B[, , 1],
                 post2$posterior$B[, , 1],)

expect_identical(post2$posterior$B[, , 1],
                 post3$posterior$B[, , 1],)

expect_identical(post1$last_draw$B[, , 1],
                 post2$last_draw$B[, , 1],)

expect_identical(post2$last_draw$B[, , 1],
                 post3$last_draw$B[, , 1],)
