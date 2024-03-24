
test_that(
  "minimum example of estimate.BSVARSIGN works", 
  {
    data(oil)
    spec = specify_bsvarSIGN$new(oil)
    burn = estimate(spec, S = 10)
    expect_identical(
      class(burn)[1],
      "PosteriorBSVARSIGN"
    )
  }
)

test_that(
  "minimum example of estimate.PosteriorBSVARSIGN works", 
  {
    data(oil)
    spec = specify_bsvarSIGN$new(oil)
    burn = estimate(spec, S = 10)
    post = estimate(burn, S = 10)
    expect_identical(
      class(post)[1],
      "PosteriorBSVARSIGN"
    )
  }
)