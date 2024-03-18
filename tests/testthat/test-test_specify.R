data(oil)
spec = specify_bsvarSIGN$new(oil)

test_that(
  "minimum example of specify_bsvarSIGN works", 
  {
    expect_identical(
      class(spec)[1],
      "BSVARSIGN"
    )
  }
)
