
test_that(
  "minimum example of specify_bsvarSIGN works", 
  {
    data(oil)
    spec = specify_bsvarSIGN$new(oil)
    expect_identical(
      class(spec)[1],
      "BSVARSIGN"
    )
  }
)
