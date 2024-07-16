
# minimum example of specify_bsvarSIGN works

data(optimism)

spec = specify_bsvarSIGN$new(optimism)

expect_identical(class(spec)[1],
                 "BSVARSIGN")

expect_identical(dim(spec$data_matrices$Y)[2],
                 dim(spec$data_matrices$X)[2])

expect_identical(length(spec$data_matrices$get_data_matrices()),
                 2L)

expect_identical(dim(spec$prior$A)[2], dim(spec$data_matrices$X)[1])

expect_identical(class(spec$prior$get_prior()),
                 "list")

expect_true(det(spec$prior$V) > 0)

expect_true(det(spec$prior$S) > 0)

expect_identical(class(spec$identification)[1],
                 "IdentificationBSVARSIGN")

expect_identical(class(spec$identification$get_identification()),
                 "list")

expect_identical(class(spec$starting_values)[1],
                 "StartingValuesBSVAR")


# "example specifying a reproduction of Antolín-Díaz & Rubio-Ramírez (2018, AER)",

data(optimism)

sign_irf       = matrix(c(0, 1, rep(NA, 23)), 5, 5)

set.seed(123)

spec           = specify_bsvarSIGN$new(optimism,
                                       p              = 12,
                                       sign_irf       = sign_irf)

expect_identical(class(spec)[1],
                 "BSVARSIGN")

expect_identical(dim(spec$data_matrices$Y)[2],
                 dim(spec$data_matrices$X)[2])

expect_identical(length(spec$data_matrices$get_data_matrices()),
                 2L)

expect_identical(dim(spec$prior$A)[2], dim(spec$data_matrices$X)[1])

expect_identical(class(spec$prior$get_prior()),
                 "list")

expect_true(det(spec$prior$V) > 0)

expect_true(det(spec$prior$S) > 0)

expect_identical(class(spec$identification)[1],
                 "IdentificationBSVARSIGN")

expect_identical(class(spec$identification$get_identification()),
                 "list")

expect_identical(class(spec$starting_values)[1],
                 "StartingValuesBSVAR")

