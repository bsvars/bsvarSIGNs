
# minimum example of specify_bsvarSIGN works

data(oil)

spec = specify_bsvarSIGN$new(oil)

expect_identical(class(spec)[1],
                 "BSVARSIGN")

expect_identical(dim(spec$data_matrices$Y)[2],
                 dim(spec$data_matrices$X)[2])

expect_identical(length(spec$data_matrices$get_data_matrices()),
                 2L)

expect_identical(dim(spec$prior$A)[2],
                 dim(spec$data_matrices$X)[1])

expect_identical(class(spec$prior$get_prior()),
                 "list")

expect_true(det(spec$prior$B_V_inv) > 0)

expect_true(det(spec$prior$A_V_inv) > 0)

expect_identical(class(spec$identification)[1],
                 "IdentificationBSVARSIGN")

expect_identical(class(spec$identification$get_identification()),
                 "list")

expect_identical(class(spec$starting_values)[1],
                 "StartingValuesBSVAR")


# "example specifying a reproduction of Antolín-Díaz & Rubio-Ramírez (2018, AER)",

data(oil)
sign_narrative = matrix(c(2, -1, 3, 2, 236, 0), ncol = 6)

sign_irf       = array(matrix(c(-1, -1, 1, 1, 1, 1, 1, -1, 1), nrow = 3),
                       dim = c(3, 3, 1))

set.seed(123)

spec           = specify_bsvarSIGN$new(oil,
                                       p              = 12,
                                       sign_irf       = sign_irf,
                                       sign_narrative = sign_narrative)

expect_identical(class(spec)[1],
                 "BSVARSIGN")

expect_identical(dim(spec$data_matrices$Y)[2],
                 dim(spec$data_matrices$X)[2])

expect_identical(length(spec$data_matrices$get_data_matrices()),
                 2L)

expect_identical(dim(spec$prior$A)[2],
                 dim(spec$data_matrices$X)[1])

expect_identical(class(spec$prior$get_prior()),
                 "list")

expect_true(det(spec$prior$B_V_inv) > 0)

expect_true(det(spec$prior$A_V_inv) > 0)

expect_identical(class(spec$identification)[1],
                 "IdentificationBSVARSIGN")

expect_identical(class(spec$identification$get_identification()),
                 "list")

expect_identical(class(spec$starting_values)[1],
                 "StartingValuesBSVAR")
