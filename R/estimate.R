
#' @title Bayesian estimation of a Structural Vector Autoregression
#' with traditional and narrative sign restrictions via Gibbs sampler
#'
#' @description Estimates Bayesian Structural Vector Autoregression model
#' using the Gibbs sampler proposed by Waggoner & Zha (2003) with traditional sign restrictions 
#' following Rubio-Ramírez, Waggoner & Zha (2010) and narrative sign restrictions 
#' following Antolín-Díaz & Rubio-Ramírez (2018). Additionally, the parameter matrices \eqn{A} and \eqn{B}
#' follow a Minnesota prior and generalised-normal prior distributions respectively with the matrix-specific
#' overall shrinkage parameters estimated using a hierarchical prior distribution. 
#' 
#' Given sign restrictions, in each Gibbs sampler iteration, the sampler draws rotation matrix 
#' \eqn{Q} uniformly from the space of \code{NxN} orthogonal matrices and checks if the sign restrictions
#' are satisfied. If a valid \eqn{Q} is found within \code{max_tries} (defined in \code{specify_bsvarSIGN}),
#' the sampler saves the current \eqn{A} and \eqn{B} draw and proceeds to the next iteration.
#' Otherwise, the sampler then proceeds to next iteration without saving the current \eqn{A} and \eqn{B} draw.
#' If a narrative sign restriction is given, the posterior
#' draws are resampled with \code{algorithm 1} in Antolín-Díaz & Rubio-Ramírez (2018).
#' 
#' See section \bold{Details} for the model equations.
#' 
#' @details 
#' The Structural VAR model is given by the reduced form equation:
#' \deqn{Y = AX + E}
#' where \eqn{Y} is an \code{NxT} matrix of dependent variables, \eqn{X} is a \code{KxT} matrix of explanatory variables, 
#' \eqn{E} is an \code{NxT} matrix of reduced form error terms, and \eqn{A} is an \code{NxK} matrix of
#' autoregressive slope coefficients and parameters on deterministic terms in \eqn{X}.
#' 
#' The structural equation is given by
#' \deqn{BE = U}
#' where \eqn{U} is an \code{NxT} matrix of structural form error terms, and
#' \eqn{B} is an \code{NxN} matrix of contemporaneous relationships. More specifically,
#' \deqn{B = Q'P}
#' where \eqn{Q} is an \code{NxN} rotation matrix and \eqn{P} is an \code{NxN} lower triangular matrix.
#' 
#' Finally, the structural shocks, \code{U}, are temporally and contemporaneously independent and jointly normally distributed with zero mean and unit variances.
#' 
#' @param specification an object of class \code{BSVARSIGN} generated using the \code{specify_bsvarSIGN$new()} function.
#' @param S a positive integer, the number of posterior draws to be generated
#' @param thin a positive integer, specifying the frequency of MCMC output thinning
#' @param show_progress a logical value, if \code{TRUE} the estimation progress bar is visible
#' 
#' @return An object of class \code{PosteriorBSVARSIGN} containing the Bayesian estimation output and containing two elements:
#' 
#'  \code{posterior} a list with a collection of \code{S} draws from the posterior distribution generated via Gibbs sampler containing:
#'  \describe{
#'  \item{A}{an \code{NxKxS} array with the posterior draws for matrix \eqn{A}}
#'  \item{B}{an \code{NxNxS} array with the posterior draws for matrix \eqn{B}}
#'  \item{hyper}{a \code{5xS} matrix with the posterior draws for the hyper-parameters of the hierarchical prior distribution}
#'  \item{skipped}{an integer of the total skipped iterations,
#'  the Gibbs sampler performs a total of S+skipped iterations,
#'  when the sampler does not find a valid rotation matrix \code{Q} within \code{max_tries},
#'  the current iteration is skipped (i.e. the current draw of \code{A,B} is not saved).
#'  A message is shown when skipped/(skipped+S/thin) > 0.05, where S/thin is the 
#'  total number of draws returned.
#'  }
#' }
#' 
#' \code{last_draw} an object of class \code{BSVARSIGN} with the last draw of 
#' the current MCMC run as the starting value to be passed to the continuation 
#' of the MCMC estimation using \code{estimate()}. 
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}, Xiaolei Wang \email{adamwang15@gmail.com}
#' 
#' @references 
#' 
#'  Antolín-Díaz & Rubio-Ramírez (2018) Narrative Sign Restrictions for SVARs, American Economic Review, 108(10), 2802-29, <doi:10.1257/aer.20161852>.
#'  
#'  Arias, Rubio-Ramírez, & Waggoner (2018), Inference Based on Structural Vector Autoregressions Identified With Sign and Zero Restrictions: Theory and Applications, Econometrica, 86(2), 685-720, <doi:10.3982/ECTA14468>.
#'  
#'  Giannone, Lenza, Primiceri (2015) Prior Selection for Vector Autoregressions, Review of Economics and Statistics, 97(2), 436-451 <doi:10.1162/REST_a_00483>.
#'  
#'  Rubio-Ramírez, Waggoner & Zha (2010) Structural Vector Autoregressions: Theory of Identification and Algorithms for Inference, The Review of Economic Studies, 77(2), 665-696, <doi:10.1111/j.1467-937X.2009.00578.x>.
#'  
#' @method estimate BSVARSIGN
#' 
#' @examples
#' # investigate the effects of the optimism shock
#' data(optimism)
#'
#' # specify identifying restrictions:
#' # + no effect on productivity (zero restriction)
#' # + positive effect on stock prices (positive sign restriction) 
#' sign_irf       = matrix(c(0, 1, rep(NA, 23)), 5, 5)
#' 
#' # specify the model and set seed
#' set.seed(123)
#' specification  = specify_bsvarSIGN$new(optimism * 100,
#'                                        p        = 12,
#'                                        sign_irf = sign_irf)
#'                                        
#' # estimate the model
#' posterior      = estimate(specification, S = 10)
#' 
#' @export
estimate.BSVARSIGN = function(specification, S, thin = 1, show_progress = TRUE) {
  
  # get the inputs to estimation
  # prior               = specification$last_draw$prior$get_prior()
  prior               = specification$prior$get_prior()
  identification      = specification$identification$get_identification()
  max_tries           = identification$max_tries
  max_tries           = ifelse(max_tries == Inf, 0, max_tries)
  data_matrices       = specification$data_matrices$get_data_matrices()
  p                   = specification$p
  
  prior$B             = t(prior$A)
  prior$Ysoc          = t(prior$Ysoc)
  prior$Xsoc          = t(prior$Xsoc)
  prior$Ysur          = t(prior$Ysur)
  prior$Xsur          = t(prior$Xsur)
  Y                   = t(data_matrices$Y)
  X                   = t(data_matrices$X)
  
  Z                   = get_Z(identification$sign_irf)
  sign                = identification$sign_irf
  sign[is.na(sign)]   = 0
  
  n_narratives        = length(identification$sign_narrative)
  get_type            = list("S" = 1, "A" = 2, "B" = 3)
  if (n_narratives > 0) {
    narrative         = matrix(NA, n_narratives, 6)
    for (i in 1:n_narratives) {
      narrative_list  = identification$sign_narrative[[i]]
      narrative[i, 1] = get_type[[narrative_list$type]]
      narrative[i, 2] = narrative_list$sign
      narrative[i, 3] = narrative_list$var
      narrative[i, 4] = narrative_list$shock
      narrative[i, 5] = narrative_list$start - p
      narrative[i, 6] = narrative_list$periods - 1
    }
  } else {
    narrative         = t(c(0, 1, 1, 1, 1, 1))
  }
  struc               = identification$sign_structural
  struc[is.na(struc)] = 0

  # estimation
  qqq                 = .Call(`_bsvarSIGNs_bsvar_sign_cpp`, S, p, Y, X, 
                              sign, narrative, struc, Z, prior, 
                              show_progress, thin, max_tries)
  
  specification$starting_values$set_starting_values(qqq$last_draw)
  output              = specify_posterior_bsvarSIGN$new(specification, qqq$posterior)
  output              = importance_sampling(output)
  
  return(output)
}


