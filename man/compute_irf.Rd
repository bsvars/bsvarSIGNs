% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{compute_irf}
\alias{compute_irf}
\title{Compute impulse response functions}
\usage{
compute_irf(posterior, horizon = 24)
}
\arguments{
\item{posterior}{posterior estimation outcome}

\item{horizon}{a positive integer number denoting the forecast horizon for the impulse responses computations.}
}
\description{
Given posterior draws compute impulse response functions with bsvars
}
