% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{monetary}
\alias{monetary}
\title{A 6-variable US monetary policy data, from 1965 Jan to 2007 Aug}
\format{
A matrix and a \code{ts} object with time series of over two hundred observations on 5 variables:
\describe{
  \item{gdpc1}{monthly real gross domestic product}
  \item{gdpdef}{monthly gross domestic product: implicit price deflator}
  \item{cprindex}{monthly consumer price index}
  \item{totresns}{monthly reserves of depository institutions}
  \item{bognonbr}{monthly non-borrowed reserves of depository institutions}
  \item{fedfunds}{monthly federal funds effective rate}
}
}
\source{
Replication package, \url{https://www.aeaweb.org/articles?id=10.1257/aer.20161852}
}
\usage{
data(monetary)
}
\description{
A sample data to identify monetary policy shock.
}
\references{
Antolín-Díaz & Rubio-Ramírez (2018) Narrative Sign Restrictions for SVARs, American Economic Review, 108(10), 2802-29, <doi:10.1257/aer.20161852>.
}
\keyword{datasets}
