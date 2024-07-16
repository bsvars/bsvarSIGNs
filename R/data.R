

#' @title A 5-variable US business cycle data, from 1955 Q1 to 2004 Q4
#'
#' @description A sample data to identify optimism shock.
#'
#' @usage data(optimism)
#'
#' @format A matrix and a \code{ts} object with time series of over two hundred observations on 5 variables:
#' \describe{
#'   \item{productivity}{quarterly factor-utilization-adjusted total factor productivity}
#'   \item{stock_prices}{quarterly end-of-period S&P 500 divided by CPI}
#'   \item{consumption}{quarterly real consumption expenditures on nondurable goods and services}
#'   \item{real_interest_rate}{quarterly real interest rate}
#'   \item{hours_worked}{quarterly hours of all persons in the non-farm business sector}
#' }
#'
#' The series are as described by Beaudry, Nam and Wang (2011) in section 2.2.
#'
#' @references
#' Arias, Jonas E., Juan F. Rubio‐Ramírez, and Daniel F. Waggoner. "Inference based on structural vector autoregressions identified with sign and zero restrictions: Theory and applications." Econometrica 86, no. 2 (2018): 685-720. <doi:10.3982/ECTA14468>
#' 
#' Beaudry, Paul, Deokwoo Nam, and Jian Wang. Do mood swings drive business cycles and is it rational?. No. w17651. National Bureau of Economic Research, 2011. <doi:10.3386/w17651>
#'
#' @source
#' Replication package, \url{https://www.econometricsociety.org/publications/econometrica/2018/03/01/inference-based-structural-vector-autoregressions-identified}
"optimism"


#' @title A 6-variable US monetary policy data, from 1965 Jan to 2007 Aug
#'
#' @description A sample data to identify monetary policy shock.
#'
#' @usage data(monetary)
#'
#' @format A matrix and a \code{ts} object with time series of over two hundred observations on 5 variables:
#' \describe{
#'   \item{gdpc1}{monthly real gross domestic product}
#'   \item{gdpdef}{monthly gross domestic product: implicit price deflator}
#'   \item{cprindex}{monthly consumer price index}
#'   \item{totresns}{monthly reserves of depository institutions}
#'   \item{bognonbr}{monthly non-borrowed reserves of depository institutions}
#'   \item{fedfunds}{monthly federal funds effective rate}
#' }
#'
#' @references
#' Antolín-Díaz & Rubio-Ramírez (2018) Narrative Sign Restrictions for SVARs, American Economic Review, 108(10), 2802-29, <doi:10.1257/aer.20161852>.
#'
#' @source
#' Replication package, \url{https://www.aeaweb.org/articles?id=10.1257/aer.20161852}
"monetary"

