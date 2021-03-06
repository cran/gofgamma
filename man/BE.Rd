% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/teststatistics.R
\name{BE}
\alias{BE}
\title{statistic of the Betsch-Ebner test}
\usage{
BE(data, k_estimator, a)
}
\arguments{
\item{data}{a vector of positive numbers. NOTE: \code{data} has to be the rescaled data, i.e. devided by the estimated \code{scale} parameter!}

\item{k_estimator}{value of the estimated \code{shape} parameter.}

\item{a}{positive tuning parameter.}
}
\value{
value of the test statistic
}
\description{
This function computes the statistic of the goodness-of-fit test for the gamma family due to Betsch and Ebner (2019).
}
\details{
The test is of weighted \eqn{L^2} type and uses a characterization of the distribution function of the gamma distribution. Values of \code{k_estimator} are found by \code{\link{gamma_est}}.
}
\examples{
X=stats::rgamma(20,3,6)
BE(X,k_estimator=gamma_est(X)[1],a=2)

}
\references{
Betsch, S., Ebner, B. (2019) "A new characterization of the Gamma distribution and associated goodness of fit tests", Metrika, 82(7):779-806. \href{https://doi.org/10.1007/s00184-019-00708-7}{DOI}
}
