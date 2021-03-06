% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/teststatistics.R
\name{HME1}
\alias{HME1}
\title{statistic of the first Henze-Meintanis-Ebner goodness-of-fit test for the gamma family}
\usage{
HME1(data, k_estimator, a = 1)
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
This function computes the goodness-of-fit test statistic for the gamma family due to the first test in Henze, Meintanis and Ebner (2012).
}
\details{
The test statistic is of weighted \eqn{L^2} type and uses a characterization of the distribution function of the gamma distribution.
}
\examples{
X=stats::rgamma(20,3,6)
HME1(X,k_estimator=gamma_est(X)[1],a=1)

}
\references{
Henze, N., Meintanis, S.G., Ebner, B. (2012) "Goodness-of-fit tests for the Gamma distribution based on the empirical Laplace transform". Communications in Statistics - Theory and Methods, 41(9): 1543-1556. \href{https://doi.org/10.1080/03610926.2010.542851}{DOI}
}
