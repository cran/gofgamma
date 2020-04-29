#' statistic of the Betsch-Ebner test
#'
#' @description
#' This function computes the statistic of the goodness-of-fit test for the gamma family due to Betsch and Ebner (2019).
#'
#' @param data           a vector of positive numbers. NOTE: \code{data} has to be the rescaled data, i.e. devided by the estimated \code{scale} parameter!
#' @param k_estimator    value of the estimated \code{shape} parameter.
#' @param a              positive tuning parameter.
#'
#' @return value of the test statistic
#'
#'
#' @details
#' The test is of weighted \eqn{L^2} type and uses a characterization of the distribution function of the gamma distribution. Values of \code{k_estimator} are found by \code{\link{gamma_est}}.
#'
#' @references
#' Betsch, S., Ebner, B. (2019) "A new characterization of the Gamma distribution and associated goodness of fit tests", Metrika, 82(7):779-806. \href{https://doi.org/10.1007/s00184-019-00708-7}{DOI}
#'
#' @examples
#' X=stats::rgamma(20,3,6)
#' BE(X,k_estimator=gamma_est(X)[1],a=2)
#'
#' @export
BE <- function(data, k_estimator,a)
{
  n = length(data)
  data = sort(data)

  B = -(k_estimator - 1)/data + 1
  Stati = 2/n * sum( sapply(1:(n-1), function(j){sum( sapply((j+1):n, function(l){ exp(-a*data[l])/a *(data[j] - k_estimator) * ( -1/a * B[l] - 1 ) + B[j]*B[l]*2/(a^3) + exp(-a*data[j])/a * B[l] * ( (k_estimator - 2 - data[j])/a - B[j]*2/(a^2) - data[j] ) } ) )} ) )
  Stati = Stati + 1/n * sum(exp(-a*data)/a * ( 2*(k_estimator - 1) - 2*data + 1 + B^2 * (-2*data/a - 2/(a^2))) + 2*B^2 /(a^3))
  return(Stati)
}

#' statistic of the Kolmogorov-Smirnov goodness-of-fit test for the gamma family
#'
#' @description
#' This function computes the goodness-of-fit test statistic for the gamma family in the spirit of Kolmogorov and Smirnov. Note that this tests the composite hypothesis of fit to the family of gamma distributions, i.e. a bootstrap procedure is implemented to perform the test, see \code{\link{crit.values}}.
#'
#' @param data           a vector of positive numbers. NOTE: \code{data} has to be the rescaled data, i.e. devided by the estimated \code{scale} parameter!
#' @param k_estimator    value of the estimated \code{shape} parameter.
#'
#' @return value of the test statistic
#'
#' @details
#' The Kolmogorov-Smirnov test is computed as described in Henze et. al. (2012). Values of \code{k_estimator} are found by \code{\link{gamma_est}}.
#'
#' @references
#' Henze, N., Meintanis, S.G., Ebner, B. (2012) "Goodness-of-fit tests for the Gamma distribution based on the empirical Laplace transform". Communications in Statistics - Theory and Methods, 41(9): 1543-1556. \href{https://doi.org/10.1080/03610926.2010.542851}{DOI}
#'
#' @examples
#' X=stats::rgamma(20,3,6)
#' KS(X,k_estimator=gamma_est(X)[1])
#'
#' @export
KS <- function(data,k_estimator){
  n=length(data)
  D <- rep(0,n)
  y <- stats::pgamma(sort(data),k_estimator,1)
  for ( i in 1:n) {
    D[i] <- max(i/n-y[i],y[i]-(i-1)/n)}
  return(max(D))
}



#' statistic of the Cramer-von Mises goodness-of-fit test for the gamma family
#'
#' @description
#' This function computes the goodness-of-fit test statistic for the gamma family in the spirit of Cramer and von Mises. Note that this tests the composite hypothesis of fit to the family of gamma distributions, i.e. a bootstrap procedure is implemented to perform the test, see \code{\link{crit.values}}.
#'
#' @param data           a vector of positive numbers. NOTE: \code{data} has to be the rescaled data, i.e. devided by the estimated \code{scale} parameter!
#' @param k_estimator    value of the estimated \code{shape} parameter.
#'
#' @return value of the test statistic
#'
#' @details
#' The Cramér-von Mises test is computed as described in Henze et. al. (2012). Values of \code{k_estimator} are found by \code{\link{gamma_est}}.
#'
#' @references
#' Henze, N., Meintanis, S.G., Ebner, B. (2012) "Goodness-of-fit tests for the Gamma distribution based on the empirical Laplace transform". Communications in Statistics - Theory and Methods, 41(9): 1543-1556. \href{https://doi.org/10.1080/03610926.2010.542851}{DOI}
#'
#' @examples
#' X=stats::rgamma(20,3,6)
#' CM(X,k_estimator=gamma_est(X)[1])
#'
#' @export
CM <- function(data,k_estimator){
  A <- 0
  n=length(data)
  y <- sort(data)
  for(j in 1:n){
    A=A+(stats::pgamma(y[j],k_estimator,1)-(2*j -1)/(2*n))^2
  }
  return(1/(12*n)+A)
}


#' statistic of the Anderson-Darling goodness-of-fit test for the gamma family
#'
#' @description
#' This function computes the goodness-of-fit test statistic for the gamma family in the spirit of Anderson and Darling. Note that this tests the composite hypothesis of fit to the family of gamma distributions, i.e. a bootstrap procedure is implemented to perform the test, see \code{\link{crit.values}}.
#'
#' @param data           a vector of positive numbers. NOTE: \code{data} has to be the rescaled data, i.e. devided by the estimated \code{scale} parameter!
#' @param k_estimator    value of the estimated \code{shape} parameter.
#'
#' @return value of the test statistic
#'
#' @details
#' The Anderson-Darling  test is computed as described in Henze et. al. (2012). Values of \code{k_estimator} are found by \code{\link{gamma_est}}.
#'
#' @references
#' Henze, N., Meintanis, S.G., Ebner, B. (2012) "Goodness-of-fit tests for the Gamma distribution based on the empirical Laplace transform". Communications in Statistics - Theory and Methods, 41(9): 1543-1556. \href{https://doi.org/10.1080/03610926.2010.542851}{DOI}
#'
#' @examples
#' X=stats::rgamma(20,3,6)
#' AD(X,k_estimator=gamma_est(X)[1])
#'
#' @export
AD <- function(data,k_estimator){
  A <- 0
  n=length(data)
  y <- sort(data)
  for(j in 1:n){
    A=A+(2*j-1)*log(stats::pgamma(y[j],k_estimator,1))+(2*(n-j)+1)*log(1-stats::pgamma(y[j],k_estimator,1))
  }
  return(-n-1/n*A)
}



#' statistic of the Watson goodness-of-fit test for the gamma family
#'
#' @description
#' This function computes the goodness-of-fit test statistic for the gamma family in the spirit of Watson. Note that this tests the composite hypothesis of fit to the family of gamma distributions, i.e. a bootstrap procedure is implemented to perform the test, see \code{\link{crit.values}}.
#'
#' @param data           a vector of positive numbers. NOTE: \code{data} has to be the rescaled data, i.e. devided by the estimated \code{scale} parameter!
#' @param k_estimator    value of the estimated \code{shape} parameter.
#'
#' @return value of the test statistic
#'
#' @details
#' The Watson test is computed as described in Henze et. al. (2012). Values of \code{k_estimator} are found by \code{\link{gamma_est}}.
#'
#' @references
#' Henze, N., Meintanis, S.G., Ebner, B. (2012) "Goodness-of-fit tests for the Gamma distribution based on the empirical Laplace transform". Communications in Statistics - Theory and Methods, 41(9): 1543-1556. \href{https://doi.org/10.1080/03610926.2010.542851}{DOI}
#'
#' @examples
#' X=stats::rgamma(20,3,6)
#' WA(X,k_estimator=gamma_est(X)[1])
#'
#' @export
WA <- function(data,k_estimator){
  A <- 0
  n=length(data)
  y <- sort(data)
  for (j in 1:n){
    A=A+stats::pgamma(y[j],k_estimator ,1)/n
  }
  return(CM(data,k_estimator)-n*(A -0.5) ^2)
}



#' statistic of the first Henze-Meintanis-Ebner goodness-of-fit test for the gamma family
#'
#' @description
#' This function computes the goodness-of-fit test statistic for the gamma family due to the first test in Henze, Meintanis and Ebner (2012).
#'
#' @param data           a vector of positive numbers. NOTE: \code{data} has to be the rescaled data, i.e. devided by the estimated \code{scale} parameter!
#' @param k_estimator    value of the estimated \code{shape} parameter.
#' @param a              positive tuning parameter.
#'
#' @return value of the test statistic
#'
#' @details
#' The test statistic is of weighted \eqn{L^2} type and uses a characterization of the distribution function of the gamma distribution.
#'
#' @references
#' Henze, N., Meintanis, S.G., Ebner, B. (2012) "Goodness-of-fit tests for the Gamma distribution based on the empirical Laplace transform". Communications in Statistics - Theory and Methods, 41(9): 1543-1556. \href{https://doi.org/10.1080/03610926.2010.542851}{DOI}
#'
#' @examples
#' X=stats::rgamma(20,3,6)
#' HME1(X,k_estimator=gamma_est(X)[1],a=1)
#'
#' @export
HME1 <-function(data,k_estimator,a=1){
  n=length(data)
  y=data
  Ts <- 0
  for (i in 1:n){
    for (j in 1:n){
      Ts=Ts+(y[i]*y[j]-k_estimator*(y[i]+y[j])+k_estimator^2)/(y[i]+y[j]+a)+(2*y[i]*y[j]-k_estimator*(y[i]+y[j]))/((y[i]+y[j]+a)^2)+(2*y[i]*y[j])/((y[i]+y[j]+a)^3)
    }
  }
  Tn=Ts/n
  return(Tn)
}



#' statistic of the second Henze-Meintanis-Ebner goodness-of-fit test for the gamma family
#'
#' @description
#' This function computes the goodness-of-fit test statistic for the gamma family due to the second test in Henze, Meintanis and Ebner (2012).
#'
#' @param data           a vector of positive numbers. NOTE: \code{data} has to be the rescaled data, i.e. devided by the estimated \code{scale} parameter!
#' @param k_estimator    value of the estimated \code{shape} parameter.
#' @param a              positive tuning parameter.
#'
#' @return value of the test statistic
#'
#' @details
#' The test statistic is of weighted \eqn{L^2} type and uses a characterization of the distribution function of the gamma distribution.
#'
#' @references
#' Henze, N., Meintanis, S.G., Ebner, B. (2012) "Goodness-of-fit tests for the Gamma distribution based on the empirical Laplace transform". Communications in Statistics - Theory and Methods, 41(9): 1543-1556. \href{https://doi.org/10.1080/03610926.2010.542851}{DOI}
#'
#' @examples
#' X=stats::rgamma(20,3,6)
#' HME2(X,k_estimator=gamma_est(X)[1],a=1)
#'
#' @export
HME2 <-function(data,k_estimator,a=4){
  n=length(data)
  Y=data
  Ts <- 0
  for (i in 1:n){
    for (j in 1:n){
      Ts=Ts+1/(2*n)*sqrt(pi/a)*(Y[i]*Y[j]+k_estimator^2-k_estimator*(Y[i]+Y[j]))*varphi(a,Y[i],Y[j])+1/(4*a*n)*(2*Y[i]*Y[j]-k_estimator*(Y[i]+Y[j]))*(2-sqrt(pi/a)*(Y[i]+Y[j])*varphi(a,Y[i],Y[j]))+1/(8*n*a^2)*Y[i]*Y[j]*((sqrt(pi/a)*(Y[i]+Y[j])^2+2*sqrt(pi*a))*varphi(a,Y[i],Y[j])-2*(Y[i]+Y[j]))
    }
  }
  return(Ts)
}


