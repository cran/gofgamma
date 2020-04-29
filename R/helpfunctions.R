#needed functions

erf <- function(x){2*stats::pnorm(sqrt(2)*x) - 1}

varphi <- function(a,j,k){
  if (is.nan ((1- erf ((j+k)/(2*sqrt(a))))*exp (((j+k)^2)/(4*a)))==TRUE) return(0)
  else (1-erf((j+k)/(2*sqrt(a))))*exp (((j+k)^2)/(4*a))
}

#' Print method for tests of Gamma distribution
#' @description
#' Printing objects of class "gofgamma".
#'
#' @param x object of class "gofgamma".
#' @param ... further arguments to be passed to or from methods.
#'
#' @details
#' A \code{gofgamma} object is a named list of numbers and character string, supplemented with \code{test} (the name of the teststatistic). \code{test} is displayed as a title.
#' The remaining elements are given in an aligned "name = value" format.
#'
#' @return
#' the argument x, invisibly, as for all \link{print} methods.
#'
#' @examples
#' \donttest{print(test.BE(rgamma(20,1)))}
#'
#' @export
print.gofgamma <- function(x, ...){
  cat("\n \n")
  cat("------------------------------------------------------------------------- \n")
  cat("\n")
  cat("         One-sample test for Gamma distribution with the", x$Test, " teststatistic.\n"  )
  cat("\n")
  if(!is.null(x$parameter)){cat("tuning parameter            =", x$param," \n")}
  cat( x$Test, "                         = ", x$T.value, " \n")
  cat("critical value              =  ", x$cv, " \n")
  cat("estimated parameters        = ", x$par.est, " \n")
  cat("significance level          = ", x$sig.level, " \n")
  cat("number of bootstrap samples = ", x$boot, " \n \n")
  if (x$Decision == TRUE) {
    cat("Gamma-hypothesis has to be rejected  \n")
  } else {
    cat("Gamma-hypothesis cannot be rejected  \n")
  }
  cat("\n")
  cat("\n")
  cat("------------------------------------------------------------------------- \n")
}

