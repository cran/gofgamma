#' Maximum-likelihood estimation of parameters for the gamma distribution
#'
#' @description
#' The maximum-likelihood estimators for the \code{shape} and \code{scale} parameters of a gamma distribution are computed due to the method of Bhattacharya (2001).
#'
#' @param data vector of positive valued observations
#'
#' @return returns a bivariate vector containing (\code{shape},\code{scale}) estimated parameter vector.
#'
#' @references
#' Bhattacharya, B. (2001) "Testing equality of scale parameters against restricted alternatives for \eqn{m \ge 3} gamma distributions with unknown common shape parameter". Journal of Statistical Computation and Simulations, 69(4):353-368, \href{https://doi.org/10.1080/00949650108812100}{DOI}
#'
#' @examples
#' gamma_est(stats::rgamma(100,shape=3,scale=6))
#'
#' @export
gamma_est <- function(data){
  para = c(0,0)
  if (sum((data<=0))>0){
    warning("negative values in the data!")
  }
  else {
  Rn = log(mean(data)) - mean(log(data))
  if (Rn > 0 && Rn <= 0.5772){
    para[1] = (0.5000876 + 0.1648852*Rn - 0.0544274*(Rn^2)) / Rn
  }
  else if (Rn > 0.5772 && Rn <= 17){
    para[1] = (8.898919 + 9.059950*Rn + 0.9775373*(Rn^2)) / (Rn * (17.79728 + 11.968477*Rn + Rn^2))
  }
  else if (Rn > 17){
    para[1] = 1/Rn
  }

  para[2] = mean(data)/para[1]
  return(para)
  }
}

