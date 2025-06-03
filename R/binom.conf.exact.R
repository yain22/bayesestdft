#' Exact Binomial Confidence Interval Using Clopper-Pearson Method
#'
#' Computes an exact two-sided confidence interval for a binomial proportion using the Clopper-Pearson method, 
#' based on the F-distribution quantiles. This method guarantees coverage by inverting the binomial test and 
#' is commonly used for small sample sizes or when exact methods are preferred over asymptotic approximations.
#'
#' @param x Integer. Number of observed successes.
#' @param n Integer. Total number of trials.
#' @param alpha Numeric. Significance level for the confidence interval (default is 0.05 for a 95\% CI).
#'
#' @return A named numeric vector with three elements:
#' \describe{
#'   \item{PointEst}{The observed proportion (\code{x/n}).}
#'   \item{Lower}{The lower bound of the exact two-sided confidence interval.}
#'   \item{Upper}{The upper bound of the exact two-sided confidence interval.}
#' }
#'
#' @examples
#' binom.conf.exact(x = 8, n = 10)
#' binom.conf.exact(x = 50, n = 100, alpha = 0.01)
#'
#' @export

binom.conf.exact = function (x, n, alpha = 0.05){
  # Lower bound
  nu1 <- 2 * (n - x + 1)
  nu2 <- 2 * x
  ll = qf(1 - alpha/2, nu1, nu2)
  lb = x/(x + ll * (n - x + 1))
  if (ll == "NaN"){
    lb = 0
  }
  
  # Upper bound
  nu1p <- nu2 + 2
  nu2p <- nu1 - 2
  pp = qf(1 - alpha/2, nu1p, nu2p)
  ub <- ((x + 1) * pp)/(n - x + (x + 1) * pp)
  if (pp == "NaN"){
    ub = 1
  }
  
  res <- c(x/n, lb, ub)
  names(res) = c("PointEst" , "Lower" , "Upper")
  return(res)
}
