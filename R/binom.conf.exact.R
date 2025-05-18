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
