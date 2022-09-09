BayesGA = function(y, ini.nu = 1 , S = 1000, delta = 0.001, a = 1, b = 0.1){

  # Sample size
  N = length(y)

  # Make a room
  nu = rep(0,S)
  eta = rep(0,S) # eta = log(nu)

  # Initial value
  nu[1] = ini.nu
  eta[1] = log(nu[1])

  for (s in 1:(S-1)){

    # A. Change of variable
    {
      eta[s] = log(nu[s])
    }
    # B. MH algorithm
    {
      # a . Define a criterion function :
      alpha = function(eta.new, eta.old){
        e = exp(1)

        # Likelihood ratio part
        f1 = lgamma( (e^eta.new + 1)/ 2 ) - ( eta.new/2 + log(pi)/2 + lgamma( (e^eta.new)/ 2 ))
        f2 = lgamma( (e^eta.old + 1)/ 2 ) - ( eta.old/2 + log(pi)/2 + lgamma( (e^eta.old)/ 2 ))

        f3 = ((exp(eta.new)+1)/2)*sum(log(1 + (y^2)/exp(eta.new)))
        f4 = ((exp(eta.old)+1)/2)*sum(log(1 + (y^2)/exp(eta.old)))

        # Prior ratio part
        p1 = a*(eta.new - eta.old)
        p2 = b*(exp(eta.new) - exp(eta.old))

        # Result
        res = min(exp(
          N*(f1 - f2) - (f3 - f4) +
            ( p1 - p2 )
        ),1)
        return(res)
      }

      # b. Choose a threshold:
      u = runif(n = 1, min = 0, max = 1)

      # c. Draw an initial proposal :
      eta.star = rnorm(n = 1, mean = eta[s], sd = sqrt(2*delta))

      # d. MH core step
      if (u < alpha(eta.new = eta.star, eta.old = eta[s])){
        eta[s+1] = eta.star
      } else {
        eta[s+1] = eta[s]
      }

    }

    # C. Change of variable
    {
      nu[s+1] = exp(eta[s+1])
    }

  }

  res = list(nu = nu)

  return(res)

}
