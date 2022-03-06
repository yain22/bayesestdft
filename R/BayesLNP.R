BayesLNP = function(y, ini.nu = 1 , S = 1000, mu = 1, sigma.sq = 1){
  # Number of sample size
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
    # B. ESS
    {
      # a. Choose an ellipse centered at mu:
      rho = rnorm(n = 1, mean = mu, sd = sqrt(sigma.sq))

      # c. Define a criterion function :
      alpha = function(eta.new, eta.old){
        e = exp(1)
        #f1 = log( gamma( (e^eta.new + 1)/ 2 ) / (sqrt(e^eta.new * pi) * gamma( (e^eta.new)/2) ) )
        #f2 = log( gamma( (e^eta.old + 1)/ 2 ) / (sqrt(e^eta.old * pi) * gamma( (e^eta.old)/2) ) )
        f1 = lgamma( (e^eta.new + 1)/ 2 ) - ( eta.new/2 + log(pi)/2 + lgamma( (e^eta.new)/ 2 ))
        f2 = lgamma( (e^eta.old + 1)/ 2 ) - ( eta.old/2 + log(pi)/2 + lgamma( (e^eta.old)/ 2 ))
        f3 = ((exp(eta.new)+1)/2)*sum(log(1 + (y^2)/exp(eta.new)))
        f4 = ((exp(eta.old)+1)/2)*sum(log(1 + (y^2)/exp(eta.old)))
        res = min(exp(N*(f1 - f2) -(f3 - f4)),1)
        return(res)
      }

      # c. Choose a threshold and fix :
      u = runif(n = 1, min = 0, max = 1)

      # c. Draw an initial proposal :
      phi = runif(n = 1, min = -pi, max = pi)
      eta.star = ( eta[s] - mu) * cos(phi) + ( rho - mu) * sin(phi) + mu

      # d. ESS core step
      if (u < alpha(eta.new = eta.star, eta.old = eta[s])){
        eta[s+1] = eta.star
      } else {
        # Define a bracket:
        phi.min = -pi ; phi.max = pi

        while(u >= alpha(eta.new = eta.star, eta.old = eta[s])){
          # Shrink the braket and try a new point:
          if (phi>0) {phi.max = phi} else {phi.min = phi}
          phi = runif(n = 1, min = -pi, max = pi)
          eta.star = ( eta[s] - mu) * cos(phi) + ( rho - mu) * sin(phi) + mu
        }
        eta[s+1] = eta.star
      }

    }

    # C. Change of variable
    {
      nu[s+1] = exp(eta[s+1])
    }

  }

  res = nu

  return(res)

}
