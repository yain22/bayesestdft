library(numDeriv)

BayesJeffreys = function(y, ini.nu = 1 , S = 1000, delta = 0.001){
  # y : Number of samples
  # ini.nu : initial value of MCMC chain
  # S : Number of MCMC sample
  # delta : Step of proposal

  # MCMC sampler
  {
    # Sample size
    N = length(y)

    # Make a room
    nu = rep(0,S)
    eta = rep(0,S) # eta = log(nu)

    # Initial value
    nu[1] = 1
    eta[1] = log(nu[1])

    for (s in 1:(S-1)){

      # A. Change of variable
      {
        eta[s] = log(nu[s])
      }
      # B. MH algorithm
      {
        # a . Define a criterion function :
        gamma_target = function(eta){
          e = exp(1)

          log_likelihood_part = N * (lgamma( (e^eta + 1)/ 2 ) - ( eta/2 + log(pi)/2 + lgamma( (e^eta)/ 2 )) ) -
            ((exp(eta)+1)/2)*sum(log(1 + (y^2)/exp(eta)))

          log_prior_part = (1/2) * (log( (e^eta) / ( e^eta + 3 ) ) + log(trigamma( (e^eta)/2 ) - trigamma( (e^eta + 1)/2 ) -2 * (e^eta + 3) / ((e^eta) * (e^eta + 1)^2 ))   )

          res = - (log_likelihood_part + log_prior_part + eta)
          return(res)
        }

        grad_gamma_target = function(eta){
          res = grad(func = gamma_target, x = eta)
          return(res)
        }

        alpha = function(eta.new, eta.old){
          e = exp(1)

          # Likelihood ratio part
          f1 = lgamma( (e^eta.new + 1)/ 2 ) - ( eta.new/2 + log(pi)/2 + lgamma( (e^eta.new)/ 2 ))
          f2 = lgamma( (e^eta.old + 1)/ 2 ) - ( eta.old/2 + log(pi)/2 + lgamma( (e^eta.old)/ 2 ))

          f3 = ((exp(eta.new)+1)/2)*sum(log(1 + (y^2)/exp(eta.new)))
          f4 = ((exp(eta.old)+1)/2)*sum(log(1 + (y^2)/exp(eta.old)))

          # Prior ratio part
          p1 = log( (e^eta.new) / ( e^eta.new + 3 ) )
          p2 = log( (e^eta.new) / ( e^eta.old + 3 ) )
          p3 = log(trigamma( (e^eta.new)/2 ) - trigamma( (e^eta.new + 1)/2 ) -2 * (e^eta.new + 3) / ((e^eta.new) * (e^eta.new + 1)^2 ))
          p4 = log(trigamma( (e^eta.old)/2 ) - trigamma( (e^eta.old + 1)/2 ) -2 * (e^eta.old + 3) / ((e^eta.old) * (e^eta.old + 1)^2 ))

          Log_of_Likelihood.ratio.part = N*(f1 - f2) - (f3 - f4) +(1/2)*( p1 - p2 + p3 - p4) + eta.new - eta.old
          Log_of_MH.correction.part = dnorm(x = eta.old, mean = eta.new - delta*grad_gamma_target(eta.new), sd = sqrt(2*delta), log = T ) -
            dnorm(x = eta.new, mean = eta.old - delta*grad_gamma_target(eta.old), sd = sqrt(2*delta), log = T )
          res = min( exp (Log_of_Likelihood.ratio.part + Log_of_MH.correction.part) ,1)

          return(res)
        }

        # b. Choose a threshold:
        u = runif(n = 1, min = 0, max = 1)

        # c. Draw an initial proposal :
        eta.star = rnorm(n = 1, mean = eta[s] - delta*grad_gamma_target(eta[s]), sd = sqrt(2*delta))

        # d. ESS core step
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

  }


  res = nu

  return(res)

}
