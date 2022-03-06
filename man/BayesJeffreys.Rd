\name{BayesJeffreys}
\alias{BayesJeffreys}
\title{Random walk Metropolis algorithm and Metropolis-adjusted Langevin algorithm to draw posterior samples from the posterior distribution of the degrees of freedom based on the Jeffreys prior}
\usage{
BayesJeffreys(
  y,
  ini.nu = 1,
  S = 1000,
  delta = 0.001,
  sampling.alg = c("MH","MALA")
)
}
\arguments{
\item{y}{- N-dimensional vector of continuous observations}

\item{ini.nu}{- Initial value of posterior sample (Default = 1)}

\item{S}{- Number of posterior samples (Default = 1000)}

\item{delta}{- Step-size of Random walk Metropolis algorithm and Metropolis-adjusted Langevin algorithm (Default = 0.001)}

\item{sampling.alg}{- Option of algorithm specifying Random walk Metropolis algorithm (sampling.alg = "MH") and Metropolis-adjusted Langevin algorithm (sampling.alg = "MALA")}
}
\value{
\item{Output}{- S-dimensional vector of posterior samples}
}
\description{
Random walk Metropolis algorithm and Metropolis-adjusted Langevin algorithm to draw posterior samples from the posterior distribution of the degrees of freedom based on the Jeffreys prior.
}
\examples{
x = rt(n = 100, df = 0.1)
nu1 = BayesJeffreys(x, sampling.alg = "MH")
nu2 = BayesJeffreys(x, sampling.alg = "MALA")
mean(nu1)
mean(nu2)
}
