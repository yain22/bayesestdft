\name{BayesGA}
\alias{BayesGA}
\title{Random walk Metropolis algorithm to draw posterior samples from the posterior distribution of the degrees of freedom based on a gamma prior}
\usage{
BayesGA(
  y,
  ini.nu = 1,
  S = 1000,
  a = 1,
  b = 0.1,
)
}
\arguments{
\item{y}{- N-dimensional vector of continuous observations}

\item{ini.nu}{- Initial value of posterior sample (Default = 1)}

\item{S}{- Number of posterior samples (Default = 1000)}

\item{a}{- Shape parameter of the gamma prior density (Default = 1)}

\item{b}{- Rate parameter of the gamma prior density (Default = 0.1)}
}
\value{
\item{Output}{- S-dimensional vector of posterior samples}
}
\description{
Random walk Metropolis algorithm to draw posterior samples from the posterior distribution of the degrees of freedom based on a gamma prior.
}
\examples{
# Example 1 (Simulation study):

x = rt(n = 100, df = 0.1)
nu = BayesGA(x)
mean(nu)

# Example 2 (Log-return data):

library(dplyr)
data(index_return)
index_return_US <- filter(index_return, Country == "United States")
x = index_return_US$log_return_rate
nu = BayesGA(x)
mean(nu)
}
