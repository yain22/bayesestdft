\name{BayesLNP}
\alias{BayesLNP}
\title{Elliptical slice sampler to draw posterior samples from the posterior distribution of the degrees of freedom based on a log-normal prior}
\usage{
BayesLNP(
  y,
  ini.nu = 1,
  S = 1000,
  mu = 1,
  sigma.sq = 1
)
}
\arguments{
\item{y}{- N-dimensional vector of continuous observations}

\item{ini.nu}{- Initial value of posterior sample (Default = 1)}

\item{S}{- Number of posterior samples (Default = 1000)}

\item{mu}{- Mean of the log-normal prior density (Default = 1)}

\item{sigma.sq}{- Variance of the log-normal prior density (Default = 1)}
}
\value{
\item{Output}{- S-dimensional vector of posterior samples}
}
\description{
Elliptical slice sampler to draw posterior samples from the posterior distribution of the degrees of freedom provided a log-normal prior.
}
\examples{
# Example 1:
x = rt(n = 100, df = 0.1)
nu = BayesLNP(x)
mean(nu)

# Example 2:
library(dplyr)
data(index_return)
index_return_US <- index_return %>% filter(Country == "United States")
x = index_return_US$log_return_rate
nu.LNP.ESS = BayesLNP(x)
mean(nu.LNP.ESS)
}
