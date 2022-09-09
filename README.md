# bayesestdft: An R package for Bayesian estimation of the number of degrees of the freedom of the Student's t-distribution

## Contents
* [Overview](#overview)
* [Installation](#installation)
* [Jeffreys prior](#jeffreys-prior)
* [Exponential prior](#exponential-prior)
* [Gamma prior](#gamma-prior)
* [Log-normal prior](#log-normal-prior)


## Overview
An R package `bayesestdft` includes tools to implement Bayesian estimation of the number of degrees of the freedom of the Student's t-distribution, developed by [Dr. Se Yoon Lee](https://sites.google.com/view/seyoonlee) (seyoonlee.stat.math@gmail.com). The package was developed to analyze simulated and real data from the published article ["The Use of a Log-Normal Prior for the Student t-Distribution" Axioms 2022, 11, 462"](https://www.mdpi.com/2075-1680/11/9/462). Readers can see the paper for technical details about the package. At current version, the main functions are `BayesLNP` and `BayesJeffreys` that implement Markov Chain Monte Carlo algorithms to sample from the number of degrees of the freedom of the Student's t-distribution. To operatre the function `BayesJeffreys`, user needs to install `R library(numDeriv)`. See the [Slides](https://github.com/yain22/bayesestdft/blob/master/doc/Explaining%20R%20Package%20bayesestdft.pdf) that summarized the technical parts of the R package.

## Required R version
```r
R version 4.0.4 (or higher)
```

## Installation

```r
library(devtools)
devtools::install_github("yain22/bayesestdft")
library(bayesestdft)
```

## Jeffreys prior

$$ \pi(x) $$


#### Estimation of the degrees of freedom from simulated data

```r
x = rt(n = 100, df = 0.1)
nu1 = BayesJeffreys(x, sampling.alg = "MH")
nu2 = BayesJeffreys(x, sampling.alg = "MALA")
mean(nu1)
mean(nu2)
```

#### Estimation of the degrees of freedom of daily log-return rate of S&P500 index time series data 

```r
library(dplyr)
data(index_return)
index_return_US <- filter(index_return, Country == "United States")
x = index_return_US$log_return_rate
nu1 = BayesJeffreys(x, sampling.alg = "MH")
nu2 = BayesJeffreys(x, sampling.alg = "MALA")
mean(nu1)
mean(nu2)
```

## Exponential prior

#### Estimation of the degrees of freedom from simulated data

```r
x = rt(n = 100, df = 0.1)
nu = BayesGA(x, a = 1, b = 0.1)
mean(nu)
```
#### Estimation of the degrees of freedom of daily log-return rate of S&P500 index time series data 

```r
library(dplyr)
data(index_return)
index_return_US <- filter(index_return, Country == "United States")
x = index_return_US$log_return_rate
nu = BayesGA(x, a = 1, b = 0.1)
mean(nu)
```

## Gamma prior

#### Estimation of the degrees of freedom from simulated data

```r
x = rt(n = 100, df = 0.1)
nu = BayesGA(x, a = 2, b = 0.1)
mean(nu)
```
#### Estimation of the degrees of freedom of daily log-return rate of S&P500 index time series data 

```r
library(dplyr)
data(index_return)
index_return_US <- filter(index_return, Country == "United States")
x = index_return_US$log_return_rate
nu = BayesGA(x, a = 2, b = 0.1)
mean(nu)
```


## Log-normal prior

#### Estimation of the degrees of freedom from simulated data

```r
x = rt(n = 100, df = 0.1)
nu = BayesLNP(x)
mean(nu)
```
#### Estimation of the degrees of freedom of daily log-return rate of S&P500 index time series data 

```r
library(dplyr)
data(index_return)
index_return_US <- filter(index_return, Country == "United States")
x = index_return_US$log_return_rate
nu = BayesLNP(x)
mean(nu)
```

## References

[1] [Se Yoon Lee. (2022) “The Use of a Log-Normal Prior for the Student t-Distribution,” Axioms](https://www.mdpi.com/2075-1680/11/9/462)
