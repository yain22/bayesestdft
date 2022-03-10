# bayesestdft: An R package for Bayesian estimation of the number of degrees of the freedom of the Student's t-distribution

## Contents
* [Overview](#overview)
* [Installation](#installation)
* [Example 1](#example-1)

## Overview
An R package `bayesestdft` includes tools to implement Bayesian estimation of the number of degrees of the freedom of the Student's t-distribution, developed by [Dr. Se Yoon Lee](https://sites.google.com/view/seyoonlee) (seyoonlee.stat.math@gmail.com). At current version, the main functions are `BayesLNP` and `BayesJeffreys` that implement Markov Chain Monte Carlo algorithms to sample from the number of degrees of the freedom of the Student's t-distribution. To operatre the function `BayesJeffreys`, user needs to install `R library(numDeriv)`.

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

## Example 1

```r
x = rt(n = 100, df = 0.1)
nu = BayesLNP(x)
mean(nu)
```

## Example 2

```r
library(dplyr)
data(index_return)
index_return_US <- filter(index_return, Country == "United States")
x = index_return_US$log_return_rate
nu = BayesLNP(x)
mean(nu)
```

## Example 3

```r
x = rt(n = 100, df = 0.1)
nu1 = BayesJeffreys(x, sampling.alg = "MH")
nu2 = BayesJeffreys(x, sampling.alg = "MALA")
mean(nu1)
mean(nu2)
```

## Example 4

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
