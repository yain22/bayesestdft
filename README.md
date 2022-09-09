# bayesestdft: An R package for Bayesian estimation of the number of degrees of the freedom of the Student's t-distribution

## Contents
* [Overview](#overview)
* [Installation](#installation)
* [Example 1](#example-1)
* [Example 2](#example-2)
* [Example 3](#example-3)
* [Example 4](#example-4)


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

### Example 1

```r
x = rt(n = 100, df = 0.1)
nu = BayesLNP(x)
mean(nu)
```

### Example 2

```r
library(dplyr)
data(index_return)
index_return_US <- filter(index_return, Country == "United States")
x = index_return_US$log_return_rate
nu = BayesLNP(x)
mean(nu)
```

### Example 3

```r
x = rt(n = 100, df = 0.1)
nu1 = BayesJeffreys(x, sampling.alg = "MH")
nu2 = BayesJeffreys(x, sampling.alg = "MALA")
mean(nu1)
mean(nu2)
```

### Example 4

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
