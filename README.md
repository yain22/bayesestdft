# bayesestdft: An R package for Bayesian estimation of the number of degrees of the freedom of the Student's t-distribution

## Contents
* [Overview](#overview)
* [Installation](#installation)
* [Example](#example)

## Overview
An R package `bayesestdft` includes tools to implement Bayesian estimation of the number of degrees of the freedom of the Student's t-distribution, developed by [Dr. Se Yoon Lee](https://sites.google.com/view/seyoonlee) (seyoonlee.stat.math@gmail.com). At current version, the main functions are `BayesLNP` and `BayesJeffreys` that implement Markov Chain Monte Carlo algorithms to sample from the number of degrees of the freedom of the Student's t-distribution. 

## Required R version
```r
R version 4.0.4 (or higher)
```

## Installation

```r
library(devtools)
devtools::install_github("yain22/bayesestdft")
ibrary(bayesestdft)
```

## Example 1

```r
# Example 1 (Simulation study):
x = rt(n = 100, df = 0.1)
nu = BayesLNP(x)
mean(nu)
```
