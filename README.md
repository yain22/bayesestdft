# bayesestdft: An R package for Bayesian estimation of the number of degrees of the freedom of the Student's t-distribution

## Contents
* [Overview](#overview)
* [Installation](#installation)
* [Example](#example)
* [References](#References)

## Overview
An R package `bayesestdft` includes tools to implement Bayesian estimation of the number of degrees of the freedom of the Student's t-distribution, developed by [Dr. Se Yoon Lee](https://sites.google.com/view/seyoonlee) (seyoonlee.stat.math@gmail.com). 

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

## Example

R package `bayesestdft` contains COVID-19 dataset comprises (i) `time_series_data` and (ii) `design_matrix` to train the BHRM. `time_series_data` includes infection growth curve from January 22nd to May 14th for 40 global countries and `design_matrix` has 45 predictors that cover health care resources, population statistics, disease prevalence, etc. Figure 2 displays infection trajectories for eight countries (US, Russia, UK, Brazil, Germany, China, India, and South Korea), spanning from January 22nd to May 14th, which accounts for 114 days.

***Figure 2: Infection trajectories for eight countries updated on May 14th (Data source: JHU CSSE).***
![](https://github.com/StevenBoys/BHRM/blob/main/Image/infect_COVID-19.png?raw=true)


To load the COVID-19 dataset, use
```r
library(BHRM)
# load the data
data("design_matrix")
data("time_series_data")
Y = time_series_data[, -c(1:2)]; X = design_matrix[, -c(1:2)]
```
To train BHRM with the aforementioned dataset, use the R function [`BHRM_cov`](https://github.com/StevenBoys/BHRM/blob/main/R/BHRM_cov.R) as following way. The [`BHRM_cov`](https://github.com/StevenBoys/BHRM/blob/main/R/BHRM_cov.R) implements a Gibbs sampling algorithm to sample from the posterior distribution for the BHRM given the dataset. It may need at least 10 minutes to train the data, depending on CPU speed.
```r
# set the hyperparameters
seed.no = 1 ; burn = 5000 ; nmc = 5000 ; thin = 30; varrho = 0
pro.var.theta.2 = 0.0002 ; pro.var.theta.3 = 0.05; mu = 0 ; rho.sq = 1
t.values = list(); num_days = 14
for(i in 1:nrow(Y)){
  t.values[[i]] = c(1:(ncol(Y) - num_days))
}
Y = Y[, c(1:(ncol(Y) - num_days))]
# run the model
res_cov = BHRM_cov(Y = Y, X = X, t.values = t.values, seed.no = seed.no, burn = burn,   
                   nmc = nmc, thin = thin, varrho = varrho, pro.var.theta.2 = pro.var.theta.2, 
                   pro.var.theta.3 = pro.var.theta.3, mu = mu, rho.sq = rho.sq)  
```
## References

[1] [Se Yoon Lee, Bowen Lei, and Bani K. Mallick. (2020) “Estimation of COVID19 spread curves integrating global data and borrowing information,” PLOS ONE](https://journals.plos.org/plosone/article/authors?id=10.1371/journal.pone.0236860)

[2] [Se Yoon Lee and Bani K. Mallick. (2021) “Bayesian Hierarchical modeling: application towards production results in the Eagle Ford Shale of South Texas,” Sankhyā: The Indian Journal of Statistics, Series B](https://rdcu.be/ceg4p) ; [[Github]](https://github.com/yain22/SWM)

[3] [Davidian, M., and Giltinan, D. M. (1995). Nonlinear models for repeated measurement data (Vol. 62). CRC press.](https://books.google.com/books?hl=en&lr=&id=0eSIBPAL4qsC&oi=fnd&pg=IA7&dq=nonlinear+mixed+effect+model+giltnan&ots=9frDPH3F4J&sig=L5Wz91waGu447OdyYHQ8Vp5ckQc#v=onepage&q=nonlinear%20mixed%20effect%20model%20giltnan&f=false)
