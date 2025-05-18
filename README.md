# Install R package
```r
devtools::install_github("yain22/winratiosim")
```

# Load R package
```r
library(winratiosim)
```

# Implement the power analysis simulation 
## Note: This simulation is optimized for Windows due to the parallel computing setup.
```r
# Install Library
library(plyr)
library(dplyr)
library(parallel)
library(kableExtra)
library(pbmcapply)
library(foreach)
library(doSNOW)

# Simulaiton under power
set.seed(20250518)
power.design_parameters = list(
  nsim = 10000,
  N = 400, # 100, 200, 300, 400, 500
  Randomization.ratio = c(1, 1),
  alpha.JFM = 0, theta.JFM = 1,
  lambda_trt = 0.13, lambda_ctl = 0.15,
  ann.icr_trt = 0.32, ann.icr_ctl = 0.55,
  xbase_trt = 45, xfinal_trt = 45+7.5, sd.delta.x_trt = 20,
  xbase_ctl = 45, xfinal_ctl = 45, sd.delta.x_ctl = 20,
  censorrate_trt = 0.2, censorrate_ctl = 0.2,
  nc = 10
)

## Running simulation
power.sim_res = winratiosim(nsim = power.design_parameters$nsim, N = power.design_parameters$N, Randomization.ratio = power.design_parameters$Randomization.ratio,
                            alpha.JFM = power.design_parameters$alpha.JFM, theta.JFM = power.design_parameters$theta.JFM,
                            lambda_trt = power.design_parameters$lambda_trt, lambda_ctl = power.design_parameters$lambda_ctl,
                            ann.icr_trt = power.design_parameters$ann.icr_trt, ann.icr_ctl = power.design_parameters$ann.icr_ctl,
                            xbase_trt = power.design_parameters$xbase_trt, xfinal_trt = power.design_parameters$xfinal_trt,
                            xbase_ctl = power.design_parameters$xbase_ctl, xfinal_ctl = power.design_parameters$xfinal_ctl,
                            sd.delta.x_trt = power.design_parameters$sd.delta.x_trt, sd.delta.x_ctl = power.design_parameters$sd.delta.x_ctl,
                            censorrate_trt = power.design_parameters$censorrate_trt, censorrate_ctl = power.design_parameters$censorrate_ctl,
                            nc = power.design_parameters$nc)

## FS test
Power_one_sided_superiority_FS_Permutation = mean(power.sim_res$df_FS.analysis.summary$p_value_FS < 0.025, na.rm = T)
Power_binom_CI_one_sided_FS_Permutation = binom.conf.exact(x=sum(power.sim_res$df_FS.analysis.summary$p_value_FS < 0.025, na.rm = T), n=nrow(power.sim_res$df_FS.analysis.summary) - sum(is.na(power.sim_res$df_FS.analysis.summary$p_value_FS)))

## YG test
Power_one_sided_superiority_WR_Ron_Yu = mean(power.sim_res$df_WR.analysis.summary$LB_R_w > 1, na.rm = T)
Power_binom_CI_one_sided_WR_Ron_Yu = binom.conf.exact(x=sum(power.sim_res$df_WR.analysis.summary$LB_R_w > 1, na.rm = T), n=nrow(power.sim_res$df_WR.analysis.summary) - sum(is.na(power.sim_res$df_WR.analysis.summary$LB_R_w)))

# Simulaiton under Type I Error scenario
t1e.design_parameters = list(
  nsim = power.design_parameters$nsim,
  N = power.design_parameters$N,
  Randomization.ratio = power.design_parameters$Randomization.ratio,
  alpha.JFM = power.design_parameters$alpha.JFM, theta.JFM = power.design_parameters$theta.JFM,
  lambda_trt = power.design_parameters$lambda_ctl, lambda_ctl = power.design_parameters$lambda_ctl,
  ann.icr_trt = power.design_parameters$ann.icr_ctl, ann.icr_ctl = power.design_parameters$ann.icr_ctl,
  xbase_trt = power.design_parameters$xbase_ctl, xfinal_trt = power.design_parameters$xfinal_ctl, sd.delta.x_trt = power.design_parameters$sd.delta.x_trt,
  xbase_ctl = power.design_parameters$xbase_ctl, xfinal_ctl = power.design_parameters$xfinal_ctl, sd.delta.x_ctl = power.design_parameters$sd.delta.x_ctl,
  censorrate_trt = power.design_parameters$censorrate_trt, censorrate_ctl = power.design_parameters$censorrate_ctl,
  nc = 10)

# Running simulation
t1e.sim_res = winratiosim(nsim = t1e.design_parameters$nsim, N = t1e.design_parameters$N, Randomization.ratio = t1e.design_parameters$Randomization.ratio,
                          alpha.JFM = t1e.design_parameters$alpha.JFM, theta.JFM = t1e.design_parameters$theta.JFM,
                          lambda_trt = t1e.design_parameters$lambda_trt, lambda_ctl = t1e.design_parameters$lambda_ctl,
                          ann.icr_trt = t1e.design_parameters$ann.icr_trt, ann.icr_ctl = t1e.design_parameters$ann.icr_ctl,
                          xbase_trt = t1e.design_parameters$xbase_trt, xfinal_trt = t1e.design_parameters$xfinal_trt,
                          xbase_ctl = t1e.design_parameters$xbase_ctl, xfinal_ctl = t1e.design_parameters$xfinal_ctl,
                          sd.delta.x_trt = t1e.design_parameters$sd.delta.x_trt, sd.delta.x_ctl = t1e.design_parameters$sd.delta.x_ctl,
                          censorrate_trt = t1e.design_parameters$censorrate_trt, censorrate_ctl = t1e.design_parameters$censorrate_trt,
                          nc = t1e.design_parameters$nc)

## FS test
t1e_one_sided_superiority_FS_Permutation = mean(t1e.sim_res$df_FS.analysis.summary$p_value_FS < 0.025, na.rm = T)
t1e_binom_CI_one_sided_FS_Permutation = binom.conf.exact(x=sum(t1e.sim_res$df_FS.analysis.summary$p_value_FS < 0.025, na.rm = T),n=nrow(t1e.sim_res$df_FS.analysis.summary) - sum(is.na(t1e.sim_res$df_FS.analysis.summary$p_value_FS)))

## YG test
t1e_one_sided_superiority_WR_Ron_Yu = mean(t1e.sim_res$df_WR.analysis.summary$LB_R_w > 1, na.rm = T)
t1e_binom_CI_one_sided_WR_Ron_Yu = binom.conf.exact(x=sum(t1e.sim_res$df_WR.analysis.summary$LB_R_w > 1, na.rm = T),n=nrow(t1e.sim_res$df_WR.analysis.summary) - sum(is.na(t1e.sim_res$df_WR.analysis.summary$LB_R_w)))

# Print out results
## Construct compact Power and Type I Error results
df.power.type1 <- data.frame(
  Method = c("FS test", "YG test"),
  
  Power = paste(
    round(c(Power_binom_CI_one_sided_FS_Permutation[1], Power_binom_CI_one_sided_WR_Ron_Yu[1]), 3),
    "(",
    round(c(Power_binom_CI_one_sided_FS_Permutation[2], Power_binom_CI_one_sided_WR_Ron_Yu[2]), 3),
    ", ",
    round(c(Power_binom_CI_one_sided_FS_Permutation[3], Power_binom_CI_one_sided_WR_Ron_Yu[3]), 3),
    ")",
    sep = ""
  ),
  
  Type_I_Error = paste(
    round(c(t1e_binom_CI_one_sided_FS_Permutation[1], t1e_binom_CI_one_sided_WR_Ron_Yu[1]), 3),
    "(",
    round(c(t1e_binom_CI_one_sided_FS_Permutation[2], t1e_binom_CI_one_sided_WR_Ron_Yu[2]), 3),
    ", ",
    round(c(t1e_binom_CI_one_sided_FS_Permutation[3], t1e_binom_CI_one_sided_WR_Ron_Yu[3]), 3),
    ")",
    sep = ""
  )
)

## Median variance values
df.variance <- data.frame(
  `Median Variance under Power` = c(
    median(power.sim_res$df_WR.analysis.summary$variance_log_R_w_permutation, na.rm = TRUE),
    median(power.sim_res$df_WR.analysis.summary$Var_logR_w, na.rm = TRUE)
  ),
  `Median Variance under Type I Error` = c(
    median(t1e.sim_res$df_WR.analysis.summary$variance_log_R_w_permutation, na.rm = TRUE),
    median(t1e.sim_res$df_WR.analysis.summary$Var_logR_w, na.rm = TRUE)
  )
)

## Round and combine everything
df.combined <- cbind(df.power.type1, round(df.variance, 4))

## Generate the summary string from df.kable (already created)
df_kable_summary <- paste0(
  "Simulation assumptions:\n",
  "  Treatment Arm — Mortality = ", power.design_parameters$lambda_trt,
  ", Annualized HFH = ", power.design_parameters$ann.icr_trt,
  ", ΔKCCQ = ", power.design_parameters$xfinal_trt - power.design_parameters$xbase_trt,
  ", SD(ΔKCCQ) = ", power.design_parameters$sd.delta.x_trt,
  ", Attrition = ", power.design_parameters$censorrate_trt, "\n",
  
  "  Control Arm — Mortality = ", power.design_parameters$lambda_ctl,
  ", Annualized HFH Rate = ", power.design_parameters$ann.icr_ctl,
  ", ΔKCCQ = ", power.design_parameters$xfinal_ctl - power.design_parameters$xbase_ctl,
  ", SD(ΔKCCQ) = ", power.design_parameters$sd.delta.x_ctl,
  ", Attrition = ", power.design_parameters$censorrate_ctl, "\n",
  
  "Randomization Ratio (T:C) = ", power.design_parameters$Randomization.ratio[1], ":", power.design_parameters$Randomization.ratio[2], "\n",
  "Sample Size (N) = ", power.design_parameters$N,
  ", Number of simulations = ", power.design_parameters$nsim, "\n",
  "Joint Frailty Model Parameters: alpha = ", power.design_parameters$alpha.JFM,
  ", theta = ", power.design_parameters$theta.JFM
)

## Final table
df_kable_summary
kable(df.combined, caption = "Simulation Results")
print(paste0("Median value of win ratio is ", round(median(power.sim_res$df_WR.analysis.summary$R_w), 3)))
print(paste0("Median value of probability of tie is ", round(median(power.sim_res$df_Total_probability[, 2]), 3)))
```
