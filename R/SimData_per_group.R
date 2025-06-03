#' Simulate Individual-Level Trial Data for a Treatment Group
#'
#' Generates individual-level simulated data for a treatment or control group in a hierarchical win ratio trial.
#' The simulation includes frailty-adjusted time-to-death, recurrent hospitalization events, censoring times,
#' and longitudinal quality-of-life scores (e.g., KCCQ). It outputs a subject-level summary suitable for further analysis.
#'
#' This function is called internally by \code{winratiosim()} and assigns the output as a global variable
#' (\code{surv_1} or \code{surv_0}) based on the treatment group.
#'
#' @param treatment Integer. Treatment group indicator (1 for treatment, 0 for control).
#' @param ngroup Integer. Number of subjects to simulate in this group.
#' @param alpha.JFM Numeric. Alpha parameter for the joint frailty model for mortality.
#' @param theta.JFM Numeric. Theta parameter (frailty variance) for the joint frailty model.
#' @param lambda Numeric. Annual mortality rate.
#' @param ann.icr Numeric. Annual incidence rate of recurrent events (e.g., hospitalizations).
#' @param censorrate Numeric. Annual censoring rate.
#' @param xbase Numeric. Baseline quality-of-life score (e.g., KCCQ).
#' @param xfinal Numeric. Expected final quality-of-life score after follow-up (if no death).
#' @param sd.delta.x Numeric. Standard deviation for change in QoL score.
#'
#' @return This function does not return an object directly. Instead, it assigns a data frame named \code{surv_1}
#' or \code{surv_0} to the global environment depending on the group. The data frame contains the following columns:
#' \describe{
#'   \item{subjid}{Subject ID}
#'   \item{treatment}{Treatment group indicator}
#'   \item{mortd}{Time to death (in days)}
#'   \item{censortime}{Observed censoring or death time}
#'   \item{death}{Death indicator (1 if died before censoring)}
#'   \item{HFH}{Number of recurrent events (e.g., hospitalizations)}
#'   \item{kccq}{Quality-of-life score (missing if death occurred before day 360)}
#' }
#'
#' @details
#' \itemize{
#'   \item Mortality is modeled via a frailty-adjusted exponential distribution.
#'   \item Recurrent events (gap times) are modeled as exponential with frailty.
#'   \item Follow-up is administratively censored at 360 days or lost to follow-up.
#'   \item KCCQ scores are simulated based on a normal distribution bounded between -100 and 100.
#' }
#'
#' @examples
#' \dontrun{
#' SimData_per_group(
#'   treatment = 1, ngroup = 100,
#'   alpha.JFM = 0, theta.JFM = 1,
#'   lambda = 0.13, ann.icr = 0.32,
#'   censorrate = 0.2, xbase = 45, xfinal = 52.5, sd.delta.x = 20
#' )
#' head(surv_1)
#' }
#'
#' @export

SimData_per_group <- function(treatment, ngroup, alpha.JFM, theta.JFM, lambda, ann.icr,
                              censorrate, xbase, xfinal, sd.delta.x) {
  
  df_enroll_follow <- data.frame(group = treatment, followup = 12, subjid = 1:ngroup)
  n <- ngroup
  
  # 1. Generate HFH Gap Time
  frailty <- rgamma(n, shape = 1 / theta.JFM, rate = 1 / theta.JFM)
  df_hfh_gap <- as.data.frame(matrix(rep(df_enroll_follow$subjid, each = 100), byrow = FALSE, ncol = 1))
  colnames(df_hfh_gap) <- "subjid"
  df_hfh_gap$treatment <- treatment
  df_hfh_gap$frailty <- rep(frailty, each = 100) 
  betaHF <- ann.icr / 360 
  rate_hfh <- df_hfh_gap$frailty * betaHF
  df_hfh_gap$z <- rexp(nrow(df_hfh_gap), rate = rate_hfh) 
  df_hfh_gap[] <- sapply(df_hfh_gap[], as.numeric)
  df_hfh_gap <- df_hfh_gap[, c("subjid", "treatment", "z")]
  df_hfh_gap <- ddply(df_hfh_gap, .(subjid), transform,
                      treatment = treatment,
                      z = z,
                      startm = cumsum(c(0, z[1:(length(z) - 1)])),
                      stopm = cumsum(z))
  
  # 2. Generate time-to-all-cause death
  df_mortality <- as.data.frame(matrix(c(df_enroll_follow$subjid), byrow = FALSE, ncol = 1))
  df_mortality <- as.data.frame(df_mortality[order(df_mortality$V1), ])
  colnames(df_mortality) <- "subjid"
  df_mortality$treatment <- treatment
  df_mortality$frailty <- frailty ^ alpha.JFM 
  betaMor <- -log(1 - lambda) / 360 
  df_mortality$r <- df_mortality$frailty * betaMor
  df_mortality$mortd <- rexp(n, rate = df_mortality$r) 
  df_mortality[] <- sapply(df_mortality[], as.numeric)
  df_mortality <- df_mortality[order(df_mortality$subjid, df_mortality$mortd), ]
  
  # 3. Generate time to lost-to-follow-up
  df_censoring <- as.data.frame(matrix(c(df_enroll_follow$subjid), byrow = FALSE, ncol = 1))
  df_censoring <- as.data.frame(df_censoring[order(df_censoring$V1), ])
  colnames(df_censoring) <- "subjid"
  df_censoring$ftime <- 360
  cr <- -log(1 - censorrate) / 360
  df_censoring$ctime <- rexp(n, rate = cr) 
  df_censoring[] <- sapply(df_censoring[], as.numeric)
  df_censoring <- df_censoring[order(df_censoring$subjid, df_censoring$ctime), ]
  
  # 4. Combine Mortality and HFH
  df_combined <- full_join(df_mortality, df_hfh_gap, by = c("subjid", "treatment"))
  df_combined <- full_join(df_combined, df_enroll_follow, by = "subjid")
  df_combined <- full_join(df_combined, df_censoring, by = "subjid")
  
  df_event <- df_combined
  df_event$followdays <- df_event$followup * 30
  df_event$fudays <- if_else(df_event$followdays <= df_event$ctime, df_event$followdays, df_event$ctime) 
  df_event$censortime <- pmin(df_event$fudays, df_event$mortd) 
  df_event$finalstopm <- if_else(df_event$stopm > df_event$censortime, df_event$censortime, df_event$stopm)
  df_event$censor_new <- if_else(df_event$stopm > df_event$censortime, 0, 1)
  
  df_event$deathdays <- if_else(df_event$mortd > df_event$fudays, floor(df_event$fudays), floor(df_event$mortd))
  df_event$death <- if_else(df_event$mortd > df_event$fudays, 0, 1)
  df_event <- df_event[!df_event$startm > df_event$finalstopm, ]
  df_event <- df_event[order(df_event$subjid), ]
  
  df_hfh_summary <- ddply(df_event, .(subjid), summarize, HFH = sum(censor_new))
  df_final <- full_join(df_event, df_hfh_summary, by = "subjid")
  df_final <- df_final[!duplicated(df_final$subjid, fromLast = TRUE), ]
  df_final <- df_final[, c("subjid", "treatment", "mortd", "followup", "group",
                           "ftime", "ctime", "followdays", "fudays",
                           "censortime", "deathdays", "death", "HFH")]
  
  # 5. Add KCCQ Data 
  df_kccq <- df_final
  mean_x_change <- xfinal - xbase
  kccq_change <- pmax(pmin(rnorm(n, mean = mean_x_change, sd = sd.delta.x), 100), -100)
  df_kccq$kccq <- if_else(df_kccq$deathdays < 360, NA_real_, kccq_change)
  
  # Final output
  surv <- df_kccq
  
  if (treatment == 1) {
    assign("surv_1", surv, envir = .GlobalEnv)
  } else if (treatment == 0) {
    assign("surv_0", surv, envir = .GlobalEnv)
  }
}
