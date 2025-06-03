#' Simulate Hierarchical Win Ratio Trials
#'
#' Simulates replicated trials using a hierarchical win ratio analysis framework with three prioritized outcome layers:
#' time-to-event (e.g., death), a continuous hospitalization frequency variable, and a continuous quality-of-life score.
#' The simulation incorporates randomization, joint frailty model parameters, censoring, and covariate progression over time.
#' Parallel computation is supported to speed up large-scale simulations.
#'
#' @param nsim Integer. Number of simulated trials to run.
#' @param N Integer. Total number of subjects in a trial (across both arms).
#' @param Randomization.ratio Numeric vector of length 2 indicating the treatment-to-control allocation ratio (e.g., \code{c(1, 1)} for 1:1 randomization).
#' @param alpha.JFM Numeric. Alpha parameter for the joint frailty model.
#' @param theta.JFM Numeric. Theta parameter for the joint frailty model.
#' @param lambda_trt Numeric. Event rate for the treatment group.
#' @param lambda_ctl Numeric. Event rate for the control group.
#' @param ann.icr_trt Numeric. Annual intercurrent event rate for the treatment group.
#' @param ann.icr_ctl Numeric. Annual intercurrent event rate for the control group.
#' @param xbase_trt Numeric. Baseline value of the continuous outcome for the treatment group.
#' @param xfinal_trt Numeric. Final value of the continuous outcome for the treatment group.
#' @param xbase_ctl Numeric. Baseline value of the continuous outcome for the control group.
#' @param xfinal_ctl Numeric. Final value of the continuous outcome for the control group.
#' @param sd.delta.x_trt Numeric. Standard deviation for the change in the continuous outcome in the treatment group.
#' @param sd.delta.x_ctl Numeric. Standard deviation for the change in the continuous outcome in the control group.
#' @param censorrate_trt Numeric. Censoring rate for the treatment group.
#' @param censorrate_ctl Numeric. Censoring rate for the control group.
#' @param nc Integer. Number of CPU cores to use for parallel processing. Default is 1.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{df_FS.analysis.summary}{A data frame summarizing functional score win results across simulations.}
#'   \item{df_WR.analysis.summary}{A data frame summarizing win ratio results across simulations.}
#'   \item{df_sample.size.summary}{A data frame with the sample sizes used in each arm per simulation.}
#'   \item{df_Total_probability}{A data frame showing average win/tie/loss probabilities across all simulations.}
#'   \item{df_Total_count}{A data frame showing average win/tie/loss counts across all simulations.}
#' }
#'
#' @examples
#' \dontrun{
#' power.design_parameters <- list(
#'   nsim = 10000
#'   N = 400,
#'   Randomization.ratio = c(1, 1),
#'   alpha.JFM = 0, theta.JFM = 1,
#'   lambda_trt = 0.13, lambda_ctl = 0.15,
#'   ann.icr_trt = 0.32, ann.icr_ctl = 0.55,
#'   xbase_trt = 45, xfinal_trt = 45 + 7.5, sd.delta.x_trt = 20,
#'   xbase_ctl = 45, xfinal_ctl = 45, sd.delta.x_ctl = 20,
#'   censorrate_trt = 0.2, censorrate_ctl = 0.2,
#'   nc = 1
#' )
#'
#' result <- do.call(winratiosim, power.design_parameters)
#' result$df_WR.analysis.summary
#' }
#'
#' @references
#' Lee, S.Y. (2025). A note on the sample size formula for a win ratio endpoint. \emph{Statistics in Medicine}.
#'
#' @export

winratiosim <- function(nsim, N, Randomization.ratio, alpha.JFM, theta.JFM,
                        lambda_trt, lambda_ctl, ann.icr_trt, ann.icr_ctl,
                        xbase_trt, xfinal_trt, xbase_ctl, xfinal_ctl,
                        sd.delta.x_trt, sd.delta.x_ctl,
                        censorrate_trt, censorrate_ctl, nc = 1) {
  
  # Replication function
  simulate_one_trial <- function(r) {
    set.seed(r*1234)
    
    allocation_ratio <- Randomization.ratio[1] / sum(Randomization.ratio)
    n_treatment <- rbinom(n = 1, size = N, prob = allocation_ratio)
    n_control <- N - n_treatment
    
    # Simulate treatment arm
    SimData_per_group(
      treatment = 1, ngroup = n_treatment,
      alpha.JFM = alpha.JFM, theta.JFM = theta.JFM,
      ann.icr = ann.icr_trt, lambda = lambda_trt, censorrate = censorrate_trt,
      xbase = xbase_trt, xfinal = xfinal_trt, sd.delta.x = sd.delta.x_trt
    )
    
    # Simulate control arm
    SimData_per_group(
      treatment = 0, ngroup = n_control,
      alpha.JFM = alpha.JFM, theta.JFM = theta.JFM,
      ann.icr = ann.icr_ctl, lambda = lambda_ctl, censorrate = censorrate_ctl,
      xbase = xbase_ctl, xfinal = xfinal_ctl, sd.delta.x = sd.delta.x_ctl
    )
    
    # Combine trial data
    df_trial <- rbind(surv_1, surv_0)
    df_trial$subjid <- if_else(df_trial$treatment == 0, df_trial$subjid + 1000, df_trial$subjid)
    df_trial$HFH_Annual <- (df_trial$HFH / df_trial$censortime) * 360
    colnames(df_trial)[1] <- "usubjid"
    df_trial <- df_trial[order(df_trial$usubjid), ]
    
    # FS score data preparation
    df_base <- df_trial[rep(seq_len(nrow(df_trial)), each = nrow(df_trial)), ]
    colnames(df_base) <- paste0(colnames(df_base), "1")
    
    df_compare <- df_trial[rep(seq_len(nrow(df_trial)), times = nrow(df_trial)), ]
    colnames(df_compare) <- paste0(colnames(df_compare), "2")
    
    df_FS_input <- cbind(df_base, df_compare)
    df_FS_input$score <- NA_real_
    df_FS_input$WR_cat <- ""
    
    # 3-layer win framework
    df_outcome1 <- Scoring_TTE(
      dataset = df_FS_input,
      var1 = "deathdays1", var2 = "deathdays2",
      censor1 = "death1", censor2 = "death2"
    )
    
    df_outcome2 <- Scoring_Conti(
      dataset = df_outcome1,
      higher_better = "No",
      var1 = "HFH_Annual1", var2 = "HFH_Annual2"
    )
    
    df_outcome3 <- Scoring_Conti(
      dataset = df_outcome2,
      higher_better = "Yes",
      var1 = "kccq1", var2 = "kccq2"
    )
    
    # Layered datasets
    df_layer1 <- df_outcome1 %>% dplyr::select(usubjid1, treatment1, usubjid2, treatment2, score)
    df_layer2 <- df_outcome2 %>% dplyr::select(usubjid1, treatment1, usubjid2, treatment2, score)
    df_layer3 <- df_outcome3 %>% dplyr::select(usubjid1, treatment1, usubjid2, treatment2, score)
    
    # Win ratio analysis
    df_result <- WR_analysis(
      dataset1 = df_layer1,
      dataset2 = df_layer2,
      dataset3 = df_layer3
    )
    
    return(df_result)
  }
  
  # For MAC User
  #res_par_com = pbmclapply(X = 1:nsim, FUN = rep_fcn, mc.cores = nc)
  
  # Parallel Setup (Window)
  cl <- makeCluster(rep("localhost", nc))
  clusterExport(cl, c("SimData_per_group", "Scoring_TTE", "Scoring_Conti", "WR_analysis"))
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  opts <- list(progress = function(n) setTxtProgressBar(pb, n))
  
  sim_results <- foreach(r = 1:nsim, .options.snow = opts, .packages = c("dplyr", "plyr")) %dopar% {
    simulate_one_trial(r)
  }
  
  close(pb)
  stopCluster(cl)
  
  # Aggregate Results
  df_FS_summary <- do.call(rbind, lapply(sim_results, function(x) x$FS.analysis.summary))
  df_WR_summary <- do.call(rbind, lapply(sim_results, function(x) x$WR.analysis.summary))
  df_sample_size <- do.call(rbind, lapply(sim_results, function(x) x$sample.size.summary))
  
  df_prob_summary <- do.call(rbind, lapply(sim_results, function(x) x$win.losses.count.summary$Total_probability))
  df_count_summary <- do.call(rbind, lapply(sim_results, function(x) x$win.losses.count.summary$Total_count))
  
  # Set names after combining
  colnames(df_prob_summary) <- c("Prob_of_Win_trt", "Prob_of_tie", "Prob_of_Win_ctl", "ALL")
  colnames(df_count_summary) <- c("Prob_of_Win_trt", "Prob_of_tie", "Prob_of_Win_ctl", "ALL")
  
  # Return as named list
  result_list <- list(
    df_FS.analysis.summary = df_FS_summary,
    df_WR.analysis.summary = df_WR_summary,
    df_sample.size.summary = df_sample_size,
    df_Total_probability = df_prob_summary,
    df_Total_count = df_count_summary
  )
  
  return(result_list)
}
