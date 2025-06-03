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
