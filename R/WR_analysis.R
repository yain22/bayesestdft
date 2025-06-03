#' Perform Hierarchical Win Ratio and Functional Score Analysis
#'
#' Analyzes matched treatment-control pairs across three prioritized outcome layers using a hierarchical win ratio framework.
#' The function computes win/loss/tie counts at each layer, win ratio test statistics, and functional score statistics as described in
#' Finkelstein and Schoenfeld (1999), Pocock et al. (2012), and Yu & Ganju (2022).
#'
#' @param dataset1 Data frame containing matched pairs and win indicators for the first (highest priority) layer (typically time-to-event).
#' @param dataset2 Data frame for the second layer (e.g., hospitalization or recurrent event score).
#' @param dataset3 Data frame for the third layer (e.g., quality-of-life or patient-reported outcomes).
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{win.losses.count.summary}{A data frame with counts and proportions of wins, ties, and losses by layer and in total.}
#'   \item{sample.size.summary}{Sample sizes for treatment, control, and total number of pairwise comparisons.}
#'   \item{FS.analysis.summary}{Test statistic, variance, z-score, and p-value for the functional score test.}
#'   \item{WR.analysis.summary}{Win ratio, variance estimates, confidence intervals, and p-value for the win ratio test.}
#' }
#'
#' @examples
#' \dontrun{
#' # Assume dataset1, dataset2, dataset3 are prepared as matched pair data sets
#' result <- WR_analysis(dataset1, dataset2, dataset3)
#' result$WR.analysis.summary
#' }
#'
#' @references
#' Finkelstein, D. M., & Schoenfeld, D. A. (1999). Combining mortality and longitudinal measures in clinical trials. \emph{Statistics in Medicine}, 18(11), 1341–1354.  
#' Pocock, S. J., Ariti, C. A., Collier, T. J., & Wang, D. (2012). The win ratio: a new approach to the analysis of composite endpoints in clinical trials based on clinical priorities. \emph{European Heart Journal}, 33(2), 176–182.  
#' Yu, R. X., & Ganju, J. (2022). Sample size formula for a win ratio endpoint. \emph{Statistics in Medicine}, 41(6), 950–963.
#'
#' @export

WR_analysis = function(dataset1,dataset2,dataset3){
  
  N_trt = length(unique(dataset1[dataset1$treatment1 == 1,]$usubjid1)) # Number of patients in treatment group
  N_ctl = length(unique(dataset1[dataset1$treatment1 == 0,]$usubjid1)) # Number of patients in control group
  N = N_trt + N_ctl
  N_comparison_win_ratio = N_trt*N_ctl # Number of comparison for win ratio computation
  dataset1.WR = dataset1[ (dataset1$treatment1 == 1) & (dataset1$treatment2 == 0),] # Data for WR at Layer 1
  dataset2.WR = dataset2[ (dataset2$treatment1 == 1) & (dataset2$treatment2 == 0),] # Data for WR at Layer 2
  dataset3.WR = dataset3[ (dataset3$treatment1 == 1) & (dataset3$treatment2 == 0),] # Data for WR at Layer 3
  
  # Layer 1
  W_T_layer1 = sum(dataset1.WR$score == 1,na.rm = T)   # Number of winner in treatment group at Layer 1
  W_C_layer1 = sum(dataset1.WR$score == -1,na.rm = T)  # Number of winner in control group at Layer 1
  W_0_layer1 = N_comparison_win_ratio - (W_T_layer1 + W_C_layer1) # Number of ties at Layer 1
  
  # Layer 2
  W_T_layer2 = sum(dataset2.WR$score == 1,na.rm = T) - W_T_layer1   # Number of winner in treatment group at Layer 2
  W_C_layer2 = sum(dataset2.WR$score == -1,na.rm = T) - W_C_layer1  # Number of winner in control group at Layer 2
  W_0_layer2 = W_0_layer1 - (W_T_layer2 + W_C_layer2) # Number of ties at Layer 2
  
  # Layer 3
  W_T_layer3 = sum(dataset3.WR$score == 1,na.rm = T) - (W_T_layer1 + W_T_layer2)  # Number of winner in treatment group at Layer 3
  W_C_layer3 = sum(dataset3.WR$score == -1,na.rm = T) - (W_C_layer1 + W_C_layer2) # Number of winner in control group at Layer 3
  W_0_layer3 = W_0_layer2 - (W_T_layer3 + W_C_layer3) # Number of ties at Layer 3
  
  # Total
  W_T = W_T_layer1 + W_T_layer2 + W_T_layer3
  W_C = W_C_layer1 + W_C_layer2 + W_C_layer3
  W_0 = N_comparison_win_ratio - (W_T + W_C) 
  
  # FS test statistics [Finkelstein, D. M., & Schoenfeld, D. A. (1999). Combining mortality and longitudinal measures in clinical trials. Statistics in medicine, 18(11), 1341-1354.]
  {
    T.stat = W_T - W_C
    temp.Ui = dataset3 %>% group_by(usubjid1) %>% summarise(Ui = sum(score, na.rm = T))
    V = ((N_trt*N_ctl)/((N_trt + N_ctl)*(N_trt + N_ctl - 1)))*sum((temp.Ui$Ui)^2, na.rm = T) # Permutation Variance
    z = T.stat/sqrt(V)  
    p_value_FS = pnorm(z, lower.tail = F)
  }
  
  # Win ratio test statistics [Pocock, S. J., Ariti, C. A., Collier, T. J., & Wang, D. (2012). The win ratio: a new approach to the analysis of composite endpoints in clinical trials based on clinical priorities. European heart journal, 33(2), 176-182.]
  {
    R_w = W_T/W_C
    variance_log_R_w_permutation = (log(R_w)/z)^2
    UB_R_w_permutation = exp(log(R_w) + 1.96*sqrt(variance_log_R_w_permutation)) 
    LB_R_w_permutation = exp(log(R_w) - 1.96*sqrt(variance_log_R_w_permutation)) 
  }
  
  # Win ratio test statistics [Yu, R. X., & Ganju, J. (2022). Sample size formula for a win ratio endpoint. Statistics in medicine, 41(6), 950-963.]
  {
    P0 = W_0/N_comparison_win_ratio
    k = N_trt/(N_trt + N_ctl)
    nomi = 4*(1 + P0)                
    denom = 3*k*(1-k)*(1 - P0)
    sigma = sqrt(nomi/denom) 
    variance_log_R_w = (sigma^2)/N
    p_value_R_w = 1 - pnorm(q = log(R_w)*sqrt((3*N_trt*N_ctl*(1-P0))/(4*N*(1+P0)))) 
    UB_R_w = exp(log(R_w) + 1.96*sqrt(variance_log_R_w)) 
    LB_R_w = exp(log(R_w) - 1.96*sqrt(variance_log_R_w)) 
  }
  
  win.losses.count.summary = data.frame(Count = c("Number of wins in Treatment Group","Number of ties","Number of losses in Treatment Group","Sum"),
                                        First_Layer = c(W_T_layer1,W_0_layer1,W_C_layer1, sum(c(W_T_layer1,W_0_layer1,W_C_layer1))),
                                        Second_Layer = c(W_T_layer2,W_0_layer2,W_C_layer2, sum(c(W_T_layer2,W_0_layer2,W_C_layer2))),
                                        Third_Layer = c(W_T_layer3,W_0_layer3,W_C_layer3, sum(c(W_T_layer3,W_0_layer3,W_C_layer3))),
                                        Total_count = c(W_T, W_0, W_C, sum(c(W_T, W_0, W_C))),
                                        Total_probability = round(c(W_T, W_0, W_C, sum(c(W_T, W_0, W_C)))/sum(c(W_T, W_0, W_C)),3) 
  )  
  
  sample.size.summary = data.frame(N = N, N_trt = N_trt, N_ctl = N_ctl,N_comparison_win_ratio = N_comparison_win_ratio)
  FS.analysis.summary = data.frame(T = T.stat, V = V, z = z, p_value_FS = p_value_FS)
  WR.analysis.summary = data.frame(R_w = R_w, logR_w = log(R_w), variance_log_R_w_permutation = variance_log_R_w_permutation,
                                   LB_R_w_95p_permutation = LB_R_w_permutation, UB_R_w_95p_permutation = UB_R_w_permutation, 
                                   Var_logR_w = variance_log_R_w, UB_R_w = UB_R_w, LB_R_w = LB_R_w,
                                   p_value_R_w = p_value_R_w)
  
  res.list = list(win.losses.count.summary = win.losses.count.summary,
                  sample.size.summary = sample.size.summary, 
                  FS.analysis.summary = FS.analysis.summary, 
                  WR.analysis.summary = WR.analysis.summary
  )
  return(res.list)
} 

