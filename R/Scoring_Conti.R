#' Score Continuous Outcomes for Win Ratio Comparison
#'
#' Assigns win/loss/tie scores to matched subject pairs based on continuous outcomes, such as hospitalization rate or quality-of-life score.
#' This function is typically used in the second or third layer of a hierarchical win ratio framework,
#' following initial scoring based on time-to-event outcomes.
#'
#' @param dataset A data frame containing matched subject pairs and prior win/loss indicators (i.e., column \code{score}).
#' @param higher_better Character. Set to \code{"Yes"} if higher values of the outcome indicate better health (e.g., QoL),
#' or \code{"No"} if lower values are preferred (e.g., hospitalization count).
#' @param var1 Character. Name of the variable for subject 1 (e.g., \code{"kccq1"} or \code{"HFH_Annual1"}).
#' @param var2 Character. Name of the variable for subject 2 (e.g., \code{"kccq2"} or \code{"HFH_Annual2"}).
#'
#' @return A modified version of the input \code{dataset}, where:
#' \describe{
#'   \item{score}{Updated to 1 (subject 1 wins), -1 (subject 2 wins), 0 (tie), or \code{NA} (missing values prevent comparison).}
#'   \item{WR_cat}{Updated win/loss explanation at this layer (e.g., "kccq winner").}
#' }
#'
#' @details
#' The function only assigns scores to pairs where the current \code{score} is \code{NA} or 0,
#' indicating that prior layers (e.g., death) did not resolve the comparison.
#'
#' A small tolerance (\eqn{10^{-10}}) is used to define ties for nearly equal values.
#'
#' @examples
#' \dontrun{
#' df_scored <- Scoring_Conti(
#'   dataset = df_after_TTE,
#'   higher_better = "Yes",
#'   var1 = "kccq1", var2 = "kccq2"
#' )
#' table(df_scored$score)
#' }
#'
#' @export

Scoring_Conti <- function (dataset, higher_better, var1, var2){
  
  temp <- dataset[dataset$score %in% c(NA_real_, 0), ]
  rest <- dataset[!dataset$score %in% c(NA_real_, 0), ]
  
  # numeric & the higher, the better
  if (higher_better == "Yes"){
    
    temp$score = if_else(temp[, var1] - temp[, var2] >= 10^(-10), 1, temp$score)
    temp$score = if_else(abs(temp[, var1] - temp[, var2]) < 10^(-10), 0, temp$score)
    temp$score = if_else(temp[, var2] - temp[, var1] >= 10^(-10), -1, temp$score)
    temp$score = if_else(temp[, var1] %in% NA_real_ | temp[, var2] %in% NA_real_, NA_real_, temp$score)
    
    # numeric & the lower, the better
  } else if (higher_better == "No"){
    temp$score = if_else(temp[, var2] - temp[, var1] >= 10^(-10), 1, temp$score)
    temp$score = if_else(abs(temp[, var1] - temp[, var2]) < 10^(-10), 0, temp$score)
    temp$score = if_else(temp[, var1] - temp[, var2] >= 10^(-10), -1, temp$score)
    temp$score = if_else(temp[, var1] %in% NA_real_ | temp[, var2] %in% NA_real_, NA_real_, temp$score)
  }
  
  # Add reason for win ratio results
  temp$WR_cat = if_else(temp$WR_cat %in% "" & temp$score %in% 1, paste(gsub('[[:digit:]]+', '', var1), "winner", sep = " "), temp$WR_cat)
  temp$WR_cat = if_else(temp$WR_cat %in% "" & temp$score %in% -1, paste(gsub('[[:digit:]]+', '', var1), "loser", sep = " "), temp$WR_cat)
  
  temp = rbind(temp, rest)
  temp = temp[order(temp$usubjid1, temp$usubjid2), ]
  
  return(as.data.frame(temp))
}
