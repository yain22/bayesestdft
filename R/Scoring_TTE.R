#' Score Time-to-Event Pairwise Comparisons for Win Ratio
#'
#' This function assigns win/loss/tie scores to matched subject pairs based on their time-to-event outcomes
#' (e.g., death or other hard endpoints). It is designed to be used as the first layer in a hierarchical win ratio framework.
#' The scoring accounts for censoring and ties due to simultaneous or unknown outcomes.
#'
#' @param dataset A data frame containing matched pairs of subjects with outcome and censoring indicators.
#' @param var1 Character. Name of the column representing the time-to-event variable for the subject in position 1 (e.g., \code{"deathdays1"}).
#' @param var2 Character. Name of the column representing the time-to-event variable for the subject in position 2 (e.g., \code{"deathdays2"}).
#' @param censor1 Character. Name of the censoring indicator variable (1 = event occurred, 0 = censored) for subject 1.
#' @param censor2 Character. Name of the censoring indicator variable for subject 2.
#'
#' @return A modified data frame identical to the input \code{dataset} but with updated values for:
#' \describe{
#'   \item{score}{Win/loss/tie indicator: 1 if subject 1 wins, -1 if subject 2 wins, NA for tie or unresolvable comparison.}
#'   \item{WR_cat}{Character column indicating whether subject 1 is a "winner" or "loser" due to the time-to-event outcome.}
#' }
#'
#' @details
#' \itemize{
#'   \item A win is declared if subject 1 survives significantly longer than subject 2.
#'   \item A loss is declared if subject 2 survives significantly longer than subject 1.
#'   \item Ties or unresolved scores occur if censoring interferes with outcome comparison or both are censored.
#' }
#'
#' @examples
#' \dontrun{
#' df_scored <- Scoring_TTE(
#'   dataset = df_FS_input,
#'   var1 = "deathdays1", var2 = "deathdays2",
#'   censor1 = "death1", censor2 = "death2"
#' )
#' head(df_scored$score)
#' }
#'
#' @export

Scoring_TTE <- function (dataset, var1, var2, censor1, censor2){
  
  temp <- dataset[dataset$score %in% c(NA_real_, 0), ]
  rest <- dataset[!dataset$score %in% c(NA_real_, 0), ]
  
  # The event happened later - Won/Lose decided
  temp$score <- if_else(temp[, var1] - temp[, var2] > 10^(-10), 1, temp$score)
  temp$score <- if_else(temp[, var2] - temp[, var1] > 10^(-10), -1, temp$score)
  
  # One of the patients had censor - The event happened after the censor - Tie
  temp$score <- if_else(temp[, censor1] %in% 1 & temp[, censor2] %in% 0 & temp[, var1] > temp[, var2], NA_real_, temp$score)
  temp$score <- if_else(temp[, censor1] %in% 0 & temp[, censor2] %in% 1 & temp[, var1] < temp[, var2], NA_real_, temp$score)
  
  # Both patients have censors - Tie
  temp$score <- if_else(temp[, censor1] %in% 0 & temp[, censor2] %in% 0, NA_real_, temp$score)
  
  # Add reason for win ratio results
  temp$WR_cat <- if_else(temp$score %in% 1, paste(gsub("[[:digit:]]+", "", var1), "winner", sep = " "), temp$WR_cat)
  temp$WR_cat <- if_else(temp$score %in% -1, paste(gsub("[[:digit:]]+", "", var1), "loser", sep = " "), temp$WR_cat)
  
  temp <- rbind(temp, rest)
  temp <- temp[order(temp$usubjid1, temp$usubjid2), ]
  
  return(as.data.frame(temp))
}

  
  
