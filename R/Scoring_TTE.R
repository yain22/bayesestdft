Scoring_TTE <- function (dataset, var1, var2, censor1, censor2){
  
  temp <- dataset[dataset$score %in% c(NA_real_, 0), ]
  rest <- dataset[!dataset$score %in% c(NA_real_, 0), ]
  
  # The event happened later - Won/Lose decided
  temp$score <- dplyr::if_else(temp[, var1] - temp[, var2] > 10^(-10), 1, temp$score)
  temp$score <- dplyr::if_else(temp[, var2] - temp[, var1] > 10^(-10), -1, temp$score)
  
  # One of the patients had censor - The event happened after the censor - Tie
  temp$score <- dplyr::if_else(temp[, censor1] %in% 1 & temp[, censor2] %in% 0 & temp[, var1] > temp[, var2], NA_real_, temp$score)
  temp$score <- dplyr::if_else(temp[, censor1] %in% 0 & temp[, censor2] %in% 1 & temp[, var1] < temp[, var2], NA_real_, temp$score)
  
  # Both patients have censors - Tie
  temp$score <- dplyr::if_else(temp[, censor1] %in% 0 & temp[, censor2] %in% 0, NA_real_, temp$score)
  
  # Add reason for win ratio results
  temp$WR_cat <- dplyr::if_else(temp$score %in% 1, paste(gsub("[[:digit:]]+", "", var1), "winner", sep = " "), temp$WR_cat)
  temp$WR_cat <- dplyr::if_else(temp$score %in% -1, paste(gsub("[[:digit:]]+", "", var1), "loser", sep = " "), temp$WR_cat)
  
  temp <- rbind(temp, rest)
  temp <- temp[order(temp$usubjid1, temp$usubjid2), ]
  
  return(as.data.frame(temp))
}

  
  