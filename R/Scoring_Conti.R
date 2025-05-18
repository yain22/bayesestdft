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
