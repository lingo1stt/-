#################################################################
#genetic/physical distance relationship with cross validation
#################################################################
library(dplyr)
library(tidyverse)
library(caret)
##read in data
total <- data.frame(matrix(NA,ncol = 6,nrow = 1))
colnames(total) <- c("set","map","mkr","phys","gen","vld")
colnames(total) <- col
for (i in 1:12) {
  data <- read.csv(paste0('C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation_hotspot_result/','chr',i,'/local/genetic_distance_chr',i,'.csv'))
  #選取基數行
  data <- data[seq(1,nrow(data),by = 2),]
  assign(paste0("gen_distance_",i),data,envir = .GlobalEnv)
  total <- rbind(total,data)
}

###############################################################################12染色體分別選取一條當作train，其餘當作test dataset. 以"三次方"建模
CV_result <- matrix(NA,nrow = 12,ncol = 3)
for (i in 1:12) {
  ###create test data and train data
  train_index <- which(total$map == paste0("chromosome ",i))
  train_data <- total[train_index,]
  test_data <- total[-train_index,]
  #model construction
  model <- lm(gen~poly(phys,3), data = train_data)
  predictions <- model %>% predict(test_data) 
  actual <- as.vector(test_data$gen)
  ##create dataframe for actual value and predict value
  p <- data.frame(p = predictions, a = actual)
  p <- na.omit(p)
  
  #calculate the quality index for the model
  ## R squared
  R2 <- round(R2(actual,predictions,na.rm = T),3)
  ##adj R squared
  n=nrow(p)
  pp <- 2
  R_squared_adj = round( 1 - (1 - R2) * ((n - 1)/(n-pp-1)),3)
  ##RMSE
  RMSE <- RMSE(actual,predictions, na.rm = T)
  ###interprete the result into matrix
  CV_result[i,1:3] <- c(RMSE, R2,R_squared_adj )
}
CV_result <- as.data.frame(CV_result)
colnames(CV_result) <- c("RMSE","R squared", "adj-R squared ")
apply(CV_result,2, function(x){round(mean(x),3)})


############################################################################12染色體分別選取一條當作train，其餘當作test dataset. 以"二次方"建模
CV_result_2 <- matrix(NA,nrow = 12,ncol = 3)
for (i in 1:12) {
  ###create test data and train data
  train_index <- which(total$map == paste0("chromosome ",i))
  train_data <- total[train_index,]
  test_data <- total[-train_index,]
  #model construction
  model <- lm(gen~poly(phys,2), data = train_data)
  predictions <- model %>% predict(test_data) 
  actual <- as.vector(test_data$gen)
  ##create dataframe for actual value and predict value
  p <- data.frame(p = predictions, a = actual)
  p <- na.omit(p)
  
  #calculate the quality index for the model
  ## R squared
  R2 <- round(R2(actual,predictions,na.rm = T),3)
  ##adj R squared
  n=nrow(p)
  pp <- 2
  R_squared_adj = round( 1 - (1 - R2) * ((n - 1)/(n-pp-1)),3)
  ##RMSE
  RMSE <- RMSE(actual,predictions, na.rm = T)
  ###interprete the result into matrix
  CV_result_2[i,1:3] <- c(RMSE, R2,R_squared_adj )
}
CV_result_2 <- as.data.frame(CV_result_2)
colnames(CV_result_2) <- c("RMSE","R squared", "adj-R squared ")
apply(CV_result_2,2, function(x){round(mean(x),3)})

################################################################################12染色體分別選取一條當作train，其餘當作test dataset. 以"一次方"建模
CV_result_1 <- matrix(NA,nrow = 12,ncol = 3)
for (i in 1:12) {
  ###create test data and train data
  train_index <- which(total$map == paste0("chromosome ",i))
  train_data <- total[train_index,]
  test_data <- total[-train_index,]
  #model construction
  model <- lm(gen~phys, data = train_data)
  predictions <- model %>% predict(test_data) 
  actual <- as.vector(test_data$gen)
  ##create dataframe for actual value and predict value
  p <- data.frame(p = predictions, a = actual)
  p <- na.omit(p)
  
  #calculate the quality index for the model
  ## R squared
  R2 <- round(R2(actual,predictions,na.rm = T),3)
  ##adj R squared
  n=nrow(p)
  pp <- 2
  R_squared_adj = round( 1 - (1 - R2) * ((n - 1)/(n-pp-1)),3)
  ##RMSE
  RMSE <- RMSE(actual,predictions, na.rm = T)
  ###interprete the result into matrix
  CV_result_1[i,1:3] <- c(RMSE, R2,R_squared_adj )
}
CV_result_1 <- as.data.frame(CV_result_1)
colnames(CV_result_1) <- c("RMSE","R squared", "adj-R squared ")
apply(CV_result_1,2, function(x){round(mean(x),3)})

################################################################################12染色體分別選取一條當作train，其餘當作test dataset. 以"四次方"建模
CV_result_4 <- matrix(NA,nrow = 12,ncol = 3)
for (i in 1:12) {
  ###create test data and train data
  train_index <- which(total$map == paste0("chromosome ",i))
  train_data <- total[train_index,]
  test_data <- total[-train_index,]
  #model construction
  model <- lm(gen~poly(phys,4), data = train_data)
  predictions <- model %>% predict(test_data) 
  actual <- as.vector(test_data$gen)
  ##create dataframe for actual value and predict value
  p <- data.frame(p = predictions, a = actual)
  p <- na.omit(p)
  
  #calculate the quality index for the model
  ## R squared
  R2 <- round(R2(actual,predictions,na.rm = T),3)
  ##adj R squared
  n=nrow(p)
  pp <- 2
  R_squared_adj = round( 1 - (1 - R2) * ((n - 1)/(n-pp-1)),3)
  ##RMSE
  RMSE <- RMSE(actual,predictions, na.rm = T)
  ###interprete the result into matrix
  CV_result_4[i,1:3] <- c(RMSE, R2,R_squared_adj )
}
CV_result_4 <- as.data.frame(CV_result_4)
colnames(CV_result_4) <- c("RMSE","R squared", "adj-R squared ")
apply(CV_result_4,2, function(x){round(mean(x),3)})



##################儲存檔案
write.csv(CV_result,file = "C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/correlation/CV result_3.csv",row.names = F)
write.csv(CV_result_1,file = "C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/correlation/CV result_1.csv",row.names = F)
write.csv(CV_result_2,file = "C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/correlation/CV result_2.csv",row.names = F)
write.csv(CV_result_4,file = "C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/correlation/CV result_4.csv",row.names = F)

#############################################
r_1 <- c()
r_2 <- c()
r_3 <- c()
for (i in 1:12) {
  index <- which(total$map == paste0("chromosome ",i))
  data <- total[index,]
  model <- lm(gen~poly(phys,3), data = data)
  sm <- summary(model)
  r_3 <- c(r_3,round(sm$adj.r.squared,3))
 
}
