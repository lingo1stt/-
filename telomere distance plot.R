################
#telomere距離與重組率關係
################
install.packages("tidyverse")
library(tidyverse)
library(dplyr)
##匯入local recombination rate資料
local_recomb <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/local_recomb_hotspot/gbs_local_recomb.txt",header = T,sep = " ")
###先找出
chr_len <- c(43260640,35954074,36189985,35489479,29733216,30731386,29643843,28434680,22692709,22683701,28357783,27561960)
####add centromere region
cent_min <- c(15.445668,12.625206,17.883934,7.882,11.15,13.201,9.104,11.98,0.993,
              7.62,11.34,11.06) * 1000000
cent_max <- c(18.05,15.48,20.51,10.06,13.54,17.84,12.71,14.41,3.93,8.69,13.5,
              12.2)* 1000000
##################染色體尾端作圖 (尾端到centromere進行分群)
#先定義0-50%的區間

category <- matrix(nrow = 12,ncol = 11)
for (i in 1:12) {
  per <- seq(cent_max[i],chr_len[i],length = 11)
  category[i,1:11] <- per
}

category <- as.data.frame(category)


###使用cut函數將local recomb rate分類 (依照該marker相對於染色體尾端的距離)
final_df <- matrix(NA,ncol = 8,nrow = 0) %>% as.data.frame()
colnames(final_df) <- c(colnames(local_recomb),"telo_pos")
for (i in 1:12) {
  temp <- local_recomb %>% filter(map == paste0("chr",i))
  temp$telo_pos <- cut(temp$phys,breaks =  category[i,],labels =c("45-50","40-45","35-40","30-35","25-30","20-25","15-20","10-15","5-10","0-5") )
  final_df <- rbind(final_df,temp)
  
}
final_df[which(final_df$spline < 0),"spline"] <- 0
###plot
p1 <- ggplot(data = final_df,aes(x = telo_pos,y = spline))
p1 + geom_boxplot(outliers = F)+ stat_summary(fun.y=mean, geom="point", shape=16, size=4,col = "red") +  
  scale_x_discrete(limits=c("0-5","5-10","10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50")) +
  xlab("Distance to closest telomere (% chromosome length)") + ylab("Local recombination rate (cM)")

final_df %>% group_by(telo_pos) %>% summarise(mean = mean(spline))
model <- aov(data = final_df, formula = spline~telo_pos)
summary(model)
TukeyHSD(model)

#####################染色體前端作圖
#先定義0-50%的區間


category <- matrix(nrow = 12,ncol = 11)
for (i in 1:12) {
  per <- seq(0,cent_min[i],length = 11)
  category[i,1:11] <- per
}

category <- as.data.frame(category)


###使用cut函數將local recomb rate分類 (依照該marker相對於染色體尾端的距離)
final_df <- matrix(NA,ncol = 8,nrow = 0) %>% as.data.frame()
colnames(final_df) <- c(colnames(local_recomb),"telo_pos")
for (i in 1:12) {
  temp <- local_recomb %>% filter(map == paste0("chr",i))
  temp$telo_pos <- cut(temp$phys,breaks =  category[i,],labels =rev(c("45-50","40-45","35-40","30-35","25-30","20-25","15-20","10-15","5-10","0-5")) )
  final_df <- rbind(final_df,temp)
  
}
final_df[which(final_df$spline < 0),"spline"] <- 0
####boxplot
p2 <- ggplot(data = final_df,aes(x = telo_pos,y = spline))
p2+ geom_boxplot(outliers = F)+ stat_summary(fun.y=mean, geom="point", shape=16, size=4,col = "red") +  
  scale_x_discrete(limits=c("0-5","5-10","10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50")) +
  xlab("Distance to closest telomere (% chromosome length)") + ylab("Local recombination rate (cM)")


###統計平均及檢定
final_df %>% group_by(telo_pos) %>% summarise(mean = mean(spline))
model <- aov(data = final_df, formula = spline~telo_pos)
summary(model)
f <- TukeyHSD(model)
f <- as.data.frame(f)
f[[1]][which(f[[1]][,4] < 0.05),4]

####每項之marker數量
table(final_df$telo_pos)

#############################################################################################
#綜合
per <- seq(0,1,by = 0.05)

category <- matrix(nrow = 12,ncol = 21)
for (i in 1:12) {
  category[i,1:21] <- chr_len[i]*per
}
###使用cut函數將local recomb rate分類 (依照該marker相對於染色體尾端的距離)
final_df <- matrix(NA,ncol = 8,nrow = 0) %>% as.data.frame()
colnames(final_df) <- c(colnames(local_recomb),"telo_pos")
labels =c("45-50","40-45","35-40","30-35","25-30","20-25","15-20","10-15","5-10","0-5")
for (i in 1:12) {
  temp <- local_recomb %>% filter(map == paste0("chr",i))
  temp$telo_pos <- cut(temp$phys,breaks =  category[i,],labels =c(rev(labels),'50-55',"55-60","60-65","65-70","70-75",'75-80','80-85','85-90','90-95','95-100')) 
  final_df <- rbind(final_df,temp)
  
}
final_df <- final_df[-1,]
p1 <- ggplot(data = final_df,aes(x = telo_pos,y = spline))
p1 + geom_boxplot(outliers = F)+ stat_summary(fun.y=mean, geom="point", shape=16, size=4,col = "red") +  
  scale_x_discrete(limits=c("0-5","5-10","10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50",
                            '50-55',"55-60","60-65","65-70","70-75",'75-80','80-85','85-90','90-95','95-100')) +
  xlab("Distance to closest telomere (% chromosome length)") + ylab("Local recombination rate (cM)")

final_df %>% group_by(telo_pos) %>% summarise(mean = mean(spline))
model <- aov(data = final_df, formula = spline~telo_pos)
summary(model)
TukeyHSD(model)

p <- final_df %>% filter(telo_pos == "50-55")
