library(dplyr)
library(tidyverse)
####以第一條染色體為例，總共發生70次重組，根據paper將window定義為0.1mb，
#因此平均每個window會發生0.16次重組。接著每個window計算裡面的重組次數，並轉換成機率，如果小於0.05便是hotspot。



###################################################
#poisson hotspot (結果寫入hotspot dataframe中)
######################################################

####### 使用 gbs_breakpoint.csv
data <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_breakpoint.csv")

chr_len <- c(43260640,35954074,36189985,35489479,29733216,30731386,29643843,28434680,22692709,22683701,28357783,27561960)
sum(chr_len)/100000
####

hotspot <- data.frame(matrix(ncol = 4,nrow = 0))
colnames(hotspot) <- c("chr","event_number","start","end")

windowsize = 100000
p_value_threshold = 0.05/3707
dpois(4,)
for (chromosome in 1:12) {
  tmp1 <- data[data$Chr== chromosome,]
  #lamda = length(tmp1$pos)/(chr_len[chromosome]/windowsize)
  #lamda = 6980/(sum(chr_len)/windowsize)
  lamda = nrow(tmp1) / (chr_len[chromosome]/windowsize)
  
  interval <- seq(1,chr_len[chromosome],by = windowsize)
  ###
  for (i in 1:(length(interval)-1)) {
    if( dpois(length(which(between(tmp1$Pos,interval[i],interval[i+1]))),lamda) < p_value_threshold ){
      temp <- data.frame(chr = chromosome,
                         event_number = length(which(between(tmp1$Pos,interval[i],interval[i+1]))),
                                               start = interval[i],end =interval[i+1] )
      hotspot <- rbind(hotspot, temp)
    }else{
      next
    }
  }
}

write.csv(hotspot,file = "gbs_all_co.csv")

###比對
wgs_co <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/all_co_event.csv",header = T)
which(wgs_co$start %in% hotspot$start)
for (i in 1:12) {
  a <- wgs_co[wgs_co$chr == i,]
  b <- hotspot[hotspot$chr == i,]
  print(a[which(a$start %in% b$start),])
}


####################################### 
# local recomb hotspot 
########################################
all_recomb <- read.table("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/local_recomb/gbs_local_recomb.txt",header = T)
all_recomb_0.5 <- read.table("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/gbs_local_recomb_spar_0.5.txt",header = T)
q1 <- quantile(all_recomb$spline, 0.25,na.rm = T)
q3 <- quantile(all_recomb$spline, 0.75,na.rm = T)
q1 <- quantile(all_recomb_0.5$spline, 0.25,na.rm = T)
q3 <- quantile(all_recomb_0.5$spline, 0.75,na.rm = T)
upper <- (q3-q1)*1.5 + q3
three_fold_upper <- mean(all_recomb$spline)*3
three_fold_upper_0.5 <- mean(all_recomb_0.5$spline)*3

all_recomb_0.5[all_recomb_0.5$spline> upper,]
temp_0.1 <- all_recomb[all_recomb$spline> three_fold_upper,]
temp_0.5 <- all_recomb_0.5[all_recomb_0.5$spline> three_fold_upper_0.5,]

#############################################
# overlap hotspot
#############################################
"C:/Users\lingo1st\Dropbox\林冠瑜\gbs_dataset_result\local_recomb\local_hotspot_spar_0.1.csv"
poisson_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_adjusted_poisson_hotspot.csv",header = T)
local_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/local_recomb/local_hotspot_spar_0.1.csv",header = T)


##使用LOCAL中心點，觀察是否落在POISSON區間前後50000bp內
result <- apply(local_hotspot,1,function(x){
  overlapp <- data.frame(matrix(nrow = 0,ncol = 3))
  part_df <- poisson_hotspot %>% filter(chr ==  as.numeric(substr(x[3],4,nchar(x[3]))) )
  for (i in 1:nrow(part_df)) {
    if(between(as.numeric(x[5]),part_df$start[i]-50000,part_df$end[i]+50000) ){
      print(i)
      temp <- intersect(c((as.numeric(x[5])-50000):(as.numeric(x[5])+50000)), c(part_df$start[i]:part_df$end[i]))
      temp_df <- data.frame(chr = as.numeric(substr(x[3],4,nchar(x[3]))),min = min(temp),max = max(temp))
      overlapp <- rbind(overlapp,temp_df)
    }
    
  }
  return(overlapp)
})
###多出來的三個column為空行所造成?
merged_result <- bind_rows(result)
merged_result <- merged_result[,1:3]


################## distance btw two breakpoint ####################
dis <- c()
###data需要排序
sort_data <- data %>% group_by(Chr) %>%
  arrange(Pos,.by_group = T)
for (i in 1:(nrow(sort_data)-1)) {
  dis[i] <- sort_data$Pos[i+1]- sort_data$Pos[i]
}

dis[dis>0 & dis == min(dis[dis>0])]
which.min(dis)
dis[which.max(dis)]
sort_data[3831:3832,]




##########################
#plot
#############################
poisson_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_adjusted_poisson_hotspot.csv",header = T)
breakpoint <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_breakpoint.csv",header = T)
chr_len <- c(43260640,35954074,36189985,35489479,29733216,30731386,29643843,28434680,22692709,22683701,28357783,27561960)
p1 <- ggplot(data = breakpoint,aes(x = Chr,y=Pos/1000000))
for (i in 1:12) {
  p1<- p1+
    geom_point(size = 0.01)
  #geom_segment(x = i , y = 1, xend = i, yend = chr_len[i]/1000000)
}
p1 <- p1+scale_x_discrete(limits = c(1:12))+ylab("Position (Mb)")+xlab("Chromosome")
p1
####try to add violin plot
p1 <- p1+ geom_violin(aes(x = Chr,y = Pos/1000000,group = Chr),fill = "orange",alpha = 0.2)
####add centromere region
cent_min <- c(15.445668,12.625206,17.883934,7.882,11.15,13.201,9.104,11.98,0.993,
              7.62,11.34,11.06)
cent_max <- c(18.05,15.48,20.51,10.06,13.54,17.84,12.71,14.41,3.93,8.69,13.5,
              12.2)

cent <- data.frame(xmin = seq(0.8,11.8,by = 1),xmax = seq(1.2,12.2,by = 1),ymin = cent_min,ymax = cent_max)
p1 + geom_rect(aes(xmin = 0.8, xmax = 1.2, ymin = 15445668/1000000, ymax = 18052668/1000000),alpha = 0.5)
for (i in 1:12) {
  p1<- p1+
    geom_rect(xmin = cent$xmin[i], xmax = cent$xmax[i], ymin = cent$ymin[i], ymax = cent$ymax[i],alpha = 0.5)
}
p1 
##### add hotspot
for (i in 1:nrow(poisson_hotspot)) {
  p1<- p1+
    geom_rect(xmin =poisson_hotspot$chr[i]-0.1 ,xmax =poisson_hotspot$chr[i]+0.1 , ymin = (poisson_hotspot$start[i]/1000000), ymax = (poisson_hotspot$end[i]/1000000),alpha = 0.1,fill = "red",alpha = 0.1)
}
p1

#############################################################################################
poisson_model <- glm(recomb_counts ~ 1, family = poisson)
summary(poisson_model)
negbin_model <- glm.nb(recomb_counts ~ 1)
summary(negbin_model)
zip_model <- zeroinfl(recomb_counts ~ 1 | 1, dist = "poisson")
summary(zip_model)
zinb_model <- zeroinfl(recomb_counts ~ 1 | 1, dist = "negbin")
summary(zinb_model)
aic_values <- c(poisson = AIC(poisson_model),
                negbin = AIC(negbin_model),
                zip = AIC(zip_model),
                zinb = AIC(zinb_model))
print(aic_values)
outliers <- boxplot.stats(hotspot$co_num)$out
recomb_counts <- data[!data %in% outliers]
