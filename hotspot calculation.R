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
poisson_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/poisson hotspot/gbs_poisson_hotspot.csv",header = T)
breakpoint <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_breakpoint.csv",header = T)
chr_len <- c(43260640,35954074,36189985,35489479,29733216,30731386,29643843,28434680,22692709,22683701,28357783,27561960)
breakpoint_2 <- breakpoint
breakpoint_2$Chr <-breakpoint_2$Chr *2 -1
p1 <- ggplot(data = breakpoint_2,aes(x = Chr,y=Pos/1000000))
for (i in 1:12) {
  p1<- p1+
    geom_point(size = 0.01,position = position_dodge(width = 2))
  #geom_segment(x = i , y = 1, xend = i, yend = chr_len[i]/1000000)
}
p1 <- p1 +  scale_x_continuous(
  breaks = seq(1, 23, by = 2),  # 原始 x 軸位置
  labels = seq(1, 12, by = 1)   # 新的 x 軸標籤
)+ylab("Position (Mb)")+xlab("Chromosome")
####try to add violin plot
p1 <- p1+ geom_violin(aes(x = Chr,y = Pos/1000000,group = Chr),fill = "orange",alpha = 0.2)
####add centromere region
cent_min <- c(15.445668,12.625206,17.883934,7.882,11.15,13.201,9.104,11.98,0.993,
              7.62,11.34,11.06)
cent_max <- c(18.05,15.48,20.51,10.06,13.54,17.84,12.71,14.41,3.93,8.69,13.5,
              12.2)

cent <- data.frame(xmin = seq(1,23,by = 2)-0.3,xmax = seq(1,23,by = 2)+0.3,ymin = cent_min,ymax = cent_max)
for (i in 1:12) {
  p1<- p1+
    geom_rect(xmin = cent$xmin[i], xmax = cent$xmax[i], ymin = cent$ymin[i], ymax = cent$ymax[i],alpha = 0.5)
}
p1 
##### add poisson hotspot
poisson_hotspot$chr <- poisson_hotspot$chr*2-1
poisson_hotspot$pos <- (poisson_hotspot$start+poisson_hotspot$end) / 2
for (i in 1:nrow(poisson_hotspot)) {
  p1<- p1+
    geom_rect(xmin =poisson_hotspot$chr[i]-0.25 ,xmax =poisson_hotspot$chr[i]+0.25 , ymin = (poisson_hotspot$start[i]/1000000), ymax = (poisson_hotspot$end[i]/1000000),alpha = 0.1,fill = "red",alpha = 0.1)
}

p1

#########add local recomb rate
install.packages("ggridges")
library(ggridges)

local_recomb <- read.table("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/local_recomb_hotspot/gbs_local_recomb.txt",sep = " ",header = T)
local_recomb$pos <- as.numeric(substr(local_recomb$map,4,nchar(local_recomb$map)) ) *2 -2
p1 + geom_vridgeline(data = local_recomb %>%
                       filter(pos==0),
                     aes(x = pos,y = phys/1000000,width = spline/60),fill="lightblue",linewidth = 0.5)
for (i in unique(local_recomb$pos)) {
  p1 <- p1 + geom_vridgeline(data = local_recomb %>%
                               filter(pos==i),
                             aes(x = pos,y = phys/1000000,width = spline/60),linewidth = 0.3,fill = "lightblue")
}
p1
###add local hotspot
local_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/local_recomb_hotspot/gbs_nonoverlap_local_hotspot.csv",header = T)
local_hotspot$chr <- local_hotspot$chr*2-1
local_hotspot$pos <- (local_hotspot$start+local_hotspot$end)/2
######geom_rect形式
for (i in 1:nrow(local_hotspot)) {
  p1<- p1+
    geom_rect(xmin =local_hotspot$chr[i]-0.25 ,xmax =local_hotspot$chr[i]+0.25 , ymin = (local_hotspot$start[i]/1000000), ymax = (local_hotspot$end[i]/1000000),alpha = 0.1,fill = "red",alpha = 0.1,size = 0.8)
}

#########overlap hotspot
overlap <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_overlap_hotspot.csv",header = T)
overlap$chr <- overlap$chr *2-1
for (i in 1:nrow(overlap)) {
  p1<- p1+
    geom_rect(xmin =overlap$chr[i]-0.25 ,xmax =overlap$chr[i]+0.25 , ymin = (overlap$start[i]/1000000), ymax = (overlap$end[i]/1000000),alpha = 0.1,fill = "red",alpha = 0.1,size = 0.8)
}
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
