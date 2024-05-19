#mareymap必須package
install.packages("MareyMap")
install.packages("tcltk")
install.packages("tkrplot")

install.packages("ggpubr")
library(tcltk)
library(MareyMap)
library(tkrplot)
library(tidyverse)
library(ggpubr)
library(circlize)
library(dplyr)
startMareyMapGUI()

####針對marey map格式進行調整
#column分別為 "set" "map" "mkr" "phys" "gen" "vld"
genetic_map <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_whole_genetic_distance.csv")

###全部可以弄再一起，先決定每條染色體的marker數量

set <- data.frame(set = rep("Oryza Sativa",nrow(genetic_map)) )
vld <- data.frame(vld = rep("TRUE",nrow(genetic_map)))
genetic_map <- cbind(set,genetic_map[,-1],vld)
genetic_map <- genetic_map %>% select(everything()[c(1,2,4,5,3)])
genetic_map <- cbind(genetic_map,vld)
colnames(genetic_map) <- c("set" ,"map", "mkr" ,"phys", "gen", "vld")

write.table(genetic_map, file = "data_for_mareymap.txt",quote = T,row.names = F,sep = " ")

 





########################### plot ####################### 
install.packages("circlize")
library("circlize")
data <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/local_recomb/local_hotspot_spar_0.1.csv",header = T)
########建立circular plot 外圈
#建立新資料集
data$map <- gsub("chromosome ","chr",data$map)
cytoband_df <- data.frame(chrom = data$map,chromStart = data$phys,chromoEnd = data$phys,value1 = data$spline)



#########################作圖主體###########################################
circos.par("start.degree" = 90,gap.degree = c(1,1,1,1,1,1,1,1,1,1,1,16))
circos.initializeWithIdeogram(cytoband_df,plotType = c("axis", "labels"),axis.labels.cex = 0.3)
circos.track(sectors =cytoband_df$chrom , x = cytoband_df$chromStart, y = cytoband_df$value1,
             panel.fun = function(x, y) {
               circos.lines(x, y,area = T)
               circos.rect(hotspot_region$start[i], -5 ,hotspot_region$end[i], 40,col = "red",lwd = 1,sector.index = hotspot_region$map[i],border = NA)
             })

for (i in 1:nrow(df_hotspot_interval)) {
  circos.rect(hotspot_region$start[i], -5 ,hotspot_region$end[i], 40,col = "red",lwd = 1,sector.index = hotspot_region$map[i],border = NA)
}
for (i in 1:nrow(data_CO)) {
  circos.rect(data_CO$start[i], -5 ,data_CO$end[i], 40,col = "blue",lwd = 2,sector.index = data_CO$chr[i],border = NA,lty = 1)
}

circos.yaxis(side = "left", sector.index = "chr1",labels.cex = 0.3)
circos.genomicDensity(data = cytoband_df,col = "brown3")
circos.yaxis(side = "left", sector.index = "chr1",labels.cex = 0.15,track.index = 3,labels.niceFacing = T)
circos.clear()

i=1
###

q1 <- quantile(data$spline, 0.25,na.rm = T)
q3 <- quantile(data$spline, 0.75,na.rm = T)
upper <- (q3-q1)*1.5 + q3
hotspot_region <- data[data$spline> upper,c(1:5,8)]
write.csv(hotspot_region,"C:/Users/lingo1st/Dropbox/林冠瑜/noisymputer result_2/imputation result/local_hs_region.csv")

a <- read.table(file.choose())
