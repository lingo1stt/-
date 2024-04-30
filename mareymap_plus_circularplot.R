install.packages("MareyMap")

install.packages("tkrplot")
install.packages("tools")
install.packages("ggpubr")
library(MareyMap)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(circlize)
startMareyMapGUI()

####建立物理距離跟遺傳距離關係圖
name <- c("x3","ax+b")
r_square <- data.frame(matrix(nrow = 0,ncol = length(name))) 
colnames(r_square) <- name
for (i in 1:12) {
  data <- read.csv(paste0('C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation_hotspot_result/','chr',i,'/local/genetic_distance_chr',i,'.csv')) %>% na.omit()
  data <- data[seq(1,nrow(data),by = 2),]
  assign(paste0("data_for_mareymap_",i),data,envir = .GlobalEnv)
  #write.table(data,paste0('C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation result/','chr',i,'/data_for_mareymap_',i,'.txt'),sep = " ", quote = TRUE, row.names = FALSE)
  m1 <- summary(lm(gen~poly(phys,3),data = data))$adj.r.squared
  m2 <- summary(lm(gen~phys,data = data))$adj.r.squared
  r_square[i,1:2] <- c(m1,m2)
  assign(paste0("p",i),ggplot(data = data,aes(x = phys/1000000, y= gen))+ geom_point()+geom_smooth(method = lm,formula = y~poly(x,3))+xlab("physical position(mb)")+ylab("genetic position(cM)")+annotate("text",x = (max(data$phys)/1000000)-5,y = min(data$gen),label = paste0("R square = ",round(m1,3))))
}
###plot
p1
plot <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,labels = LETTERS[1:12],ncol = 6,nrow = 2)
pdf(file = 'C:/Users/lingo1st/Dropbox/林冠瑜/noisymputer result_2/imputation result/correlation_2.pdf',width = 30)
plot
dev.off()

write.csv(r_square,paste0('C:/Users/lingo1st/Dropbox/林冠瑜/noisymputer result_2/imputation result/correlation comparision.csv'))

###
m1 <- lm(gen~poly(phys,3),data = data_for_mareymap_2)
P3 <- ggplot(data = data_for_mareymap_2,aes(x = phys/1000000, y= gen)) + geom_point()+geom_smooth(method = lm,formula = y~poly(x,3))
P3
t <- summary(m1)
t
m2 <- lm(gen~phys,data = data_for_mareymap_2)
summary(m2)

#### test for local recomb hotspothttp://127.0.0.1:20937/graphics/plot_zoom_png?width=1200&height=900
test <- read.table("C:/Users/lingo1st/Dropbox/林冠瑜/noisymputer result_2/imputation result/chr1/local_recomb_chr1.txt",header = T)
boxplot(test$spline)
boxplot(test$loess)

q1 <- quantile(test$loess, 0.25,na.rm = T)
q3 <- quantile(test$loess, 0.75,na.rm = T)
q1 <- quantile(test$spline, 0.25,na.rm = T)
q3 <- quantile(test$spline, 0.75,na.rm = T)
upper <- (q3-q1)*1.5 + q3
test[test$spline> upper]
test[test$loess> upper,]

### read the local_recombination file 
install.packages("circlize")
library("circlize")
data <- read.table("C:/Users/lingo1st/Dropbox/林冠瑜/noisymputer result_2/imputation result/all_local_recomb.txt",header = T)
########建立circular plot 外圈
#建立新資料集
data$map <- gsub("chromosome ","chr",data$map)
cytoband_df <- data.frame(chrom = data$map,chromStart = data$phys,chromoEnd = data$phys,value1 = data$spline)

#畫出重組發生位點之barplot
cal <- function(chr){
  data <-  read.csv(paste0("C:/Users/lingo1st/Dropbox/林冠瑜/noisymputer result_2/imputation result/",chr,"/final_sample_plus_parent.",chr,"_imputed_chunks.txt")
                    ,sep = " ")
  ###data 前處理
  p <- strsplit(data$genotypes,"")
  for (i in 1:length(p)) {
    for (j in 1:length(p[[1]])) {
      data[i,j+2] <- p[[i]][j]
    }
  }
  row_name <- data$X.
  rownames(data) <- row_name
  data <- data[,-1:-2]  %>% as.data.frame()
  
  ###
  breakpoint <- list()
  start <- c()
  col_name <- rownames(data)
  for (i in 1:ncol(data)) {
    start <- data[1,i]
    breakpoint[[i]] <- NaN*seq(nrow(data))
    for (j in 1:nrow(data)) {
      if( data[j,i] == start ){
        next
      }else{
        breakpoint[[i]][j] <- (as.numeric(substr(col_name[j],5, nchar(col_name[j]))) + as.numeric(substr(col_name[j-1],5, nchar(col_name[j-1]))))/2
        start <- data[j,i]
      }
    }
  }
  ####
  for (i in 1:24) {
    breakpoint[[i]] <- breakpoint[[i]][which(!is.na(breakpoint[[i]]))]
  }
  ###
  p <- unlist(breakpoint)
  breakpoint[[" total without 3&24"]] <- p[order(p)]
  capture.output(breakpoint, file = paste0("C:/Users/lingo1st/Dropbox/林冠瑜/noisymputer result_2/imputation result/",chr,"/","breakpoint","_",chr,".csv"))
  ###
  plot_df <- data.frame(chr = rep(as.numeric(substr(chr,4,6)),length(breakpoint[[25]])),pos = unlist(breakpoint[[25]]))
  assign(paste0("plot_",chr), plot_df, envir=.GlobalEnv)
  
}
df_name <- paste0("plot_chr",1:12)
plot_chr_list <- mget(df_name) 
data_for_coevent <- do.call(rbind,plot_chr_list)
data_for_coevent$chr <- paste0("chr",data_for_coevent$chr)
data_for_coevent$chromEnd <- data_for_coevent$pos
colnames(data_for_coevent) <- c("chrom","chromStart","ch")


###找出hotspot周圍兩marker的位置作為區間
hotspot_interval <- list()
for (i in 1:nrow(hotspot_region)) {
  hotspot_interval[[i]] <- c(hotspot_region$map[i],cytoband_df$chromStart[ (as.numeric(rownames(hotspot_region)[i])-1) ],cytoband_df$chromStart[ as.numeric(rownames(hotspot_region)[i])+1 ])
}

df_hotspot_interval <- as.data.frame(matrix(unlist(hotspot_interval),ncol = 3,byrow = T))
colnames(df_hotspot_interval) <- c("chr","start","end")
df_hotspot_interval$start <- as.numeric(df_hotspot_interval$start)
df_hotspot_interval$end <- as.numeric(df_hotspot_interval$end)
###將起始/結束位點column結合到hotspot_region dataframe中
hotspot_region$start <- df_hotspot_interval$start
hotspot_region$end <- df_hotspot_interval$end
write.csv(hotspot_region,"C:/Users/lingo1st/Dropbox/林冠瑜/noisymputer result_2/imputation result/hotspot_interval.csv")

### intersection btw local hotspot & CO event interval(因為區間都是100kb因此只顯示起始點)
CO_start <- c(1200001,26400001,28900001,32100001,40200001,8100001,28100001,28500001,32200001,35900001,19800001,29700001,
              20500001,27400001,27800001,1500001,2100001,2700001,24700001,27900001,700001,24600001,11900001,16700001,
              21200001,25100001,2800001,19600001)
CO_end <- CO_start+100000
chr <- c(rep("chr1",5),rep("chr2",2),rep("chr3",3),rep("chr4",2),rep("chr5",2),rep("chr6",1),rep("chr7",5),rep("chr8",2),rep("chr10",2),
         rep("chr11",2),rep("chr12",2))
data_CO <- data.frame(chr = chr,start = CO_start,end = CO_end)
chr5_CO_interval <- c(20500001,27400001)
chr7_CO_interval <- c(1500001,2100001,2700001,24700001,27900001)
chr12_CO_interval <- c(2800001,19600001)
##首先，起始點不能在第一區間的結束與第二區間的起始點之間

for (i in 1:nrow(df_hotspot_interval)) {
  a <- get(paste0(df_hotspot_interval$chr[i],"_CO_interval"))
  target <- a[which.min(abs(a-df_hotspot_interval$start[i]))]
  if(df_hotspot_interval$start[i] < target){
    df_hotspot_interval$start[i] <- target
  }else if(df_hotspot_interval$start[i]  > (target+100000)){
    df_hotspot_interval$start[i] <- NA
    next
    }
  
  if(df_hotspot_interval$end[i] > target+100000){
    df_hotspot_interval$end[i] <- target+100000
  }
  else{
    next
  }
}

###前面會將不符合規定的區間標記上Na，因此將具有NA的欄位去除
df_hotspot_interval <- na.omit(df_hotspot_interval)
write.table(df_hotspot_interval,"C:/Users/lingo1st/Dropbox/林冠瑜/noisymputer result_2/imputation result/hotspot_interval.txt")
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
