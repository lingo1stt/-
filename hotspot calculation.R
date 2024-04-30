####detection of hotspot
library(dplyr)
#####cal function#####
cal <- function(chr){
  data <-  read.csv(paste0("C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation_hotspot_result/",chr,"/imputation/final_sample_plus_parent.",chr,"_imputed_chunks.txt")
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
  #capture.output(breakpoint, file = paste0("C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation result/",chr,"/","breakpoint","_",chr,".csv"))
  ###
  plot_df <- data.frame(chr = rep(as.numeric(substr(chr,4,6)),length(breakpoint[[25]])),pos = unlist(breakpoint[[25]]))
  assign(paste0("plot_",chr), plot_df, envir=.GlobalEnv)
  
}

####以第一條染色體為例，總共發生70次重組，根據paper將window定義為0.1mb，
#因此平均每個window會發生0.16次重組。接著每個window計算裡面的重組次數，並轉換成機率，如果小於0.05便是hotspot。
###result:9個熱點 
cal("chr1")
for (i in 1:12) {
  cal(paste0("chr",i))
}
lamda = length(plot_chr1$pos)/(43260640/100000)
hotspot <- list()
interval <- seq(1,43260640,by=100000)
for (i in 1:(length(interval)-1)) {
  if( dpois(length(which(between(plot_chr1$pos,interval[i],interval[i+1]))),lamda) <0.05 ){
    tmp <- list(c(interval[i],interval[i+1]))
    hotspot <- append(hotspot,tmp)
  }else{
    next
  }
}
##using 1mb as window
#result:只有一個熱點，區間內至少要大於五次重組才能算
lamda = length(plot_chr1$pos)/(43260640/1000000)
hotspot <- list()
interval <- seq(1,43260640,by=1000000)
for (i in 1:(length(interval)-1)) {
  if( dpois(length(which(between(plot_chr1$pos,interval[i],interval[i+1]))),lamda) <0.05 ){
    tmp <- list(c(interval[i],interval[i+1]))
    hotspot <- append(hotspot,tmp)
  }else{
    next
  }
}
dpois(5,1.6)


########################################
#function (須配合cal function 使用)
########################################
chr_len <- c(43260640,35954074,36189985,35489479,29733216,30731386,29643843,28434680,22692709,22683701,28357783,27561960)
###generate plot_df
for (i in 1:12) {
  cal(paste0("chr",i))
}

####
cal_hotspot <- function(windowsize,chromosome,p_value_threshold){
  tmp1 <- get(paste0("plot_chr",chromosome))
  #lamda = length(tmp1$pos)/(chr_len[chromosome]/windowsize)
  lamda = 610/(sum(chr_len)/windowsize)
  hotspot <- list()
  interval <- seq(1,chr_len[chromosome],by = windowsize)
  ###
  for (i in 1:(length(interval)-1)) {
    if( dpois(length(which(between(tmp1$pos,interval[i],interval[i+1]))),lamda) < p_value_threshold ){
      tmp2 <- list(c(length(which(between(tmp1$pos,interval[i],interval[i+1]))),interval[i],interval[i+1]))
      hotspot <- append(hotspot,tmp2)
    }else{
      next
    }
  }
  assign(paste0("hotspot_chr",chromosome),as.data.frame(do.call(rbind,hotspot)),envir = .GlobalEnv)
  tmp3 <- get(paste0("hotspot_chr",chromosome))
  colnames(tmp3) <- c("CO event number","start(bp)","end(bp)")
  #write.csv(tmp3,file =  paste0("C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation result/","chr",chromosome,"/hotspot position_100kbwindow_chr",chromosome,".csv"))
  write.csv(tmp3,file = paste0("C:/Users/lingo1st/OneDrive/桌面/",chromosome,".csv"))
}

 for (i in 1:12) {
  cal_hotspot(100000,i,0.01)
 }

######distance btw two breakpoint
dis <- c()
for (i in 1:(nrow(total_plot_df)-1)) {
  dis[i] <- total_plot_df$pos[i+1]- total_plot_df$pos[i]
}

dis[dis>0 & dis == min(dis[dis>0])]
which.min(dis)
total_plot_df[52:53,]

####hotspot plot
#############plot#############
####run the function
for (i in 1:12) {
  cal(paste0("chr",i))
}
#####
total_plot_df <- rbind(plot_chr1,plot_chr2,plot_chr3,plot_chr4,plot_chr5,plot_chr6,plot_chr7,
                       plot_chr8,plot_chr9,plot_chr10,plot_chr11,plot_chr12)
total_plot_df$chr <- as.character(total_plot_df$chr)
chr_len <- c(43260640,35954074,36189985,35489479,29733216,30731386,29643843,28434680,22692709,22683701,28357783,27561960)
p1 <- ggplot(data = total_plot_df,aes(x = chr,y=pos/1000000))
for (i in 1:12) {
  p1<- p1+
    geom_point(size = 0.5)+ geom_jitter(width = 0.1)
    geom_segment(x = i , y = 1, xend = i, yend = chr_len[i]/1000000)
}
p1 <- p1+scale_x_discrete(limits = c(1:12))+ylab("position (Mb)")+xlab("chromosome")

####add centromere region
cent_min <- c(15.445668,12.625206,17.883934,7.882,11.15,13.201,9.104,11.98,0.993,
              7.62,11.34,11.06)
cent_max <- c(18.05,15.48,20.51,10.06,13.54,17.84,12.71,14.41,3.93,8.69,13.5,
              12.2)

cent <- data.frame(xmin = seq(0.8,11.8,by = 1),xmax = seq(1.2,12.2,by = 1),ymin = cent_min,ymax = cent_max)
p1 + geom_rect(aes(xmin = 0.8, xmax = 1.2, ymin = 15445668/1000000, ymax = 18052668/1000000),alpha = 0.1)
for (i in 1:12) {
  p1<- p1+
    geom_rect(xmin = cent$xmin[i], xmax = cent$xmax[i], ymin = cent$ymin[i], ymax = cent$ymax[i],alpha = 0.1)
}
p1 


##### add hotspot
for (i in 1:12) {
  tmp <- get(paste0("hotspot_chr",i))
  tmp$chr <- rep(i,time = nrow(tmp))
  assign(paste0("hotspot_chr",i),tmp,envir = .GlobalEnv)
}
hotspot_df <- rbind(hotspot_chr1,hotspot_chr2,hotspot_chr3,hotspot_chr4,hotspot_chr5,hotspot_chr6,hotspot_chr7,
                    hotspot_chr8,hotspot_chr9,hotspot_chr10,hotspot_chr11,hotspot_chr12)
colnames(hotspot_df) <- c("CO number","start","end","chr")

for (i in 1:nrow(hotspot_df)) {
  p1<- p1+
    geom_rect(xmin =hotspot_df$chr[i]-0.2 ,xmax =hotspot_df$chr[i]+0.2 , ymin = (hotspot_df$start[i]/1000000), ymax = (hotspot_df$end[i]/1000000),alpha = 0.1,fill = "red",alpha = 0.1)
}
p1B


