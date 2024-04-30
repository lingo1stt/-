library(dplyr)
library(tidyverse)

####### #################
#function (輸入chr+染色體數字可得到每個個體在該染色體的斷點位置)
#########################
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
  #write.csv(data,file = paste0("C:/Users/lingo1st/Dropbox/林冠瑜/noisymputer result_2/imputation result/",chr,"/",chr,"_snp matrix.csv"))
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
  for (i in 1:26) {
    breakpoint[[i]] <- breakpoint[[i]][which(!is.na(breakpoint[[i]]))]
  }
  ###
  p <- unlist(breakpoint[c(-3,-24)])
  breakpoint[[" total without 3&24"]] <- p[order(p)]
  #capture.output(breakpoint, file = paste0("C:/Users/lingo1st/Dropbox/林冠瑜/noisymputer result/",chr,"/","breakpoint","_",chr,".csv"))
  ###
  plot_df <- data.frame(chr = rep(as.numeric(substr(chr,4,6)),length(breakpoint[[27]])),pos = unlist(breakpoint[[27]]))
  assign(paste0("plot_",chr), plot_df, envir=.GlobalEnv)
  
}
cal("chr1")

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
  geom_point(size = 1)+
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
p1

########## correlation btw snp density & breakpoint with 5mb interval ###############
library(vcfR)

###test using chr1
vcf1 <- read.vcfR(file.choose(),verbose = FALSE)
pos <- getFIX(vcf1) %>% data.frame()
pos$POS <- as.numeric(pos$POS)
snp_density <- table(cut(pos$POS,seq(0,45000000,5000000)))
break_count <- table(cut(plot_chr1$pos,seq(0,45000000,5000000)))
cor_df <- data.frame(breakpoint = break_count[1:8],snp_density = snp_density[1:8])
cor_df_1 <- cor_df[,c(2,4)]

###to generate the snp number in 5mb window of each chromosome
for (i in 1:12) {
  vcf <- read.vcfR(paste0("C:/Users/lingo1st/Desktop/NOISYmputer/NOISYmputer_data/final_sample_plus_parent.chr",i,".vcf"),verbose = FALSE)
  pos <- getFIX(vcf) %>% data.frame()
  pos$POS <- as.numeric(pos$POS)
  snp_density <- table(cut(pos$POS,seq(0,45000000,5000000)))
  assign(paste0("cor_df",i),data.frame(snp_density = snp_density[1:9]))
  tmp <- get(paste0("cor_df",i))
  assign(paste0("cor_df",i),tmp[,2], envir=.GlobalEnv)
}
###generate the breakpoint number in each window
for (i in 1:12) {
  tmp <- get(paste0("plot_chr",i))
  break_count <- table(cut(tmp$pos,seq(0,45000000,5000000)))
  assign(paste0("breakpt",i),data.frame(breakpt = break_count[1:9]))
  tmp2 <- get(paste0("breakpt",i))
  assign(paste0("breakpt",i),tmp2[,2],envir=.GlobalEnv)
}

c(cor_df1,cor_df2,cor_df3,cor_df4)
append(cor_df1,cor_df2)

###combine two variable into one dataframe
total_cor <- breakpt1
total_snp <- cor_df1
for (i in 2:12) {
  total_cor <- append(total_cor,get(paste0("breakpt",i)))
  total_snp <- append(total_snp,get(paste0("cor_df",i)))
}
total_cor_df <- data.frame(snp = total_snp,breakpt = total_cor,chr = as.factor(rep(1:12,each = 9)))
###linear model
model2 <- lm(data = total_cor_df,breakpt~snp)
summary(model2)
p2 <- ggplot(data = total_cor_df,aes(x=snp,y=breakpt),group = chr)
p2+geom_point(aes(col = chr)) + geom_smooth(method = lm,formula = y~x) +
ggtitle("Correlation between SNP number and breakpoint number within 5Mb window")+
  xlab("SNP number")+ ylab("Breakpoint number")+
  scale_color_brewer(palette="Paired")

