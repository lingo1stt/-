library(dplyr)
library(tidyverse)
library(vcfR)
getwd()
###使用abh_genotype.csv 計算 
data <-  read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/gbs_abhgenotype.csv",header = T)
###output格式為matrix，ncol為單條染色體最大之co數量，其他無法填入的位置即為na
breakpoint <- data.frame(Chr = NA,Pos = NA)
col_name <- rownames(data)
#i為染色體, j為snp數量
for (i in 1:12) {
  ##選出該條染色體之子集-->dat
  dat <- data[,names(data)[which(data[1,] == i)] ]
  dat <- dat[-1,]
  rownames(dat) <- data$id[-1]
  ###使用apply，針對每個row(個體)計算co發生位置，並記錄在co_site中
  result <- apply(dat,1,function(x){
    co_site <- c()
    start <- ifelse(x[1] == "-",x[2],x[1])
    for (j in 2:length(x)) {
      if( x[j] == start | x[j] == "-" ){
        next
      }else{
        co_site <- c(co_site,(as.numeric(strsplit(names(x[j]),"_")[[1]][2]) + as.numeric(as.numeric(strsplit(names(x[j-1]),"_")[[1]][2])) )/2  )
        start <- x[j]
        
      }
    }
    return(co_site)
  }
  )
  ###result list長度為個體數，每個sublist中包含所有斷點的位置
  ###將其放在一個Dataframe中
  result <- unlist(result)
  result_df <- data.frame(Chr = rep(i,length(result)),Pos = result)
  breakpoint <- rbind(breakpoint,result_df)

}

write.csv(breakpoint,file = "breakpoint_df.csv")





###檢查用
co_event <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/CO event number.csv",row.names = T)





###################### plot #############################
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




######## CO event vs chromosome length
co <- breakpoint %>% 
  group_by(Chr) %>%
  summarise(total = n())
CO_length_df <- data.frame(len = chr_len,co = co$total[-13])
p2 <- ggplot(data = CO_length_df,aes(x = len/1000000,y=co))
p2 + geom_point() + geom_smooth(method = lm,formula = y~x) + ylab("Amount of CO event")+
  xlab("Chromosome length (Mb)")
cor(CO_length_df$len,CO_length_df$co)


#######CO event vs SNP density (1Mb)
vcf_chr1 <- read.vcfR("C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data/gbs_r.chr1.vcf",verbose = F)
chr1_fix <- as.numeric( vcf_chr1@fix[,"POS"])
##每條染色體分區間計算裡面的co以及snp數
snp_num <- c()
co_num <- c()
chr <- c()
for (i in 1:12) {
  temp_vcf <- read.vcfR(paste0("C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data/gbs_r.chr",i,".vcf"),verbose = F)
  temp_fix <-  as.numeric( temp_vcf@fix[,"POS"])
  temp_length <- chr_len[i]
  ###逐個區間計算
  for (j in seq(1,temp_length,by = 1000000)) {
    ##interval內snp數量
    snp_num <- c(snp_num,length(which(temp_fix>j & temp_fix<j+1000000)) )
    ##interval 內co數量
    temp_df <- breakpoint %>% filter(Chr == i)
    co_num <- c(co_num, length(temp_df$Pos[temp_df$Pos>j & temp_df$Pos < j+1000000]) )
  }
  chr <- c(chr,rep(i,length(seq(1,temp_length,by = 1000000))) )
  
}
###合併資料 + 繪圖
co_density_df <- data.frame(snp = snp_num,co = co_num,chr = chr)
p3 <- ggplot(data = co_density_df, aes(x = snp, y = co))
p3 + geom_point() + geom_smooth(method = lm, formula = y~x)+
  xlab("SNP number within 1Mb interval") + ylab("CO number within 1Mb interval")
cor(co_density_df$snp,co_density_df$co) 

### 儲存資料
write.csv(co_density_df, file = "co_vs_density.csv")
