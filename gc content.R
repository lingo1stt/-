################
#gc content
################
##找出熱點區間
library(seqinr)
###snp matrix
abh_geno <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_abhgenotype.csv",header = T)

############叫出bkpt
breakpoint <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_breakpoint.csv",header = T)

data <- read.fasta("C:/Users/lingo1st/Dropbox/碩論/np ir64/REF_genome.fa/OsNB1.0.fa")

#########建立熱點區間 (加入overlapped interval)
overlap_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_overlap_hotspot.csv",header = T)
poisson_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file//poisson hotspot/gbs_poisson_hotspot.csv",header = T)
local_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/local_recomb_hotspot/gbs_nonoverlap_local_hotspot.csv",header = T)

#names(local_hotspot)[names(local_hotspot) == "map"] <- "chr"
#local_hotspot$chr <- gsub("chr","",local_hotspot$chr)
###綜合成hotspot_interval
hotspot_interval <- rbind(local_hotspot[,c(1,2,3)],poisson_hotspot[,c(1,4,5)],overlap_hotspot[,c(1,2,3)])
hotspot_interval$type <- rep(c("local","poisson","overlap"),times = c(25,93,21))
hotspot_interval$chr <- as.numeric(hotspot_interval$chr)
####
cent_min <- c(15.445668,12.625206,17.883934,7.882,11.15,13.201,9.104,11.98,0.993,
              7.62,11.34,11.06)
cent_min <- cent_min*1000000
cent_max <- c(18.05,15.48,20.51,10.06,13.54,17.84,12.71,14.41,3.93,8.69,13.5,
              12.2)
cent_max <- cent_max*1000000
######################################################################################################################
#計算gc
#################local #####################################

folder_path <-"C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/fasta file/local_hotspot_fasta"

files <- list.files(folder_path)
num_files <- length(files)
local_gc <- c()
for (i in 1:num_files) {
  print(i)
  temp_fasta <- read.fasta(paste0(folder_path,"/",files[i]))
  test <- unlist(lapply(temp_fasta,FUN = function(x){
    GC(x)
  }))
  result <- mean(test)
  local_gc <- c(local_gc,result)
  
}

mean(local_gc)

##################poisson #############################
folder_path <-"C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/fasta file/poisson_hotspot_fasta"

files <- list.files(folder_path)
num_files <- length(files)
poisson_gc <- c()
for (i in 1:num_files) {
  print(i)
  temp_fasta <- read.fasta(paste0(folder_path,"/",files[i]))
  test <- unlist(lapply(temp_fasta,FUN = function(x){
    GC(x)
  }))
  result <- mean(test)
  poisson_gc <- c(poisson_gc,result)
  
}
 mean(poisson_gc)
 
 
 ####################overlap ###########################
 folder_path <-"C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/fasta file/overlapp_hotspot_fasta"
 
 files <- list.files(folder_path)
 num_files <- length(files)
 overlap_gc <- c()
 for (i in 1:num_files) {
   print(i)
   temp_fasta <- read.fasta(paste0(folder_path,"/",files[i]))
   test <- unlist(lapply(temp_fasta,FUN = function(x){
     GC(x)
   }))
   result <- mean(test)
   overlap_gc <- c(overlap_gc,result)
   
 }
 mean(overlap_gc)
 
##########whole gc ################
whole_gc <- c()
 for (i in 1:12) {
   print(i)
   whole_gc <- c(whole_gc, GC(data[[i]][-c(cent_min:cent_max)]))
   
 }
 
######plot
library(tidyverse)
 library(multcompView)
 library(datasets)
 install.packages("ggpubr")
 library(ggpubr)
 model = aov(value~type,data = plot_df)
 summary(model)
 tukey <- TukeyHSD(model)
 ###cld用來生成分組的字母
 cld <- multcompLetters4(model, tukey)
 ###tk用來把字母畫上去
 Tk <- group_by(plot_df, type) %>%
   summarise(mean=mean(value), quant = quantile(value, probs = 0.75)) %>%
   arrange(desc(mean))
 cld <- as.data.frame.list(cld$type)
 Tk$cld <- cld$Letters
 p2+geom_boxplot() + scale_x_discrete(limits = c("O","L","P","R")) + 
   xlab("Type") + ylab("Nucleotide diversity")+
   geom_point(size = 1) +  geom_text(data = Tk , 
                                     aes(x =Type, y = quant, label = cld),
                                     color = "red",vjust = -5)+ 
   theme(
     text = element_text(size = 8)
   )
 
plot_df <- data.frame(value = c(local_gc,poisson_gc,overlap_gc,whole_gc)*100,
                       Type = rep(c("L","P","O","Non-hotspot"),times = c(25,76,20,12)))

p1 <- ggplot(data = plot_df,aes(x = Type,y =value,fill = Type))
p1 + geom_boxplot() +geom_point(size = 1)+
   scale_x_discrete(limits = c("O","L","P","Non-hotspot"))+ stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("O","Non-hotspot")),method = "wilcox.test")+
  xlab("Type") + ylab("GC content")
  




