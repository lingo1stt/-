library(seqinr)
library(tidyverse)
window_size <- 10000
fas_data_path <- "C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/nuc_test"
sliding = 0
one_side_length <- 50000 
###############data#########################
overlap_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_overlap_hotspot.csv",header = T)
poisson_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file//poisson hotspot/gbs_adjusted_poisson_hotspot.csv",header = T)
local_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/local_recomb_hotspot/gbs_nonoverlap_local_hotspot.csv",header = T)

##合併
hotspot_interval <- rbind(local_hotspot[,c(1,2,3)],poisson_hotspot[,c(2,5,6)],overlap_hotspot[,c(1,2,3)])
hotspot_interval$type <- rep(c("local","poisson","overlap"),times = c(25,152,21))
hotspot_interval$chr <- as.numeric(hotspot_interval$chr)
###加入區間長度
hotspot_interval$length <- apply(hotspot_interval,1,function(x){
  as.numeric(x[3]) - as.numeric(x[2])
})


#################################################
#version 3 :input改成輸入想要的等分數
#################################################
nuc_diversity <- function(cut_number,fas_data_path = NA,single_value = F,file_index ){
  cat("############################################","\n","you are calculating nucleotide diversity...###############################","\n")
  # 取得檔案數量
  folder_path <-fas_data_path
  # 列出路徑下的所有檔案
  files <- list.files(folder_path)
  num_files <- length(files)

  plot_mat <- matrix(ncol = num_files,nrow = cut_number)  
  ###########################################
  for (q in 1:num_files) {
  fas <- read.fasta(paste0(folder_path,"/",files[q]))
  ###先曉得fasta長度
  fas_len <- length(fas[[1]])
  #前置data及apply使用之function
  name <- names(fas)
  sub <- function(x,pos){
    return(x[interval_index[pos]:interval_index[pos+1]-1])
  }
  
  freq <- function(y,ind){
    return(length(which(y == temp_list[[ind]])))
  }
  
  ##計算Pi的迴圈本體
  pi_diversity <- c()
  interval_index <- round(seq(1,fas_len,length = cut_number+1),0)
  ###區間迴圈
  for (i in 1: (length(interval_index)-1) ) {
    ##擷取每個區間的長度(temp_list)
    temp_list <- lapply(fas,FUN =  sub,pos = i)
    uni <- unique(temp_list)
    ##比對
    freq_list_2 <- list()
    for (p in 1:length(uni)) {
      sub_list <- list()
      freq_list_2[[p]] <- sub_list
    }
    ##將每個個體的名字分配給所屬的不重複序列(freq_list_2) (freq_list為不重複序列本身)
    for (l in 1:length(fas)) {
      freq_list <- lapply(uni,FUN = freq,ind = l)
      freq_list_2[[which(freq_list %in% length(temp_list[[1]]))]] <- append(freq_list_2[[which(freq_list %in% length(temp_list[[1]]))]],name[l]) 
      
    }
    ##去掉na之後的長度並轉換成比例，用比例當作freq_list_2的名字
    uni_len <- (sapply(freq_list_2,length))/sum(sapply(freq_list_2,length))
    names(freq_list_2) <- uni_len
    
    ##兩兩比較
    pi_value <- list()
    if(length(uni)==1){
      result = 0
    }else{
      for (j in 1:length(uni)) {
        for (k in (j+1):(length(uni))) {
          if(j<length(uni)){
            diff <-which(Map(`!=`,uni[[j]],uni[[k]])==T)
            pi <- length(diff)/length(uni[[j]])
            listname <- paste0(j,",",k)
            pi_value[[listname]] <- pi
            
          }else{
            break
          }
        }
      }
      ###將pi_value list名稱分開，找到各自數字代表的頻率-->freq_value
      list_name <- names(pi_value)
      freq_index <- as.numeric(unlist(sapply(list_name,function(x){strsplit(x,",")})))
      freq_value <- c()
      for (o in 1:length(freq_index)) {
        freq_value[o]<- as.numeric(names(freq_list_2[freq_index[o]]))
      }
      ##將頻率乘上pi值(xixj*pi)
      result <- mapply(function(x, y, sublist) x * y * sublist, freq_value[seq(1, length(freq_value), 2)], freq_value[seq(2, length(freq_value), 2)], pi_value)
      ###乘上n/n-1
      result <- sum(result)*(length(fas)/(length(fas)-1))
    }
    pi_diversity <- append(pi_diversity,result)
    
    
    
  }
  plot_mat[1:nrow(plot_mat),q] <- pi_diversity
  }
  
  
  assign(paste0("plot_mat"),plot_mat,envir = .GlobalEnv)
  write.csv(plot_mat,file = paste0("diverisity value.csv"))
  return("jon done")
   
}



############################################################################################################



#########local ##############
nuc_diversity(fas_data_path = "C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/fasta file/local_hotspot_fasta" ,
              cut_number = 10
              )

local_diversity_value <- plot_mat
avg_l <- apply(local_diversity_value,1,function(x){
  mean(x)
})
write.csv(local_diversity_value,file = "gbs_local_hotspot_diversity_value.csv")


####poisson###
poisson_hotspot <- poisson_hotspot[poisson_hotspot$co_num >=11,c(-1,-8,-9)]
write.csv(poisson_hotspot,file = "gbs_poisson_hotspot.csv",row.names = F)

nuc_diversity(fas_data_path = "C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/fasta file/poisson_hotspot_fasta" ,
              cut_number = 10
)
poisson_diversity_value <- plot_mat
avg_p <- apply(poisson_diversity_value,1,function(x){
  mean(x)
})
write.csv(poisson_diversity_value,file = "gbs_poisson_hotspot_diversity_value.csv")

#####overlap
nuc_diversity(fas_data_path = "C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/fasta file/overlapp_hotspot_fasta" ,
              single_value = F,
              cut_number = 10
)
overlap_diversity_value = plot_mat
avg_o <- apply(overlap_diversity_value,1,function(x){
  mean(x)
})
write.csv(overlap_diversity_value,file = "gbs_overlap_hotspot_diversity_value.csv")

###random
nuc_diversity(fas_data_path = "C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/fasta file/random_fasta" ,
              single_value = F,
              cut_number = 10
)
random_diversity_value = plot_mat
avg_r <- apply(random_diversity_value,1,function(x){
  mean(x)
})

#####################################################################################################
#plot
####################################################################################################
#trend
plot.df <- data.frame(value = c(avg_l,avg_p,avg_o,avg_r),
                      Type = rep(c("L","P","O","R"),time = c(10,10,10,10)),
                      position = rep(c("-50%","-40%","-30%","-20%","-10%",
                                       "10%","20%","30%","40%","50%"),4))
write.csv(plot.df,file = "div_plot_df.csv",row.names = F)

p1 <- ggplot(data = plot.df,aes(x = position,y = value,fill= Type))
p1 + geom_col(position = "dodge") + xlab("Relative position") + ylab("Nucleotide diversity") + 
  theme(
    text = element_text(size = 8)
  )

##boxplot
p2 <- ggplot(data = plot.df,aes(x = type,y = value))
p2 + geom_boxplot()
####anova
install.packages("multcompView")
library(multcompView)
library(datasets)
model = aov(value~type,data = plot.df)
summary(model)
tukey <- TukeyHSD(model)
###cld用來生成分組的字母
cld <- multcompLetters4(model, tukey)
###tk用來把字母畫上去
Tk <- group_by(plot.df, Type) %>%
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
