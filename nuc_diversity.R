######################################################
# nucleotide diversity
######################################################
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
##################### function ################################
###VERSION2: 改成計算整個fasta檔的single nucleotide value，或是分成window size 計算
nuc_diversity <- function(window_size=NA,sliding = 0,fas_data_path = NA, one_side_length,single_value = F,file_index ){
  cat("############################################","\n","you are calculating nucleotide diversity...###############################","\n")
  # 取得檔案數量
  folder_path <-fas_data_path
  # 列出路徑下的所有檔案
  files <- list.files(folder_path)
  #將files中的順率調整成跟hotspotinterval 中的一樣(用start位置來分)
  file_name_start_pos <- unlist(lapply(strsplit(files,split = "_"),function(x){
    x[1]
  }))
  ##file env為hotspot interval dataframe中的名字
  file_env <- as.character(hotspot_interval$start[which(hotspot_interval$type == "local")])
  order <- c()
  for (cor in 1:length(file_env)) {
    order <- c(order,which(file_env[cor] == file_name_start_pos))
  }
  
  files<- files[order]
  num_files <- length(files)
  if(sliding == 0){
    plot_mat <- matrix(ncol = num_files,nrow = length(seq(1,(2*one_side_length)-sliding,by = window_size))  )
  }else{
    plot_mat <- matrix(ncol = num_files,nrow = length(seq(1,(2*one_side_length)-sliding,by = sliding))  )
  }
  ###########################################
    fas <- read.fasta(paste0(folder_path,"/",files[file_index]))
    #前置data及apply使用之function
    name <- names(fas)
    sub <- function(x,pos){
      return(x[pos:(pos+(window_size-1))])
    }
    
    freq <- function(y,ind){
      return(length(which(y == temp_list[[ind]])))
    }
    
    ##計算Pi的迴圈本體
    pi_diversity <- c()
    for (i in seq(1,(2*one_side_length)-window_size+1,length = nrow(plot_mat))) {
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
        freq_list_2[[which(freq_list %in% window_size)]] <- append(freq_list_2[[which(freq_list %in% window_size)]],name[l]) 
        
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
    plot_mat[1:nrow(plot_mat),file_index] <- pi_diversity
    
  
  
  assign(paste0("plot_mat"),plot_mat,envir = .GlobalEnv)
  write.csv(plot_mat,file = paste0("diverisity value.csv"))
  return("jon done")
}



################poisson hotspot######################
nuc_diversity(window_size = 10000,one_side_length=50000,fas_data_path = "C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/poisson_hotspot_fasta" ,sliding = 0)
poisson_plot_mat <- as.data.frame(plot_mat) 
poisson_recomb_value <- apply(poisson_plot_mat,1,function(x){
  mean(x)
})
####################local recombination################################
nuc_diversity(window_size = round(hotspot_interval$length[f]/10,0),
               one_side_length=round(hotspot_interval$length[f]/2,0),
               fas_data_path = "C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/fasta file/local_hotspot_fasta" ,
               sliding = 0,
              file_index = 1)
local_recomb_value
locate <- c("-50000","-40000","-30000","-20000","-10000",
            "10000","20000","30000","40000","50000")
random = c(0.000009,0.000011,0.000008,0.0000075,0.000013,0.000008,0.000014,0.0000122,0.000013,0.000012)
plot_df <- data.frame(value = c(local_recomb_value,poisson_recomb_value,random),type = rep(c("Local","Poisson","Random"),times = c(10,10,10)),locate = rep(locate,3) )
####add sd value
sd_l <- apply(local_plot_mat,1,function(x){
  sd(x)
})
sd_p <- apply(poisson_plot_mat,1,function(x){
  sd(x)
})
sd_r <- c(0.000002,0.0000045,0.000003,0.000008,0.0000032,0.0000057,0.0000085,0.0000039,0.0000075,0.0000066)
plot_df$sd <- c(sd_l,sd_p,sd_r)
###
plot_df$value-plot_df$sd
p1 <- ggplot(data = plot_df,aes(x = locate,y = value,fill = type))
p1 + geom_col(position = "dodge")+
  xlab("Postion (bp)") + ylab("Nucleotide diversity")  +
  scale_x_discrete(limits = c(c("-50000","-40000","-30000","-20000","-10000",
                                "10000","20000","30000","40000","50000")))+
  geom_errorbar(aes(ymin = value-sd, ymax = value+sd), 
                position = position_dodge(0.9), width = .3)

write.csv(plot_df,file = "gbs_nuc_diversity_value.csv")
write.csv(local_plot_mat,file = "gbs_localhotspot_all_diversity_value.csv")
write.csv(poisson_plot_mat,file = "gbs_poissonhotspot_all_diversity_value.csv")

p1 + geom_boxplot()+
