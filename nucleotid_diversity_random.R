library(vcfR)
library(seqinr)
###########################################################################################################
###snp matrix
abh_geno <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_abhgenotype.csv",header = T)

############叫出bkpt
breakpoint <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_breakpoint.csv",header = T)

data <- read.fasta("C:/Users/lingo1st/Dropbox/碩論/np ir64/REF_genome.fa/OsNB1.0.fa")

#########建立熱點區間
local_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/local_recomb/local_hotspot_spar_0.1.csv",header = T)
poisson_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_adjusted_poisson_hotspot.csv",header = T)
names(local_hotspot)[names(local_hotspot) == "map"] <- "chr"
local_hotspot$chr <- gsub("chr","",local_hotspot$chr)
###綜合成hotspot_interval
hotspot_interval <- rbind(local_hotspot[,c(3,6,7)],poisson_hotspot[,c(2,5,6)])
hotspot_interval$type <- rep(c("local","poisson"),times = c(145,152))
hotspot_interval$chr <- as.numeric(hotspot_interval$chr)
hotspot_interval$start <- as.numeric(hotspot_interval$start)
hotspot_interval$end <- as.numeric(hotspot_interval$end)

chr = 1
one_side_length = 50000
output_path
vcf_reffas_path= "C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data"
vcf_name="gbs_r.chr"
position = 1050001
dir.create("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/hotspot_fasta_2")
setwd("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/hotspot_fasta_2")
###################################################################################################################
random <- function(ram_chr,one_side_length ,vcf_name =NA,ref_data, vcf_reffas_path = NA){
  
  #前置data
  data  <- ref_data
  dat <- read.vcfR(paste0(vcf_reffas_path,"/",vcf_name,ram_chr,".vcf"),verbose = F)
  dat_fix <- getFIX(dat) %>% as.data.frame()
  dat_gt <- as.data.frame(dat@gt %>% t()) 
  name <- c(colnames(dat@gt)[-1])
  dat_fix$CHROM <- as.numeric(dat_fix$CHROM)
  ###簡化版snp matrix (imputation後修正過的genotype)
  cols <- colnames(abh_geno)[abh_geno[1, ] == ram_chr]
  cols = na.omit(cols)
  true_geno <-abh_geno[,cols]
  true_geno <- true_geno[-1, ]
  ###
  geno_pos <- colnames(true_geno) 
  num_geno_pos <- as.numeric(sapply(strsplit(geno_pos,"_"),function(x){x[2]}) )
  ###
  fas_dat <- list()
  
  ##如果有breakpoint存在在隨機找的區間，就重新亂數一個數字當作中心
    repeat {
      position = as.numeric(sample(1000:length(data[[ram_chr]]),1,replace = F))
      
      #如果position的位置沒有出現在熱點區間內，就停止抽取
      if (  length(which(apply(hotspot_interval %>% filter(chr == ram_chr),1,function(x){
        between(position,as.numeric(x[2]),as.numeric(x[3]))})))==0 ) {
        cat("random site:",position,"\n")
        break
      } else {next}
    }

    ##針對position位置建立檔案及計算nucleotide diversity 
    
    ##產生fasta file
    
    fas_dat <- list()
    ##要尋找區間內存在的snp
    dat_fix$POS <- as.numeric(dat_fix$POS)
    snp_pos_within <- dat_fix$POS[dat_fix$POS < position+one_side_length & dat_fix$POS > position-one_side_length ] %>% as.numeric()
    if(length(snp_pos_within)==0){
      cat("\n","No snp within the interval","\n")
      return("ddd")
      
    }
    #個體迴圈
    for (jj in 1:length(name)) {
      assign(name[jj],toupper(data[[ram_chr]][(position-one_side_length):(position+one_side_length)]))
      ##先幫每個個體建立斷點兩側序列
      temp <- toupper(data[[ram_chr]][(position-one_side_length):(position+one_side_length)])
      ##針對每個點進行替換 (要替換成最接近true_geno位點的snp的基因型)
      for (kk in 1:length(snp_pos_within)) {
        #先確認snp的基因型以及其周圍tru_geno的基因型
        snp_vcf_type = dat_gt[jj+1,which(dat_fix$POS==snp_pos_within[kk] & dat_fix$CHROM == ram_chr)]
        vcf_alt <- dat_fix$ALT[which(dat_fix$POS==snp_pos_within[kk] & dat_fix$CHROM == ram_chr)]
        vcf_ref <- dat_fix$REF[which(dat_fix$POS==snp_pos_within[kk] & dat_fix$CHROM == ram_chr)]
        true_geno_left <- true_geno[jj,findInterval(snp_pos_within[kk],num_geno_pos)]
        true_geno_right <- true_geno [jj,findInterval(snp_pos_within[kk] +1,num_geno_pos)]
        ####找出snp_within位點的major allele
        major_allele <- dat_gt[-1,which(dat_fix$POS==snp_pos_within[kk] & dat_fix$CHROM == ram_chr)]
        major_allele <- names( which.max(table(major_allele)) )
        major_allele <- if (major_allele == "1/1") vcf_ref else vcf_alt
        ## 先透過原先vcf將該個體的核苷酸置換 (異質結合或是空值都使用major allele取代)
        temp[ (snp_pos_within[kk]-(position-one_side_length+1 )) ] <- case_when(
          snp_vcf_type == "1/1" ~ vcf_ref,
          snp_vcf_type == "0/0" ~ vcf_alt,
          TRUE ~ major_allele 
        )
        
        #接著根據truegeno中的旁邊基因型替換成應該有的核苷酸(由於temp序列只包含部分，因此其index= snp絕對位置- window size +1)
        temp[ (snp_pos_within[kk]-(position-one_side_length+1 )) ] <- case_when(
          true_geno_left == true_geno_right & true_geno_left == "A" ~ vcf_ref,
          true_geno_left == true_geno_right & true_geno_left == "B" ~ vcf_alt,
          TRUE ~ temp[ (snp_pos_within[kk]-(position-one_side_length ))+1 ] 
        )
        
        
      }
      fas_dat[[jj]] <- temp
    }
    names(fas_dat) <- name
    ######### step2: 建立一個fas_dat的list，包含所有個體的序列，使用該檔案進行計算
    
    ##前置function
    sub <- function(x,pos){
      return(x[pos:(pos+(window_size-1))])
    }
    
    freq <- function(y,ind){
      return(length(which(y == temp_list[[ind]])))
    }
    ###計算Pi
    pi_diversity <- c()
    step <- if (sliding == 0) window_size else sliding
    for (i in seq(1,(2*one_side_length)-window_size+1,by = step)) {
      ##擷取每個區間的長度(1000bp)(temp_list)
      temp_list <- lapply(fas_dat,FUN = sub,pos = i)
      uni <- unique(temp_list)
      ##比對
      freq_list_2 <- list()
      for (p in 1:length(uni)) {
        sub_list <- list()
        freq_list_2[[p]] <- sub_list
      }
      ##將每個個體的名字分配給所屬的不重複序列(freq_list_2) (freq_list為不重複序列本身)
      for (l in 1:24) {
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
        result <- sum(result)*(length(fas_dat)/(length(fas_dat)-1))
        
      }
      pi_diversity <- append(pi_diversity,result)
      
    }
  assign(paste0("diversity_","value"),pi_diversity,envir = .GlobalEnv)
  cat("#######################################################","\n")
  return(cat("nucleotide diversity value within the interval is assigned to ",paste0("diversity_value, diversity_chr and diversity_pos"),"\n") )
    }
    