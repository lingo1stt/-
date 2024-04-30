library(vcfR)
library(seqinr)
library(dplyr)


library(tidyr)
###data 前處理
for (q in 1:12) {
  dat <- read.csv(paste0("C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation result/chr",q,"/final_sample_plus_parent.chr",q,"_imputed_chunks.txt"),sep = " ")
  p <- strsplit(dat$genotypes,"")
  for (i in 1:length(p)) {
    for (j in 1:length(p[[1]])) {
      dat[i,j+2] <- p[[i]][j]
    }
  }
  row_name <- dat$X.
  rownames(dat) <- row_name
  dat <- dat[,-1:-2]  %>% as.data.frame() %>% t()
  write.csv(dat,file = paste0("C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation result/chr",q,"/",q,"_snp matrix.csv"))
  assign(paste0("snp_",q),as.data.frame(dat) ,envir = .GlobalEnv)
}


############叫出bkpt
for (i in 1:12) {
  temp <- scan(paste0("C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation result/chr",i,"/breakpoint_chr",i,".csv"),sep = "\n",what = "")
  a <- as.numeric(which(temp == "$` total without 3&24`"))
  b <- strsplit(temp[a:length(temp)],"[[:space:]]+")
  for (j in 1:length(b)) {
    b[[j]] <- as.numeric(b[[j]])
  }
  c <- na.omit(unlist(b)) %>% round(0)
  assign(paste0("bkpt_pos_chr",i),c,envir = .GlobalEnv)
}

data <- read.fasta("C:/Users/lingo1st/Dropbox/碩論/np ir64/REF_genome.fa/OsNB1.0.fa")

#########建立熱點區間
local_hotspot <- c(20522947, 20537617,20545769,
                   20556928,20672614,27359627,
                   27362882,27425162,27435631,
                   27526405,2107364,2888125,
                   2925595,19662508,19680032,
                   19849482)
local_hotspot_s <- c()
local_chr <-c(5,5,5,5,5,5,5,5,5,5,7,12,12,12,12,12)
poisson_hotspot <- c(500001, 1200001, 4700001,9000001,
                     26400001, 27700001,28900001,32100001,
                     35100001, 40400001, 42200001,300001,
                     900001,5100001,8100001,8500001,
                     21300001,28100001,29400001,30700001,
                     34100001,35200001,3800001,14400001,17100001,27700001,28500001,
                   32200001,32600001,35900001,19900001,26500001,29800001,35300001,20500001,27300001,27400001,29500001,6600001,27800001,
                     30400001,900001,1200001, 1500001,2100001,2200001,2700001, 5700001,
                     24700001,27900001,700001,800001, 2200001,20600001,24600001,25300001,
                     19100001,11900001, 16700001,100001,2700001, 5500001,21200001, 23200001,
                     25100001,2800001, 19600001)
poisson_chr <- c(rep(1:12,times = c(11,11,8,4,4,3,9,6,1,2,6,2)))

#poisson_hotspot <- c(28900001,900001,28100001,3800001,19900001,20500001,1500001,800001,21200001,19600001)
#poisson_chr <- c(1,2,2,3,4,5,7,8,11,12)

hotspot_interval <- data.frame(start =c(local_hotspot-50000,poisson_hotspot),
                               end = c(local_hotspot+50000,poisson_hotspot+100000),
                               chr = c(local_chr,poisson_chr),
                               type = c(rep("local",times = 16),rep("poisson",times = 67)))
write.csv(hotspot_interval,row.names = F,file ="C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation result/all_hotspot_interval.csv" )
if(length(which(apply(hotspot_interval %>% filter(chr == chr_check),1,function(x){
  between(position,x[1],x[2])})))==0 )
position = 8
####移除環境變數
all_vars <- ls()

# 使用正則表達式找出以 "data" 和 "snp" 開頭的變數
vars_to_keep <- grep("^data|^snp|^bkpt|^nuc", all_vars, value = TRUE)

# 移除除了 "data" 和 "snp" 開頭的其他變數
vars_to_remove <- setdiff(all_vars, vars_to_keep)
rm(list = vars_to_remove)

###################   fas_generate function  ########################################

fas_generate <- function(chr,one_side_length ,window_size ,sliding = 0,output_path){
  ##前置data
  dat <- read.vcfR(paste0("C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data/final_sample_plus_parent.chr",chr,".vcf"))
  dat_fix <- getFIX(dat) %>% as.data.frame()
  dat_gt <- as.data.frame(dat@gt)
  name <- c(colnames(dat@gt)[-1])
  dat_fix$CHROM <- as.numeric(dat_fix$CHROM)
  ###簡化版snp matrix (imputation後修正過的genotype) 
  true_geno <- get(paste0("snp_",chr))
  geno_pos <- colnames(true_geno) 
  num_geno_pos <- as.numeric(substr(geno_pos, start = 5, stop = nchar(geno_pos)) )
  # 列出路徑下的所有檔案
  folder_path <-paste0(output_path,"/",one_side_length,"_",window_size,"_",sliding,"_revised")
  path <- paste0(folder_path,"/chr",chr)
  
  files <- list.files(path)
  # 取得檔案數量
  num_files <- length(files)
  if (num_files ==0){
    #save_file <- readline(prompt = "there is no fasta file,Do you want to create files? (y/n): ")
    #if (tolower(save_file) == "y"){
    
    print("creating...")
    ###建立檔案位置
    dir.create(folder_path)
    dir.create(paste0(folder_path,"/chr",chr))
    #####前置data
    


    ##
    bkpt_pos <- get(paste0("bkpt_pos_chr",chr))
    for (ii in 1:length(bkpt_pos)) {
      print(ii)
      fas_dat <- list()
      ##要尋找區間內存在的snp
      dat_fix$POS <- as.numeric(dat_fix$POS)
      snp_pos_within <- dat_fix$POS[dat_fix$POS < bkpt_pos[ii]+one_side_length & dat_fix$POS > bkpt_pos[ii]-one_side_length ] %>% as.numeric()
      if(length(snp_pos_within)==0){next}
      
      for (jj in 1:length(name)) {
        
        ##先幫每個個體建立斷點兩側序列
        assign(name[jj],toupper(data[[chr]][(bkpt_pos[ii]-one_side_length):(bkpt_pos[ii]+one_side_length)]))
        
        ##針對每個點進行替換 (要替換成最接近true_geno位點的snp的基因型)
        temp <- get(name[jj])
        for (kk in 1:length(snp_pos_within)) {
          
          ##找出ref,alt allele
          ###找出在bkpt兩側最近的true_geno snp位置
          poss <- findInterval(bkpt_pos[ii],num_geno_pos)
          flanking <- c(num_geno_pos[poss],num_geno_pos[poss+1])
          change_pos <- as.numeric(which(dat_fix$CHROM == chr &  dat_fix$POS == snp_pos_within[kk]))
          ref <-  if (snp_pos_within[kk]<flanking[1]) dat_fix$REF[which(dat_fix$POS == flanking[1])] else dat_fix$REF[which(dat_fix$POS == flanking[2])]
          alt <- if (snp_pos_within[kk]<flanking[1]) dat_fix$ALT[which(dat_fix$POS == flanking[1])] else dat_fix$ALT[which(dat_fix$POS == flanking[2])]
          
          
          ###確認該點於該個體的genotype並且進行替換
          #將原始序列NA值替換成N
          if (is.na( temp[(snp_pos_within[kk]-(bkpt_pos[ii]-one_side_length))+1] )){temp[(snp_pos_within[kk]-(bkpt_pos[ii]-one_side_length))+1] = "N" }
          #將dat_gt在置換位點上的NA改成N
          if (is.na(dat_gt[,colnames(dat_gt)==name[jj]][change_pos])){dat_gt[,colnames(dat_gt)==name[jj]][change_pos] = "N" }
          
          
          #進行替換
          if(snp_pos_within[kk]<flanking[1]){
            if (true_geno[jj,which(num_geno_pos == flanking[1])] =="A"){
              temp[snp_pos_within[kk]-(bkpt_pos[ii]-one_side_length)+1] <- ref
            }else{
              temp[snp_pos_within[kk]-(bkpt_pos[ii]-one_side_length)+1] <- alt
            }
            
          }else if (snp_pos_within[kk]>flanking[2]){
            if (true_geno[jj,which(num_geno_pos == flanking[2])] =="A"){
              temp[snp_pos_within[kk]-(bkpt_pos[ii]-one_side_length)+1] <- ref
            }else{
              temp[snp_pos_within[kk]-(bkpt_pos[ii]-one_side_length)+1] <- alt
            }
          }
          
        }
        fas_dat[[jj]] <- temp
      }
      names(fas_dat) <- name
      
      write.fasta(fas_dat,as.string = F,names = name,file.out = paste0(folder_path,"/chr",chr,"/chr",chr,"_bkpt_",ii,".fa"))
    }
  }
}


for (i in 2:12) {
  fas_generate(i,one_side_length = 6000,window_size = 1000,sliding = 0,output_path = "C:/Users/lingo1st/Dropbox/林冠瑜/nucloetide diversity")
  
}





#################################################
#nucleotide diversity fucntion version 3
#################################################

#將變成三大功能，第一為直接提供fasta檔，第二是提供vcf檔以及reference檔的路徑，並且提供指定位置，第三是根據提供的資料隨機抽取

nuc_diversity_v3 <- function(fas_data_path=NA,chr = NA,vcf_reffas_path=NA,vcf_name=NA,ref_data,position=NA,one_side_length,window_size,sliding = 0,output_path,random_pos = F,repeat_time = 1){
  ###警示
  if(missing(one_side_length) || missing(window_size) || missing(output_path) ){
    stop("window_size or one_side_length or outputpath參數遺失")
  }
  
  
  ##### function 1: 如果有提供fas資料，就直接根據該fasta檔進行計算，每個給定的file_path會產出一個dataframe
  if(!is.na(fas_data_path)){
    ####警示
    #if(is.na(chr)){
     # stop("must add chr parameter")
   # }
    cat("############################################","\n","you are now using the first method of the function","\n","calculating nucleotide diversity...")
    # 取得檔案數量
    folder_path <-fas_data_path
    # 列出路徑下的所有檔案
    files <- list.files(folder_path)
    num_files <- length(files)
    if(sliding == 0){
      plot_mat <- matrix(ncol = num_files,nrow = length(seq(1,(2*one_side_length)-sliding,by = window_size))  )
    }else{
      plot_mat <- matrix(ncol = num_files,nrow = length(seq(1,(2*one_side_length)-sliding,by = sliding))  )
    }
    for (q in 1:num_files) {
      fas <- read.fasta(paste0(folder_path,"/",files[q]))
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
        ##擷取每個區間的長度(1000bp)(temp_list)
        temp_list <- lapply(fas,FUN =  sub,pos = i)
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
          result <- sum(result)*(length(fas)/(length(fas)-1))
        }
        pi_diversity <- append(pi_diversity,result)
        
        
        
      }
      plot_mat[1:nrow(plot_mat),q] <- pi_diversity
      
    }
    
    assign(paste0("plot_mat_",chr),plot_mat,envir = .GlobalEnv)
    dir.create(paste0(output_path,"/diversity value"))
    #dir.create(paste0(output_path,"/diversity value_revised/",one_side_length,"_",window_size,"_",sliding,"_revised"))
    write.csv(plot_mat,file = paste0(output_path,"/diversity value/","diverisity value.csv"))
    return("jon done")
  }
  
  ##### function 2: 根據給定的vcf、ref fasta、chr以及position來計算該點的diversity值
  if( random_pos ==F & is.na(fas_data_path)  ){
    ##警示
    if(is.na(vcf_reffas_path) || is.na(vcf_name) || is.null(ref_data) || is.na(position) ){
      stop("請檢查五個參數是否遺失:chr,vcf_name,vcf_reffas_path,ref_data,position")
    }
    cat("############################################","\n","you are now using the second method of the function","\n","calculating nucleotide diversity...")
    
    ##前置data
    data  <-ref_data
    dat <- read.vcfR(paste0(vcf_reffas_path,"/",vcf_name,chr,".vcf"))
    dat_fix <- getFIX(dat) %>% as.data.frame()
    dat_gt <- as.data.frame(dat@gt)
    name <- c(colnames(dat@gt)[-1])
    dat_fix$CHROM <- as.numeric(dat_fix$CHROM)
    ###簡化版snp matrix (imputation後修正過的genotype)
    true_geno <- get(paste0("snp_",chr))
    geno_pos <- colnames(true_geno) 
    num_geno_pos <- as.numeric(substr(geno_pos, start = 5, stop = nchar(geno_pos)) )
    ###
    fas_dat <- list()
    ##要尋找區間內存在的snp
    dat_fix$POS <- as.numeric(dat_fix$POS)
    snp_pos_within <- dat_fix$POS[dat_fix$POS < position+one_side_length & dat_fix$POS > position-one_side_length ] %>% as.numeric()
    if(length(snp_pos_within)==0){
      cat("\n","No snp within the interval","\n")
      step <- if (sliding == 0) window_size else sliding
      
      assign(paste0("diversity_",position),c(rep(0,length( seq(1,(2*one_side_length)-window_size+1,by = step)))),envir = .GlobalEnv)
      return(cat("nucleotide diversity value within the interval is assigned to ",paste0("diversity_",position),"\n"))
      }
    
    for (jj in 1:length(name)) {
      
      ##先幫每個個體建立斷點兩側序列
      assign(name[jj],toupper(data[[chr]][(position-one_side_length):(position+one_side_length)]))
      
      ##針對每個點進行替換 (要替換成最接近true_geno位點的snp的基因型)
      temp <- get(name[jj])
      for (kk in 1:length(snp_pos_within)) {
        
        ##找出ref,alt allele
        ###找出在bkpt兩側最近的true_geno snp位置
        poss <- findInterval(position,num_geno_pos)
        flanking <- c(num_geno_pos[poss],num_geno_pos[poss+1])
        change_pos <- as.numeric(which(dat_fix$CHROM == chr &  dat_fix$POS == snp_pos_within[kk]))
        ref <-  if (snp_pos_within[kk]<flanking[1]) dat_fix$REF[which(dat_fix$POS == flanking[1])] else dat_fix$REF[which(dat_fix$POS == flanking[2])]
        alt <- if (snp_pos_within[kk]<flanking[1]) dat_fix$ALT[which(dat_fix$POS == flanking[1])] else dat_fix$ALT[which(dat_fix$POS == flanking[2])]
        
        
        ###確認該點於該個體的genotype並且進行替換
        #將原始序列NA值替換成N
        if (is.na( temp[(snp_pos_within[kk]-(position-one_side_length))+1] )){temp[(snp_pos_within[kk]-(position-one_side_length))+1] = "N" }
        #將dat_gt在置換位點上的NA改成N
        if (is.na(dat_gt[,colnames(dat_gt)==name[jj]][change_pos])){dat_gt[,colnames(dat_gt)==name[jj]][change_pos] = "N" }
        
        
        #進行替換
        if(snp_pos_within[kk]<flanking[1]){
          if (true_geno[jj,which(num_geno_pos == flanking[1])] =="A"){
            temp[snp_pos_within[kk]-(position-one_side_length)+1] <- ref
          }else{
            temp[snp_pos_within[kk]-(position-one_side_length)+1] <- alt
          }
          
        }else if (snp_pos_within[kk]>flanking[2]){
          if (true_geno[jj,which(num_geno_pos == flanking[2])] =="A"){
            temp[snp_pos_within[kk]-(position-one_side_length)+1] <- ref
          }else{
            temp[snp_pos_within[kk]-(position-one_side_length)+1] <- alt
          }
        }
        
      }
      fas_dat[[jj]] <- temp
    }
    names(fas_dat) <- name
    #####建立一個fas_dat的list，包含所有個體的序列，使用該檔案進行計算
    
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
    assign(paste0("diversity_",position),pi_diversity,envir = .GlobalEnv)
    cat("#######################################################","\n")
    cat("nucleotide diversity value within the interval is assigned to ",paste0("diversity_",position),"-->","\n",pi_diversity,"\n")
    #ans <- readline(prompt = "do you want to save your fasta file? (y/n): ")
    #if (tolower(ans) =="y" ){
      #write.fasta(fas_dat,as.string = F,names = name,file.out =paste0(output_path,"/","diversity_",position,".fasta"))
      cat("file saved in ",output_path)
    #}else{return("end of the function")}
  }
  
  ##### function 3: 隨機在染色體上找點計算
  if( random_pos == T){
    ##警示
    if(is.na(vcf_reffas_path) || is.na(vcf_name) || is.null(ref_data)  ){
      stop("請檢查4個主要參數是否遺失:vcf_name,vcf_reffas_path,ref_fa_name")
    }
    cat("############################################","\n","you are now using the third method of the function","\n","calculating nucleotide diversity...")
    
    ##前置data
    data <-ref_data

    if(!is.na(chr)){
      dat <- read.vcfR(paste0(vcf_reffas_path,"/",vcf_name,chr,".vcf"))
      dat_fix <- getFIX(dat) %>% as.data.frame()
      dat_gt <- as.data.frame(dat@gt)
      name <- c(colnames(dat@gt)[-1])
      dat_fix$CHROM <- as.numeric(dat_fix$CHROM)
      bkt <- get(paste0("bkpt_pos_chr",chr))
      
      ###簡化版snp matrix (imputation後修正過的genotype)
      true_geno <- get(paste0("snp_",chr))
      geno_pos <- colnames(true_geno) 
      num_geno_pos <- as.numeric(substr(geno_pos, start = 5, stop = nchar(geno_pos)) )
      
    }
    path <- vcf_reffas_path
    # 列出路徑下的所有檔案
    files <- list.files(path)
    # 取得檔案數量作為染色體數量
    vcf_num <- as.numeric(length(files[grep(paste0("^",vcf_name),files)]) )
    
    #####R建立新環境以存放隨機抽取中已讀取過的vCF(命名為vcf_1~12)
    read_files_env <- new.env()
    ###
    step <- if (sliding == 0) window_size else sliding
    random_pi <- c()
    random_pi_pos <- c(rep( seq(1,(2*one_side_length)-window_size+1,by =step ),time = repeat_time))
    random_pi_chr <- c()
    count =0
    
    ###
    while(count < repeat_time){
      
      count = count +1
      chr_check <- if (is.na(chr)) NA else chr
      if(is.na(chr_check)){
        chr_check = sample(1:vcf_num,1)
        if( !exists(paste0("vcf_",chr_check),envir =read_files_env ) ){
          print("##################NOT USED YET ##################################")
          dat <- read.vcfR(paste0(vcf_reffas_path,"/",vcf_name,chr_check,".vcf"))
          assign(paste0("vcf_",chr_check),dat,envir = read_files_env)
          dat_fix <- getFIX(dat) %>% as.data.frame()
          dat_gt <- as.data.frame(dat@gt)
          name <- c(colnames(dat@gt)[-1])
          dat_fix$CHROM <- as.numeric(dat_fix$CHROM)
          bkt <- get(paste0("bkpt_pos_chr",chr_check))
          ###簡化版snp matrix (imputation後修正過的genotype)
          true_geno <- get(paste0("snp_",chr_check))
          geno_pos <- colnames(true_geno) 
          num_geno_pos <- as.numeric(substr(geno_pos, start = 5, stop = nchar(geno_pos)) )
          
        }else{
          print("#################### ALREADY USED FILE ###############################")
          dat <- read_files_env[[paste0("vcf_",chr_check)]]
          dat_fix <- getFIX(dat) %>% as.data.frame()
          dat_gt <- as.data.frame(dat@gt)
          name <- c(colnames(dat@gt)[-1])
          dat_fix$CHROM <- as.numeric(dat_fix$CHROM)
          bkt <- get(paste0("bkpt_pos_chr",chr_check))
          ###簡化版snp matrix (imputation後修正過的genotype)
          true_geno <- get(paste0("snp_",chr_check))
          geno_pos <- colnames(true_geno) 
          num_geno_pos <- as.numeric(substr(geno_pos, start = 5, stop = nchar(geno_pos)) )
        }
        
        
        
      }
      random_pi_chr <- c(random_pi_chr,chr_check)
      ##如果有breakpoint存在在隨機找的區間，就重新亂數一個數字當作中心
      repeat {
        position = as.numeric(sample(1000:length(data[[chr_check]]),1,replace = F))
        
        #如果position的位置沒有出現在熱點區間內，就停止抽取
        if (  length(which(apply(hotspot_interval %>% filter(chr == chr_check),1,function(x){
          between(position,x[1],x[2])})))==0 ) {
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
        count = count-1
        random_pi_chr <- head(random_pi_chr,-1)
        next
      }
      
      for (jj in 1:length(name)) {
        print(jj)
        
        ##先幫每個個體建立斷點兩側序列
        assign(name[jj],toupper(data[[chr_check]][(position-one_side_length):(position+one_side_length)]))
        
        ##針對每個點進行替換 (要替換成最接近true_geno位點的snp的基因型)
        temp <- get(name[jj])
        for (kk in 1:length(snp_pos_within)) {
          
          ##找出ref,alt allele
          ###找出在bkpt兩側最近的true_geno snp位置
          poss <- findInterval(position,num_geno_pos)
          flanking <- c(num_geno_pos[poss],num_geno_pos[poss+1])
          change_pos <- as.numeric(which(dat_fix$CHROM == chr_check &  dat_fix$POS == snp_pos_within[kk]))
          ref <-  if (snp_pos_within[kk]<flanking[1]) dat_fix$REF[which(dat_fix$POS == flanking[1])] else dat_fix$REF[which(dat_fix$POS == flanking[2])]
          alt <- if (snp_pos_within[kk]<flanking[1]) dat_fix$ALT[which(dat_fix$POS == flanking[1])] else dat_fix$ALT[which(dat_fix$POS == flanking[2])]
          
          dat_fix %>% filter(POS == snp_pos_within[kk])
          ###確認該點於該個體的genotype並且進行替換
          #將原始序列NA值替換成N
          if (is.na( temp[(snp_pos_within[kk]-(position-one_side_length))+1] )){temp[(snp_pos_within[kk]-(position-one_side_length))+1] = "N" }
          #將dat_gt在置換位點上的NA改成N
          if (is.na(dat_gt[,colnames(dat_gt)==name[jj]][change_pos])){dat_gt[,colnames(dat_gt)==name[jj]][change_pos] = "N" }
          
          
          #進行替換
          if(snp_pos_within[kk]<flanking[1]){
            if (true_geno[jj,which(num_geno_pos == flanking[1])] =="A"){
              temp[snp_pos_within[kk]-(position-one_side_length)+1] <- ref
            }else{
              temp[snp_pos_within[kk]-(position-one_side_length)+1] <- alt
            }
            
          }else if (snp_pos_within[kk]>flanking[2]){
            if (true_geno[jj,which(num_geno_pos == flanking[2])] =="A"){
              temp[snp_pos_within[kk]-(position-one_side_length)+1] <- ref
            }else{
              temp[snp_pos_within[kk]-(position-one_side_length)+1] <- alt
            }
          }
          
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
        random_pi <- append(random_pi,result)
        
      }
    }
    
    
    assign(paste0("diversity_","value"),random_pi,envir = .GlobalEnv)
    assign(paste0("diversity_","pos"),random_pi_pos,envir = .GlobalEnv)
    assign(paste0("diversity_","chr"),random_pi_chr,envir = .GlobalEnv)
    
    
    cat("#######################################################","\n")
    return(cat("nucleotide diversity value within the interval is assigned to ",paste0("diversity_value, diversity_chr and diversity_pos"),"\n") )
    
    
    
  }
}


####################################### 測試 #############################################

vcf_reffas_path="C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data"
vcf_name = "final_sample_plus_parent.chr"

### method 1
for (i in 2:12) {
  nuc_diversity_v3(fas_data_path = paste0("C:/Users/lingo1st/Dropbox/林冠瑜/nucloetide diversity/6000_1000_0_revised/chr",i),
                   chr = i,
                   output_path ="C:/Users/lingo1st/Dropbox/林冠瑜/nucloetide diversity" )
}

nuc_diversity_v3(fas_data_path = "C:/Users/lingo1st/Dropbox/林冠瑜/nucloetide diversity/poisson_dis fasta(100kb)/p value_0.01",
                 window_size = 10000,one_side_length = 50000,
                 output_path ="C:/Users/lingo1st/Dropbox/林冠瑜/nucloetide diversity" )
###method 2
nuc_diversity_v3(vcf_reffas_path = "C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data",
                 ref_data  = data,
                 vcf_name = "final_sample_plus_parent.chr",
                 position = 20523254,
                 chr = 5,window_size = 1000,one_side_length = 6000,
                 output_path ="C:/Users/lingo1st/Dropbox/林冠瑜/nucloetide diversity" )
#method 3
nuc_diversity_v3(vcf_reffas_path = "C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data",
                 ref_data  = data,
                 vcf_name = "final_sample_plus_parent.chr",
                 window_size = 1000,one_side_length = 5000,
                 output_path ="C:/Users/lingo1st/Dropbox/林冠瑜/nucloetide diversity",random_pos = T,repeat_time = 3 )


################# 找出local recombination rate 熱點的區間進行計算 #########################
local_hotspot <- c(20522947, 20537617,20545769,
                   20556928,20672614,27359627,
                   27362882,27425162,27435631,
                   27526405,2107364,2888125,
                   2925595,19662508,19680032,
                   19849482)
local_chr <-c(5,5,5,5,5,5,5,5,5,5,7,12,12,12,12,12)
window_size = 1000
one_side_length = 5000
output_path ="C:/Users/lingo1st/Dropbox/林冠瑜/nucloetide diversity"
vcf_reffas_path = "C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data"
ref_data  = data
i=11
for (i in 1:16) {
  nuc_diversity_v3(window_size = 10000,one_side_length = 50000,
                   output_path = output_path,vcf_reffas_path = vcf_reffas_path,
                   ref_data = data,vcf_name = "final_sample_plus_parent.chr",
                   position = local_hotspot[i],chr =local_chr[i] )
}
###將數值從環境中取出
local_hotspot_value <- c()
for (i in 1:16) {
  p <- get(paste0("diversity_",local_hotspot[i]))
  local_hotspot_value <- c(local_hotspot_value,p)
  
}
##################### 找出 physical distance 熱點的區間進行計算 #########################
poisson_hotspot <- c(28900001,900001,28100001,3800001,19900001,20500001,1500001,800001,21200001,19600001)
poisson_chr <- c(1,2,2,3,4,5,7,8,11,12)
poisson_hotspot <- poisson_hotspot +50000
window_size = 10000
one_side_length = 50000
output_path ="C:/Users/lingo1st/Dropbox/林冠瑜/nucloetide diversity"
vcf_reffas_path = "C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data"
ref_data  = data

for (i in 1:10) {
  nuc_diversity_v3(window_size = 10000,one_side_length = 50000,
                   output_path = output_path,vcf_reffas_path = vcf_reffas_path,
                   ref_data = data,vcf_name = "final_sample_plus_parent.chr",
                   position = poisson_hotspot[i],chr =poisson_chr[i] )
}
###將數值從環境中取出
poisson_hotspot_value <- c()
for (i in 1:16) {
  p <- get(paste0("diversity_",poisson_hotspot[i]))
  poisson_hotspot_value <- c(poisson_hotspot_value,p)
  
}
############################### random sampling #######################
nuc_diversity_v3(vcf_reffas_path = "C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data",
                 ref_data  = data,
                 vcf_name = "final_sample_plus_parent.chr",
                 window_size = 10000,one_side_length = 50000,
                 output_path ="C:/Users/lingo1st/Dropbox/林冠瑜/nucloetide diversity",random_pos = T,repeat_time = 20 )



random_plot_df <- data.frame(chr = rep(diversity_chr,each = 10),value = diversity_value)
write.csv(random_plot_df,file ="C:/Users/lingo1st/Dropbox/林冠瑜/nucloetide diversity/random sampling.csv")
### 計算95%ci 
l.model <- lm(diversity_value ~ 1)
confint(l.model,level = 0.95)

#####################將三區域的數值合併成dataframe################################
locate <- c("-50000","-40000","-30000","-20000","-10000",
            "10000","20000","30000","40000","50000")

final_plot_df <- data.frame(value =c(local_hotspot_value,poisson_hotspot_value,diversity_value),
                            location = rep(locate,times = (16+67+20)),
                            type = rep(c("Local recombination rate hotspot","Poisson-determined hotspot","Random-selected positions"),each = c(160,670,200)))


###################### average col_plot #############################
l_avg <- c()
p_avg <- c()
r_avg <- c()
for (i in 1:10) {
  p <- final_plot_df %>% filter(location == locate[i] & group =="Local recombination rate hotspots" ) %>% select("value") %>% unlist() %>% mean()
  l_avg <- c(l_avg,p)
  
  q <- final_plot_df %>% filter(location == locate[i]& group =="Poisson method-identified hotspots") %>% select("value") %>% unlist() %>% mean()
  p_avg <- c(p_avg,q)
  
  r<- final_plot_df %>% filter(location == locate[i]& group =="Random-selected regions") %>% select("value") %>% unlist() %>% mean()
  r_avg <-c(r_avg,r)
}
avg_plot_df <- data.frame(position = rep(locate,3),value = c(l_avg,p_avg,r_avg), Type = rep(c("Local recombination rate-\ndetermined hotspot","Poisson-determined hotspot","Random-selected positions"),each = 10) )
p2 <-ggplot(data = avg_plot_df,aes(x = position,y = value,fill = Type)) 
p2+geom_col(position = "dodge") +
  xlab("Postion (bp)") + ylab("Nucleotide diversity")  +
  scale_x_discrete(limits = c(c("-50000","-40000","-30000","-20000","-10000",
                                "10000","20000","30000","40000","50000")))



##################只取5,7,12 進行作圖 #################
p_avg_5712 <- c()
which(poisson_chr == 5 | poisson_chr == 7| poisson_chr == 12 )

for (i in 1:10) {
  
  q <- final_plot_df %>% filter(location == locate[i]& group =="Poisson method-identified hotspots") %>% select("value") 
  q <- q[c(which(poisson_chr == 5 | poisson_chr == 7| poisson_chr == 12 )),]
  q <- mean(q)
  p_avg_5712 <- c(p_avg_5712,q)
  
}
avg_plot_df_5712 <- data.frame(position = rep(locate,3),value = c(l_avg,p_avg_5712,r_avg), Type = rep(c("Local recombination rate-\ndetermined hotspot","Poisson-determined hotspot","Random-selected positions"),each = 10) )

p2 <-ggplot(data = avg_plot_df_5712,aes(x = position,y = value,fill = Type)) 
p2+geom_col(position = "dodge") +
  xlab("Postion (bp)") + ylab("Nucleotide diversity")  +
  scale_x_discrete(limits = c(c("-50000","-40000","-30000","-20000","-10000",
                                "10000","20000","30000","40000","50000")))

##############################進行常態性檢定以及變異數同質性檢定####################################
#################### 每個點分開做
p_nor_1 <- poisson_plot_df %>% filter(location == -5000) %>% select("value") %>% unlist()
l_nor_1 <-local_plot_df %>% filter(location == -5000) %>% select("value") %>% unlist()

print(shapiro.test(l_nor_1))
u <- unique(combined_plot_df$location)
for (i in 1:10) {
  p <- local_plot_df %>% filter(location == u[i]) %>% select("value") %>% unlist()
  print(shapiro.test(p))
}

###變異數同質性檢定(levene.test)
instball.packages("car")
library(car)
i=1
for (i in 1:10) {
  p1 <- local_plot_df %>% filter(location == u[i]) %>% select("value") %>% unlist()
  p2 <- poisson_plot_df %>% filter(location == u[i]) %>% select("value") %>% unlist()
  val <- c(p1,p2)
  fac <- as.factor(c( rep(c(1,2),times = c(length(p1),length(p2))) ))
  var_df <- data.frame(fac = fac,v = val)
  print(leveneTest(v~fac,data = var_df))
}

####welch anova (變異數不同質)
install.packages("userfriendlyscience")
library(userfriendlyscience)
oneway.test(data = combined_plot_df,value~group,var.equal = F)
games.howell <- function(grp, obs) {
  
  #Create combinations
  combs <- combn(unique(grp), 2)
  
  # Statistics that will be used throughout the calculations:
  # n = sample size of each group
  # groups = number of groups in data
  # Mean = means of each group sample
  # std = variance of each group sample
  n <- tapply(obs, grp, length)
  groups <- length(tapply(obs, grp, length))
  Mean <- tapply(obs, grp, mean)
  std <- tapply(obs, grp, var)
  
  statistics <- lapply(1:ncol(combs), function(x) {
    
    mean.diff <- Mean[combs[2,x]] - Mean[combs[1,x]]
    
    #t-values
    t <- abs(Mean[combs[1,x]] - Mean[combs[2,x]]) / sqrt((std[combs[1,x]] / n[combs[1,x]]) + (std[combs[2,x]] / n[combs[2,x]]))
    
    # Degrees of Freedom
    df <- (std[combs[1,x]] / n[combs[1,x]] + std[combs[2,x]] / n[combs[2,x]])^2 / # Numerator Degrees of Freedom
      ((std[combs[1,x]] / n[combs[1,x]])^2 / (n[combs[1,x]] - 1) + # Part 1 of Denominator Degrees of Freedom 
         (std[combs[2,x]] / n[combs[2,x]])^2 / (n[combs[2,x]] - 1)) # Part 2 of Denominator Degrees of Freedom
    
    #p-values
    p <- ptukey(t * sqrt(2), groups, df, lower.tail = FALSE)
    
    # Sigma standard error
    se <- sqrt(0.5 * (std[combs[1,x]] / n[combs[1,x]] + std[combs[2,x]] / n[combs[2,x]]))
    
    # Upper Confidence Limit
    upper.conf <- lapply(1:ncol(combs), function(x) {
      mean.diff + qtukey(p = 0.95, nmeans = groups, df = df) * se
    })[[1]]
    
    # Lower Confidence Limit
    lower.conf <- lapply(1:ncol(combs), function(x) {
      mean.diff - qtukey(p = 0.95, nmeans = groups, df = df) * se
    })[[1]]
    
    # Group Combinations
    grp.comb <- paste(combs[1,x], ':', combs[2,x])
    
    # Collect all statistics into list
    stats <- list(grp.comb, mean.diff, se, t, df, p, upper.conf, lower.conf)
  })
  
  # Unlist statistics collected earlier
  stats.unlisted <- lapply(statistics, function(x) {
    unlist(x)
  })
  
  # Create dataframe from flattened list
  results <- data.frame(matrix(unlist(stats.unlisted), nrow = length(stats.unlisted), byrow=TRUE))
  
  # Select columns set as factors that should be numeric and change with as.numeric
  results[c(2, 3:ncol(results))] <- round(as.numeric(as.matrix(results[c(2, 3:ncol(results))])), digits = 5)
  
  # Rename data frame columns
  colnames(results) <- c('groups', 'Mean Difference', 'Standard Error', 't', 'df', 'p', 'upper limit', 'lower limit')
  
  return(results)
}
games.howell(grp = combined_plot_df$group,obs = combined_plot_df$value)

### normal one way anova
TukeyHSD(aov(value ~ group, data = combined_plot_df), 'group')

#############################合併做一次檢定

combined_plot_df$value[1:670]
diversity_value

print(shapiro.test(local_hotspot_value_2))
shapiro.test(combined_plot_df$value[1:670])
shapiro.test(diversity_value)
###變異數同質性檢定(levene.test)
combined_plot_df$value[671:830] <-  local_hotspot_value_2

leveneTest(value~group,data = combined_plot_df)

###welch anova
oneway.test(data = combined_plot_df,value~group,var.equal = F)
games.howell(grp = combined_plot_df$group,obs = combined_plot_df$value)



################針對5 7 12的數值進行檢定
final_plot_df <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/nucloetide diversity/綜合三區域原始數據/hotpsot_random conbine.csv")
local_value <- final_plot_df %>% filter(group == "Local recombination rate hotspots") %>% select("value") %>% unlist()
temp <- final_plot_df%>% filter(group == "Random-selected regions" | group == "Local recombination rate hotspots" )
p <- c(10* which(poisson_chr == 5 | poisson_chr == 7| poisson_chr == 12 ) -9)

t <- lapply(p,function(x){
  final_plot_df %>% slice(x:(x+9))
  
})
final_5712 <- do.call(rbind,t)
final_5712 <- rbind(final_5712,temp)
oneway.test(data = final_5712,value~group,var.equal = F)
games.howell(grp = final_5712$group,obs = final_5712$value)


