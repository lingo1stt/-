####
#detach("package:vcfR",unload = T) # 解除library
library(vcfR)
library(seqinr)


###snp matrix
abh_geno <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_abhgenotype.csv",header = T)

############叫出bkpt
breakpoint <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_breakpoint.csv",header = T)

data <- read.fasta("C:/Users/lingo1st/Dropbox/碩論/np ir64/REF_genome.fa/OsNB1.0.fa")

#########建立熱點區間 (加入overlapped interval)
overlap_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_overlap_hotspot.csv",header = T)
poisson_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file//poisson hotspot/gbs_adjusted_poisson_hotspot.csv",header = T)
local_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/local_recomb_hotspot/gbs_nonoverlap_local_hotspot.csv",header = T)

#names(local_hotspot)[names(local_hotspot) == "map"] <- "chr"
#local_hotspot$chr <- gsub("chr","",local_hotspot$chr)
###綜合成hotspot_interval
hotspot_interval <- rbind(local_hotspot[,c(1,2,3)],poisson_hotspot[,c(2,5,6)],overlap_hotspot[,c(1,2,3)])
hotspot_interval$type <- rep(c("local","poisson","overlap"),times = c(25,152,21))
hotspot_interval$chr <- as.numeric(hotspot_interval$chr)
###################################################################################################################

chr = 1
one_side_length = 50000
output_path
vcf_reffas_path= "C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data"
vcf_name="gbs_r.chr"
position = 1050001
dir.create("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/hotspot_fasta_2")
setwd("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/hotspot_fasta_2")


fas_generate <- function(chr,position,one_side_length ,vcf_name =NA,ref_data, vcf_reffas_path = NA){
  
  ##前置data
  data  <- ref_data
  library(vcfR)
  dat <- read.vcfR(paste0(vcf_reffas_path,"/",vcf_name,chr,".vcf"),verbose = F)
  dat_fix <- getFIX(dat) %>% as.data.frame()
  dat_gt <- as.data.frame(dat@gt %>% t()) 
  name <- c(colnames(dat@gt)[-1])
  dat_fix$CHROM <- as.numeric(dat_fix$CHROM)
  ###簡化版snp matrix (imputation後修正過的genotype)
  cols <- colnames(abh_geno)[abh_geno[1, ] == chr]
  cols = na.omit(cols)
  true_geno <-abh_geno[,cols]
  true_geno <- true_geno[-1, ]
  ###
  geno_pos <- colnames(true_geno) 
  num_geno_pos <- as.numeric(sapply(strsplit(geno_pos,"_"),function(x){x[2]}) )
  #先解除lib vcfR，後面才可以使用seqinr的write.fasta
  detach("package:vcfR",unload = T)
  library(seqinr)
  ###
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
    assign(name[jj],toupper(data[[chr]][(position-one_side_length):(position+one_side_length)]))
    ##先幫每個個體建立斷點兩側序列
    temp <- toupper(data[[chr]][(position-one_side_length):(position+one_side_length)])
    ##針對每個點進行替換 (要替換成最接近true_geno位點的snp的基因型)
    for (kk in 1:length(snp_pos_within)) {
    #先確認snp的基因型以及其周圍tru_geno的基因型
    snp_vcf_type = dat_gt[jj+1,which(dat_fix$POS==snp_pos_within[kk] & dat_fix$CHROM == chr)]
    vcf_alt <- dat_fix$ALT[which(dat_fix$POS==snp_pos_within[kk] & dat_fix$CHROM == chr)]
    vcf_ref <- dat_fix$REF[which(dat_fix$POS==snp_pos_within[kk] & dat_fix$CHROM == chr)]
    true_geno_left <- true_geno[jj,findInterval(snp_pos_within[kk],num_geno_pos)]
    true_geno_right <- true_geno [jj,findInterval(snp_pos_within[kk] +1,num_geno_pos)]
    ####找出snp_within位點的major allele
    major_allele <- dat_gt[-1,which(dat_fix$POS==snp_pos_within[kk] & dat_fix$CHROM == chr)]
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
  #####建立一個fas_dat的list，包含所有個體的序列，使用該檔案進行計算
  write.fasta(fas_dat, as.string = F,names = name,file.out  = paste0(position-one_side_length,"_",position+one_side_length,".fasta"))
  cat("snp amount: ",snp_pos_within,"\n")
  cat("file:",paste0(position-one_side_length,"_",position+one_side_length,".fasta"),"created","\n")
}

#################running the function
fas_generate(chr = hotspot_interval$chr[1],position = hotspot_interval$start[1]+50000,one_side_length = 50000,vcf_name="gbs_r.chr",
             vcf_reffas_path= "C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data",ref_data = data)
#### local
setwd("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result")

for (i in 1:25) {
  fas_generate(chr = hotspot_interval$chr[i],
               position = round((hotspot_interval$start[i]+hotspot_interval$end[i])/2,0),
               one_side_length = round((hotspot_interval$start[i]+hotspot_interval$end[i])/2,0) - hotspot_interval$start[i],
               vcf_name="gbs_r.chr",
               vcf_reffas_path= "C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data",
               ref_data = data)
}

###overlap
setwd()
for (i in 178:198) {
  fas_generate(chr = hotspot_interval$chr[i],
               position = round((hotspot_interval$start[i]+hotspot_interval$end[i])/2,0),
               one_side_length = round((hotspot_interval$start[i]+hotspot_interval$end[i])/2,0) - hotspot_interval$start[i],
               vcf_name="gbs_r.chr",
               vcf_reffas_path= "C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data",
               ref_data = data)
}

for (i in 1:93) {
  fas_generate(chr = poisson_hotspot$chr[i],
               position = round((poisson_hotspot$start[i]+poisson_hotspot$end[i])/2,0),
               one_side_length = round((poisson_hotspot$start[i]+poisson_hotspot$end[i])/2,0) - poisson_hotspot$start[i],
               vcf_name="gbs_r.chr",
               vcf_reffas_path= "C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data",
               ref_data = data)
}

####產生random區間

ram_chr = c(sample(1:12,24,replace = T))
ram_start <- c()
for (i in 1:24) {
  repeat {
    position = as.numeric(sample(1000:length(data[[ram_chr[i]]]),1,replace = F))
    
    #如果position的位置沒有出現在熱點區間內，就停止抽取
    if (  length(which(apply(hotspot_interval %>% filter(chr == i),1,function(x){
      between(position,as.numeric(x[2]),as.numeric(x[3]))})))==0 ) {
      cat("random site:",position,"\n")
      ram_start <- c(ram_start,position)
      break
    } else {next}
  }
}

ramdom_interval <- data.frame(chr = ram_chr, start = ram_start, end = ram_start+50000)
for (i in 1:24) {
  fas_generate(chr = ramdom_interval$chr[i],
               position = round((ramdom_interval$start[i]+ramdom_interval$end[i])/2,0),
               one_side_length = round((ramdom_interval$start[i]+ramdom_interval$end[i])/2,0) - ramdom_interval$start[i],
               vcf_name="gbs_r.chr",
               vcf_reffas_path= "C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data",
               ref_data = data)
}

