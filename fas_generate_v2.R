####
###snp matrix
abh_geno <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_abhgenotype.csv",header = T)

############叫出bkpt
breakpoint <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_breakpoint.csv",header = T)

data <- read.fasta("C:/Users/lingo1st/Dropbox/碩論/np ir64/REF_genome.fa/OsNB1.0.fa")

#########建立熱點區間
local_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/local_recomb/local_hotspot_spar_0.1.csv",header = T)
poisson_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_adjusted_poisson_hotspot.csv",header = T)
combine_hotspot <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_hotspot_overlap.csv",header = T)
names(local_hotspot)[names(local_hotspot) == "map"] <- "chr"
local_hotspot$chr <- gsub("chr","",local_hotspot$chr)

###綜合成hotspot_interval
hotspot_interval <- rbind(local_hotspot[,c(3,6,7)],poisson_hotspot[,c(2,5,6)])
hotspot_interval$type <- rep(c("local","poisson"),times = c(145,152))
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
  write.fasta(fas_dat,as.string = F,names = name,file.out = paste0(position-one_side_length,"_",position+one_side_length,".fasta"))
  cat("snp amount: ",snp_pos_within,"\n")
  cat("file:",paste0(position-one_side_length,"_",position+one_side_length,".fasta"),"created","\n")
  }
fas_generate(chr = hotspot_interval$chr[1],position = hotspot_interval$start[1]+50000,one_side_length = 50000,vcf_name="gbs_r.chr",
             vcf_reffas_path= "C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data",ref_data = data)
for (i in 219:297) {
  fas_generate(chr = hotspot_interval$chr[i],position = hotspot_interval$start[i]+50000,one_side_length = 50000,vcf_name="gbs_r.chr",
               vcf_reffas_path= "C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data",ref_data = data)
}


