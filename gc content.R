install.packages("vcfR")
install.packages("seqinr")
library(vcfR)
library(seqinr)
library(dplyr)
library(tidyverse)
library(tidyr)
install.packages("LncFinder")
library(LncFinder)
###########資料匯入
local_hotspot <- c(20522947, 20537617,20545769,
                   20556928,20672614,27359627,
                   27362882,27425162,27435631,
                   27526405,2107364,2888125,
                   2925595,19662508,19680032,
                   19849482)
local_chr <-c(5,5,5,5,5,5,5,5,5,5,7,12,12,12,12,12)

poisson_hotspot <- c(500001, 1200001, 4700001,9000001,
                     26400001, 27700001,28900001,32100001,
                     35100001, 40400001, 42200001,300001,
                     900001,5100001,8100001,8500001,
                     21300001,28100001,29400001,30700001,
                     34100001,35200001,3800001,14400001,
                     17100001,27700001, 28500001,32200001,
                     32600001, 35900001,19900001,26500001,29800001,35300001,
                     20500001,27300001,27400001,29500001,6600001,27800001,
                     30400001,900001,1200001, 1500001,2100001,2200001,2700001, 5700001,
                     24700001,27900001,700001,800001, 2200001,20600001,24600001,25300001,
                     19100001,11900001, 16700001,100001,2700001, 5500001,21200001, 23200001,
                     25100001,2800001, 19600001)
poisson_chr <- c(rep(1:12,times = c(11,11,8,4,4,3,9,6,1,2,6,2)))

###
###data 前處理
for (q in 1:12) {
  
  dat <- read.csv(paste0("C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation_hotspot_result/chr",q,"/imputation/final_sample_plus_parent.chr",q,"_imputed_chunks.txt"),sep = " ")
  p <- strsplit(dat$genotypes,"")
  for (i in 1:length(p)) {
    for (j in 1:length(p[[1]])) {
      dat[i,j+2] <- p[[i]][j]
    }
  }
  row_name <- dat$X.
  rownames(dat) <- row_name
  dat <- dat[,-1:-2]  %>% as.data.frame() %>% t()
  #write.csv(dat,file = paste0("C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation result/chr",q,"/",q,"_snp matrix.csv"))
  assign(paste0("snp_",q),as.data.frame(dat) ,envir = .GlobalEnv)
}
####
data <-  read.fasta("C:/Users/lingo1st/Dropbox/碩論/np ir64/REF_genome.fa/OsNB1.0.fa")
vcf_reffas_path="C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data"
vcf_name = "final_sample_plus_parent.chr"
ref_data <-data
position = 20522947
chr = 1
one_side_length = 50000
window_size = 10000

###########
gc_content <- c()
gc_content_p <- c()
gc_content_l <-c()
nuc_diversity_v3 <- function(fas_data_path=NA,chr = NA,vcf_reffas_path=NA,vcf_name=NA,ref_data,position=NA,one_side_length,window_size,sliding = 0,output_path,random_pos = F,repeat_time = 1){
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
  if(length(snp_pos_within)==0){next}
  for (j in 1:length(name)) {
    ##先幫每個個體建立斷點兩側序列
    assign(name[j],toupper(data[[chr]][(position-one_side_length):(position+one_side_length)]))
    
    ##針對每個點進行替換
    temp <- get(name[j])
    for (k in 1:length(snp_pos_within)) {
      
      ##找出ref,alt allele
      change_pos <- as.numeric(which(dat_fix$CHROM == chr &  dat_fix$POS == snp_pos_within[k]))
      ref <- dat_fix$REF[change_pos]
      alt <- dat_fix$ALT[change_pos]
      ###確認該點於該個體的genotype並且進行替換
      #將原始序列NA值替換成N
      if (is.na( temp[(snp_pos_within[k]-(position-one_side_length))+1] )){temp[(snp_pos_within[k]-(position-one_side_length))+1] = "N" }
      #將dat_gt在置換位點上的NA改成N
      if (is.na(dat_gt[,colnames(dat_gt)==name[j]][change_pos])){dat_gt[,colnames(dat_gt)==name[j]][change_pos] = "N" }
      #如果序列上置換位點不屬於親本的任一ALLELE則發出警示
      if (temp[(snp_pos_within[k]-(position-one_side_length))+1] != ref & temp[(snp_pos_within[k]-(position-one_side_length))+1] != alt &temp[(snp_pos_within[k]-(position-one_side_length))+1] != "N" ){
        cat("warning,sth is wrong")
        break
      }
      if(dat_gt[,colnames(dat_gt)==name[j]][change_pos] =="1/1"){
        temp[snp_pos_within[k]-(position-one_side_length)+1] <- ref
      }else if (dat_gt[,colnames(dat_gt)==name[j]][change_pos] =="0/0"){
        temp[snp_pos_within[k]-(position-one_side_length)+1] <- alt
      }
      
    }
    fas_dat[[j]] <- temp
  }
  
  names(fas_dat) <- name
  #####建立一個fas_dat的list，包含所有個體的序列，使用該檔案進行計算(GC CONTENT)
  #write.fasta(fas_dat,names = name,as.string = F,file.out = paste0("C:/Users/lingo1st/Dropbox/林冠瑜/gc content","/",position,"_hotspot.fasta") )
  
  P <- sapply(fas_dat,function(x){
    GC(x,forceToLower = T)
  })
  gc_content_overlap <<- c( gc_content_overlap,mean(P) ) 
  
}

##############################################################################################################################
#####
for (i in 1:3) {
  nuc_diversity_v3(
    data <-  read.fasta("C:/Users/lingo1st/Dropbox/碩論/np ir64/REF_genome.fa/OsNB1.0.fa"),
    vcf_reffas_path="C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data",
    vcf_name = "final_sample_plus_parent.chr",
    ref_data <-data,
    position = local_hotspot[i],
    chr = local_chr[i],
    one_side_length = 5000,
    window_size = 10000
  )
}
###whole genome
gc_content_whole <- c()
for (i in 1:12) {
  g <- GC(data[[i]],forceToLower = T,exact = T)
  gc_content_whole <- c(gc_content_whole,g)
}

####計算重疊區域
gc_content_overlap <- c()
over_lap <- c(2055001,27354814,27488203,2128682,2869062,19656254)
one_side <- c(5000,45000,11798,28682,30938,43747)
over_lap_chr <- c(5,5,5,7,12,12)

for (i in 1:6) {
  nuc_diversity_v3(
    data <-  read.fasta("C:/Users/lingo1st/Dropbox/碩論/np ir64/REF_genome.fa/OsNB1.0.fa"),
    vcf_reffas_path="C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data",
    vcf_name = "final_sample_plus_parent.chr",
    ref_data <-data,
    position = over_lap[i],
    chr = over_lap_chr[i],
    one_side_length = one_side[i],
    window_size = 10000
  )
}
############################################################################################
#plot

plot_df <- data.frame(GC = c(gc_content,gc_content_p,gc_content_whole,gc_content_overlap),Type = rep(c("Local","Poisson","Whole-genome","Overlapped hotspot"),times = c(16,68,12,6)) )
plot_df <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gc content/gc content_all interval.csv",header = F)
colnames(plot_df) <- c("GC","Type")
p1 <- ggplot(data = plot_df, aes(x = Type,y = GC, fill = Type))
p1 + geom_boxplot()+stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") + scale_x_discrete(limits = c("L","P","Overlap","Non-hotspot"))+ylab("GC content")+ 
  geom_signif(comparisons = list(c("Overlap", "L"),
                                 c("Overlap","Non-hotspot")),
              map_signif_level = TRUE,
              y_position = c(0.475,0.479))
model <- aov(GC~Type,data = plot_df)
summary(model)
