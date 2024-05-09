library(qtl)
install.packages("ABHgenotypeR")
library(ABHgenotypeR)
ind_name <- read.table("name.txt")
colnames(ind_name) <- "id"
temp <- data.frame(id = NA)
ind_name <- rbind(temp,ind_name)
snp_name <- c()
###以tassel中abh genotype產出之csv檔為格式範例
#column為snp,row為個體,rownames跟colnames都需要設定為個體與snp之名稱
### column one為id，第一個值為NA，後面接個體名稱
###row one為染色體數，同樣第一個值為NA，後面接染色體條數

####主要迴圈
for (chr in 1:12) {
  data <-  read.csv(paste0("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/chr",chr,"/gbs_r.chr",chr,"_imputed_chunks.txt")
                    ,sep = "\t")
  ###
  data <- data[-(1:7),]
  
  ##把snp名稱與基因型分開 (建立一個名為p的genotype 矩陣)
  temp <- strsplit(data,split = " ")
  snp_name <- sapply(temp,function(x){
     x[1]
  }) 
  
 p <- sapply(temp,function(x){
    strsplit( x[2],split = "")
  })
 p <- unlist(p)
 p <- matrix(p,nrow = 296 ,ncol = length(p)/296,byrow = F)
 colnames(p) <- snp_name
 ###多加染色體row
 chromosome <- matrix(rep(chr,length(snp_name)),nrow = 1)
 p <- rbind(chromosome, p)
 ###
 ind_name <- cbind(ind_name,p)
}
  


  
  

write.csv(ind_name,"gbs_abhgenotype.csv",quote = TRUE,row.names = F)
##後續又再改成data_for_binmap.csv並移至imputation result資料夾中
genotype <- readABHgenotypes("gbs_abhgenotype.csv", nameA = "Nipponbare", nameB = "IR64", readPos = TRUE)
plotGenos(genotype)
