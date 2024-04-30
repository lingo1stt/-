library(qtl)
library(ABHgenotypeR)
data <-  read.csv(paste0("C:/Users/lingo1st/Dropbox/林冠瑜/noisymputer result_2/imputation result/",chr,"/final_sample_plus_parent.",chr,"_imputed_chunks.txt")
                  ,sep = " ")
####主要迴圈
for (chr in 1:12) {
  data <-  read.csv(paste0("C:/Users/lingo1st/Dropbox/林冠瑜/noisymputer result_2/imputation result/chr",chr,"/final_sample_plus_parent.chr",chr,"_imputed_chunks.txt")
                    ,sep = " ")
  ###data 前處理
  p <- strsplit(data$genotypes,"")
  for (i in 1:length(p)) {
    for (j in 1:length(p[[1]])) {
      data[i,j+2] <- p[[i]][j]
    }
  }
  row_name <- data$X.
  rownames(data) <- row_name
  data <- data[,-1:-2]  %>% as.data.frame()
  test_data <- t(data)%>% as.data.frame()
  #整理成必要格式
  rown <- rownames(data)
  test <-  data.frame(matrix(nrow = 0,ncol = length(rown)))
  test[1,] <- substr(rown,2,3)
  colnames(test) <- rown
  test2 <- rbind(test,test_data)
  assign(paste0("geno_",chr),test2,envir = .GlobalEnv)
}
total <- cbind(geno_1,geno_2,geno_3,geno_4,geno_5,geno_6,geno_7,geno_8,geno_9,geno_10,
               geno_11,geno_12)
write.csv(total,"C:/Users/lingo1st/Desktop/data.csv",quote = TRUE)
##後續又再改成data_for_binmap.csv並移至imputation result資料夾中
genotype <- readABHgenotypes("C:/Users/lingo1st/Desktop/data.csv", nameA = "NP", nameB = "IR64", readPos = TRUE)
plotGenos(genotype)