####擷取部分vcf
###可直接針對gt,fix的部分以matrix形式進行操作
library(vcfR)
data <- read.vcfR("C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data/gbs_T.chr1.vcf",verbose = F)
data_fix <- data@fix
data_gt <- data@gt

####找出np,ir64相同的snp位置並去除
rm <- which(data_gt[,"Nipponbare"] == data_gt[,"IR64"])
new_data_fix <- data_fix[-rm,]
new_data_gt <- data_gt[-rm,]
data@fix <- new_data_fix
data@gt <- new_data_gt
data_fix <- data@fix
data_gt <- data@gt

###將Np中0/0的部分換成1/1，同時ref與alt也要改變
change_index <- which(data@gt[,"Nipponbare"] =="0/0")

###先觀察員本資料
data@gt[head(change_index),]
data@fix[head(change_index),]
###進行替換
for (i in 1:length(change_index)) {
  one_one <- which(data@gt[change_index[i],] == "1/1")
  zero_zero <-which(data@gt[change_index[i],] == "0/0")
  data@gt[change_index[i], one_one] <- "0/0"
  data@gt[change_index[i],zero_zero] <- "1/1"
  ###換fix的ref,alt allele 
  REF <- data@fix[change_index[i],"REF"]
  data@fix[change_index[i],"REF"] <- data@fix[change_index[i],"ALT"]
  data@fix[change_index[i],"ALT"] <- REF
}

write.vcf(data,file = "C:/Users/lingo1st/OneDrive/桌面/gbs_r.chr1.vcf.gz")

#####改全部染色體
for (j in 2:12) {
  data <- read.vcfR(paste0("C:/Users/lingo1st/OneDrive/桌面/NOISYmputer/NOISYmputer_data/gbs_T.chr",j,".vcf"),verbose = F)
  ###去除沒有多型性的snp
  rm <- which(data@gt[,"Nipponbare"] == data@gt[,"IR64"])
  data@gt <- data@gt[-rm,]
  data@fix <- data@fix[-rm,]
  ######
  ###將Np中0/0的部分換成1/1，同時ref與alt也要改變
  change_index <- which(data@gt[,"Nipponbare"] =="0/0")
  for (i in 1:length(change_index)) {
    one_one <- which(data@gt[change_index[i],] == "1/1")
    zero_zero <-which(data@gt[change_index[i],] == "0/0")
    data@gt[change_index[i], one_one] <- "0/0"
    data@gt[change_index[i],zero_zero] <- "1/1"
    ###換fix的ref,alt allele 
    REF <- data@fix[change_index[i],"REF"]
    data@fix[change_index[i],"REF"] <- data@fix[change_index[i],"ALT"]
    data@fix[change_index[i],"ALT"] <- REF
  }
  write.vcf(data,file = paste0("C:/Users/lingo1st/OneDrive/桌面/gbs_r.chr",j,".vcf.gz"))
  
}

write.csv(colnames(data@gt),file = "name.txt",row.names = F,quote = F)
