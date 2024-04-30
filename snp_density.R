library(vcfR)
###########
largest_density_size <- c()
largest_density_pos <- c()
for (i in 4:12) {
  data <- read.vcfR(paste0("C:/Users/lingo1st/Desktop/NOISYmputer/NOISYmputer_data/final_sample_plus_parent.chr",i,".vcf"))
  data_fix <- data@fix
  data_fix <- as.data.frame(data_fix)
  data_fix$POS <- as.numeric(data_fix$POS)
  b <- c()
  for (j in seq(1,max(data_fix$POS),100000)) {
    sub_data <- data_fix[data_fix$POS<= (j+99999) & data_fix$POS >= j,]
    b <- append(b,as.numeric(nrow(sub_data)))
  }
  largest_density_pos[i] <- c(100000*which.max(b))
  largest_density_size[i] <- max(b)
}

####每條染色體密度最高的區間(0.1Mb)，以及裡面的marker數量
snp_den <- data.frame(chr = as.character(c(1:12)), size = largest_density_size, pos = largest_density_pos)

####
chr_len <- c(43260640,35954074,36189985,35489479,29733216,30731386,29643843,28434680,22692709,22683701,28357783,27561960)
snp_amount <- c(275319,217875,225028,189639,145899,190352,188189,175811,142283,167227,200622,175713)
round(sum(chr_len)/sum(snp_amount,1))

data <- read.vcfR(paste0("C:/Users/lingo1st/Desktop/NOISYmputer/NOISYmputer_data/final_sample_plus_parent.chr",i,".vcf"))
data_fix <- data@fix
data_fix <- as.data.frame(data_fix)
data_fix$POS <- as.numeric(data_fix$POS)
distance <- c()
for (i in 1:nrow(data_fix)) {
  distance[i] <- data_fix$POS[2] - data_fix$POS[1]
}
which.max(na.omit(distance))
data_fix$POS[which.max(na.omit(distance))]
