####目標是創造overlap hotspot中CO event 上下游序列的fasta
library(seqinr)
overlap <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_overlap_hotspot.csv",header = T)
breakpoint <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_breakpoint.csv",header = T)
data <- read.fasta("C:/Users/lingo1st/Dropbox/碩論/np ir64/REF_genome.fa/OsNB1.0.fa")
abh <-  read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_abhgenotype.csv",header = T)
####找出marker位置

marker <- colnames(abh)
marker <- strsplit(marker,split = "_")
marker <- sapply(marker,function(x){
  x[2]
})
marker_chr <-unlist( as.vector(abh[1,]))
marker_df <- data.frame(marker = marker,chr = marker_chr)
marker_df <- marker_df[-1,]
######## 找出有位在overlap區間內的breakpoint位置
####
co_in_hotspot <- c()
chr <- c()
for (i in 1:nrow(breakpoint)) {
  df <- overlap %>% filter(chr == breakpoint$Chr[i])
  if(nrow(df) == 0){next}else{
    for (j in 1:nrow(df)) {
      if(between(breakpoint$Pos[i],df$start[j],df$end[j]) ){
        co_in_hotspot <- c(co_in_hotspot,breakpoint$Pos[i])
        chr <- c(chr,breakpoint$Chr[i])
      }
    }
  }
}

df_for_motif <- data.frame(chr = chr, pos = co_in_hotspot) 
df_for_motif <- unique(df_for_motif)

#######找出bkpt的左右兩個marker
i=1
start <- c()
end <- c()
mark_chr <- c()
for (i in 1:nrow(df_for_motif)) {
  temp_marker <- marker_df %>% filter(chr == df_for_motif$chr[i])
  temp_pos <- as.numeric(as.vector(temp_marker$marker))
  left_marker <- findInterval(df_for_motif$pos[i],temp_pos)
  start<- c(start,temp_pos[left_marker])
  end <- c(end,temp_pos[left_marker+1])
  mark_chr <- c(mark_chr,df_for_motif$chr[i])
}
df_for_motif$start <- start
df_for_motif$end <- end


##############
#製造fasta檔
##############
###依照兩側marker區間建立fasta(但長度會不同)
i=1
fas_dat <- list()
for (i in 1:nrow(df_for_motif)) {
  temp_fas<- data[[df_for_motif$chr[i]]]
  fas <- temp_fas[df_for_motif$start[i]:df_for_motif$end[i]]
  fas_dat[[i]] <- fas
  names(fas_dat[[i]]) <- i
}
write.fasta(fas_dat,names =names(fas_dat) ,file.out = paste0("chr",df_for_motif$chr[i],"_",df_for_motif$start[i],"_",df_for_motif$end[i],".fasta"))

###依照中心點向外擴散(1kb)
fas_dat_2 <- list()
for (i in 1:nrow(df_for_motif)) {
 temp_fas<- data[[df_for_motif$chr[i]]]
 fas_dat_2[[i]] <- temp_fas[(round(df_for_motif$pos[i],0)-500):(round(df_for_motif$pos[i],0)+500)]
 names(fas_dat_2[[i]]) <- i
}
names(fas_dat_2) <- as.character(df_for_motif$pos)
write.fasta(fas_dat_2,names = names(fas_dat_2),file.out = "motif_overlap_hotspot.fasta")
