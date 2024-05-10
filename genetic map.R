install.packages("LinkageMapView")
library("LinkageMapView")
###轉換格式
##column one : group (ex:chr1)
## column two: position (ex:445566)
## column three: locus (ex:s1_445566) (snp代號)
genetic_map <- data.frame(matrix(nrow = 0,ncol = 3))
for (i in 1:12) {
  data <- read.csv(paste0('C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/','chr',i,'/genetic_map_chr',i,'.csv'),header = F)
  uni <- unique(data$V2)
  data_unique <- data[data$V2 %in% uni & !duplicated(data$V2,fromLast = T), ]

  chr = rep(paste0("chr",i),nrow(data_unique))
  temp_df <- data.frame(group = chr, position = data_unique$V2, locus = data_unique$V1)
  genetic_map <- rbind(genetic_map, temp_df)
}
genetic_map <- na.omit(genetic_map)
###
outfile = file.path('C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result', "gbs_genetic_map.pdf")
lmv.linkage.plot(genetic_map,outfile)


#####
#找出每條染色體最大的遺傳距離
genetic_map %>% group_by(group) %>% summarise(max = max(position))









col <- c('set','map','mkr','phys','gen','vld')
total <- data.frame(matrix(nrow = 0,ncol = length(col)))
colnames(total) <- col
"C:\Users\lingo1st\Dropbox\林冠瑜\gbs_dataset_result\chr1\genetic map_chr1.csv"

genetic_map <- data.frame(group = total$map,position = total$gen,locus = total$mkr)
colnames(genetic_map) <- c('group','locus','position')




i=1
