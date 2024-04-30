install.packages("LinkageMapView")
library("LinkageMapView")
col <- c('set','map','mkr','phys','gen','vld')
total <- data.frame(matrix(nrow = 0,ncol = length(col)))
colnames(total) <- col
for (i in 1:12) {
  data <- read.csv(paste0('C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation result/','chr',i,'/genetic_distance_chr',i,'.csv'))
  #選取基數行
  data <- data[seq(1,nrow(data),by = 2),]
  assign(paste0("gen_distance_",i),data,envir = .GlobalEnv)
  total <- rbind(total,data)
}

genetic_map <- data.frame(group = total$map,position = total$gen,locus = total$mkr)
colnames(genetic_map) <- c('group','locus','position')
genetic_map <- na.omit(genetic_map)
###
outfile = file.path('C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation result', "genetic_map_edit.pdf")
lmv.linkage.plot(genetic_map,outfile)



i=1
