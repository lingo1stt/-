install.packages("LinkageMapView")
library("LinkageMapView")
install.packages("ggpubr")
library(ggpubr)
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

write.csv(genetic_map,file = "gbs_whole_genetic_distance.csv")
#####
#找出每條染色體最大的遺傳距離
genetic_map %>% group_by(group) %>% summarise(max = max(position))


####遺傳距離vs 物理距離 (使用genetic map dataframe)
genetic_map$phys <- as.numeric(sapply(strsplit(genetic_map$locus,"_"),function(x){
  x[2]
}) 
)
for (i in 1:12) {
  p <- ggplot(data = genetic_map[genetic_map$group == paste0("chr",i),],aes(x = phys/1000000,y = position)) + geom_point() + 
    geom_smooth(method = lm, formula = y~x)  + ylab("Genetic distance (cM)") + labs(x = NULL) + labs(y = NULL)
  assign(paste0("p",i), p, envir = .GlobalEnv)
  ###計算一次式的r squared
  model <- lm(position~phys, data = genetic_map[genetic_map$group == paste0("chr",i),])
  temp_summary <- summary(model)
  r2<- round(temp_summary[["adj.r.squared"]], 3)
  cat("chr",i,": ",r2,"\n")
}
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,nrow = 2,ncol = 6, labels = c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6"
                                                                               ,"Chr7","Chr8","Chr9","Chr10","Ch11","Chr12"),
          vjust = 3,hjust = -1,font.label = list(size = 8))



