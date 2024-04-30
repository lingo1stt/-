library(circlize)

###data資料集 (local recombination rate)
data <- read.table("C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation result/all_local_recomb.txt",header = T)
data$map <- gsub("chromosome ","chr",data$map)
cytoband_df <- data.frame(chrom = data$map,chromStart = data$phys,chromoEnd = data$phys,value1 = data$spline)

#### 重組事件資料集
#把"$` total without 3&24`"以後的row進行處理
recomb_pos <- c()
chr <- c()
for (i in 1 :12 ) {
  data_recomb_event <- read.csv(paste0("C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/imputation_hotspot_result/chr",i,"/imputation/breakpoint_chr",i,".csv"),
                                header = FALSE, stringsAsFactors = FALSE)
  #擷取"total without 3&24後面的row(為所有個體之總和)"
  for (j in (which( data_recomb_event[,1] == "$` total without 3&24`")+1) : nrow(data_recomb_event)) {
    #將dataframe轉為vector並且將na去除
    t <- na.omit(as.numeric(data_recomb_event[j,]))
    chr <- c(chr,rep(paste0("chr",i),length(t)))
    recomb_pos <- c(recomb_pos,t)
  }
}
y <-  rnorm(n = 610, mean = 20, sd = 3)
recomb_df <- data.frame(y = y,x = recomb_pos,chr = chr)
#############建立熱點區間 (hotspot interval)
local_hotspot <- c(20522947, 20537617,20545769,
                   20556928,20672614,27359627,
                   27362882,27425162,27435631,
                   27526405,2107364,2888125,
                   2925595,19662508,19680032,
                   19849482)
local_hotspot_s <- c()
local_chr <-c(5,5,5,5,5,5,5,5,5,5,7,12,12,12,12,12)
poisson_hotspot <- c(500001, 1200001, 4700001,9000001,
                     26400001, 27700001,28900001,32100001,
                     35100001, 40400001, 42200001,300001,
                     900001,5100001,8100001,8500001,
                     21300001,28100001,29400001,30700001,
                     34100001,35200001,3800001,14400001,17100001,27700001,28500001,
                     32200001,32600001,35900001,19900001,26500001,29800001,35300001,20500001,27300001,27400001,29500001,6600001,27800001,
                     30400001,900001,1200001, 1500001,2100001,2200001,2700001, 5700001,
                     24700001,27900001,700001,800001, 2200001,20600001,24600001,25300001,
                     19100001,11900001, 16700001,100001,2700001, 5500001,21200001, 23200001,
                     25100001,2800001, 19600001)
poisson_chr <- c(rep(1:12,times = c(11,11,8,4,4,3,9,6,1,2,6,2)))
hotspot_interval$chr <- paste0("chr",hotspot_interval$chr)
str(hotspot_interval)
##centromere region
cent_min <- c(15.445668,12.625206,17.883934,7.882,11.15,13.201,9.104,11.98,0.993,
              7.62,11.34,11.06)
cent_min <- cent_min * 1000000
cent_max <- c(18.05,15.48,20.51,10.06,13.54,17.84,12.71,14.41,3.93,8.69,13.5,
              12.2)
cent_max <- cent_max * 1000000 

################作圖
circos.clear()
circos.par("start.degree" = 90,gap.degree = c(1,1,1,1,1,1,1,1,1,1,1,16))
circos.initializeWithIdeogram(cytoband_df,plotType = c("axis", "labels"),axis.labels.cex = 0.3)
###local recombination rate
circos.track(sectors =cytoband_df$chrom , x = cytoband_df$chromStart, y = cytoband_df$value1,
             panel.fun = function(x, y) {
               circos.lines(x, y,area = T,border = F,col = "ivory4")
             })
###local hotsot
for (i in 1:16) {
  circos.rect(hotspot_interval$start[i], -5 ,hotspot_interval$end[i], 40,col = "red",lwd = 2,sector.index = hotspot_interval$chr[i],border = NA)
}
### 加入centromere 區域
for(i in 1:12){
  circos.rect(cent_min[i], -5 ,cent_max[i], 39,col = "azure3",lwd = 2,sector.index = paste0("chr",i),border = NA)
  
}

###no.2圈: 重組事件分布
##先創造
circos.track(sectors = recomb_df$chr , x = recomb_df$x, y = recomb_df$y,ylim = c(0,40),
             panel.fun = function(x, y) {
               circos.points(x, y,col = "lightskyblue4",cex = 0.1,pch = 21)
             })

##poisson hotspot
for (i in 17:83) {
  circos.rect(hotspot_interval$start[i], 0 ,hotspot_interval$end[i], 40,col = "red",lwd = 3,sector.index = hotspot_interval$chr[i],border = NA)
}
##

circos.genomicDensity(data = den_df,count_by = "number",col = "thistle")
circos.yaxis(side = "left", sector.index = "chr1",labels.cex = 0.3,track.index = 2,labels.niceFacing = T)
circos.yaxis(side = "left", sector.index = "chr1",labels.cex = 0.15,track.index = 4,labels.niceFacing = T)

##########加入density plot
library(vcfR)
den <- read.vcfR(file = "C:/Users/lingo1st/Dropbox/林冠瑜/final.vcf/final_sample_plus_parent.vcf")
den_g <- as.data.frame(getFIX(den)) 
str(den_g)
den_g$POS <- as.numeric(den_g$POS)
den_g$CHROM<- as.numeric(den_g$CHROM)
den_g$CHROM <- paste0("chr",den_g$CHROM)
den_df <- data.frame(chrom = den_g$CHROM, chromStart = den_g$POS,chromEnd = den_g$POS)
den_df$chromEnd <- den_df$chromEnd +1
