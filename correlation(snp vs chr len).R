library(tidyverse)
library(dplyr)
chr_len <- c(43260640,35954074,36189985,35489479,29733216,30731386,29643843,28434680,22692709,22683701,28357783,27561960)/1000000
snp <- c(6160,4695,3630,3948,3091,1532,3479,3201,2543,3027,3749,2282)
row.names(df) <- c("Chr1","Chr2","Chr3","Chr4","Chr5","Chr6","Chr7","Chr8","Chr9","Chr10","Chr11",'Chr12')
df <- data.frame(chr_len = chr_len, amount = snp)
p1 <- ggplot(data = df,aes(x = chr_len,y = amount))
p1 + geom_point()+
  geom_smooth(method = 'lm', formula =y ~ x)+
  geom_label(
    label=rownames(df), 
    nudge_x = 0.25, nudge_y = 2
    
  )+xlab("Chromosome length (Mb)")+ ylab("SNP amount")

m1 <- lm(amount~chr_len,data = df)
summary(m1)
cor.test(df$chr_len,df$amount)
sum(cent_max - cent_min)


#### co num vs chr_len
breakpoint <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_breakpoint.csv",header = T)
co_num <- breakpoint %>% group_by(Chr) %>%
  summarise(num = n())
cor.test(co_num$num,chr_len)



######## CO event vs chromosome length
breakpoint <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_breakpoint.csv",header = T)
co <- breakpoint %>% 
  group_by(Chr) %>%
  summarise(total = n())
CO_length_df <- data.frame(len = chr_len,co = co$total[-13])
p2 <- ggplot(data = CO_length_df,aes(x = len/1000000,y=co))
p2 + geom_point() + geom_smooth(method = lm,formula = y~x) + ylab("Amount of CO event")+
  xlab("Chromosome length (Mb)")
cor(CO_length_df$len,CO_length_df$co)
m1 <- lm(total~Chr,data = co)
summary(m1)