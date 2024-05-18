####### 使用 gbs_breakpoint.csv
data <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gbs_dataset_result/csv file/gbs_breakpoint.csv")

chr_len <- c(43260640,35954074,36189985,35489479,29733216,30731386,29643843,28434680,22692709,22683701,28357783,27561960)

####

hotspot <- data.frame(matrix(ncol = 4,nrow = 0))
colnames(hotspot) <- c("chr","event_number","start","end")

windowsize = 100000
p_value <- c()
###跑第一次只需要取co_num的部分，取出來做zero-inflate建模，之後使用建模後的lamda及pi進行pvalue計算，之後再BH 校正
for (chromosome in 1:12) {
  tmp1 <- data[data$Chr== chromosome,]
  #lamda = length(tmp1$pos)/(chr_len[chromosome]/windowsize)
  #lamda = 6980/(sum(chr_len)/windowsize)
  lamda = nrow(tmp1) / (chr_len[chromosome]/windowsize)
  #lambda <- exp(coef(zip_model)[1])
  interval <- seq(1,chr_len[chromosome],by = windowsize)
  ###
  for (i in 1:(length(interval)-1)) {
    temp <- data.frame(chr = chromosome,
                       co_num = length(which(between(tmp1$Pos,interval[i],interval[i+1]))),
                       p_value = (1-pi)* dpois(length(which(between(tmp1$Pos,interval[i],interval[i+1]))),lamda),
                       
                         start = interval[i],end =interval[i+1] )
    hotspot <- rbind(hotspot, temp)
    
  }
}
recomb_counts <- hotspot$co_num
###########################zero-inflated poisson
install.packages("psc1")
library(pscl)
# 拟合ZIP模型
zip_model <- zeroinfl(recomb_counts ~ 1 | 1, dist = "poisson")

# 查看模型摘要
summary(zip_model)

# 获取拟合模型的参数
lambda <- exp(coef(zip_model)[1])  # 泊松分布的均值参数
pi <- plogis(coef(zip_model)[2])    # 零膨胀部分的概率
(1 - pi) * dpois(0, lambda)

###########################p value 調整 BH method
p_value <- hotspot$p_value
adj_p <- p.adjust(p_value,method="BH")
hotspot$adj_p_value <- adj_p
p <- (hotspot[hotspot$adj_p_value < 0.01,])
write.csv(p,file = "adjusted_poisson_hotspot.csv")


###比對
wgs_co <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/imputation_plus_hotspot/all_co_event.csv",header = T)
which(wgs_co$start %in% hotspot$start)
for (i in 1:12) {
  a <- wgs_co[wgs_co$chr == i,]
  b <- p[p$chr == i,]
  print(a[which(a$start %in% b$start),])
}

####卡方檢定
# 计算观测频数
obs_freq <- table(hotspot$co_num)

# 计算泊松分布的期望频数
lambda_est <- mean(hotspot$co_num)
exp_freq <- dpois(as.numeric(names(obs_freq)), lambda = lambda_est) * length(hotspot$co_num)

# 进行卡方检验
chisq_test <- chisq.test(obs_freq, p = exp_freq, rescale.p = TRUE)
###goodness of fit
fit <- goodfit(hotspot$co_num, type = "poisson")
summary(fit)
recomb_counts <- hotspot$co_num[order(hotspot$co_num)]
lambda_poisson <- mean(recomb_counts)
x <- 0:max(recomb_counts)
pois_fit <- dpois(x, 1)

# 查看检验结果
print(chisq_test)
install.packages("MASS")
install.packages("vcd")
library(vcd)
library(MASS)
hist(recomb_counts, breaks = max(recomb_counts) - min(recomb_counts), probability = TRUE, 
     main = "Histogram of Recombination Counts", xlab = "Recombination Count")
fit_nbinom <- fitdistr(recomb_counts, "negative binomial")
size_nbinom <- fit_nbinom$estimate["size"]
mu_nbinom <- fit_nbinom$estimate["mu"]
nbinom_fit <- dnbinom(x, size = size_nbinom, mu = mu_nbinom)
lines(x, nbinom_fit, type = "b", col = "blue", lwd = 2)

obs_freq <- table(recomb_counts)
# 计算负二项分布的期望频数
exp_freq_nbinom <- dnbinom(as.numeric(names(obs_freq)), size = size_nbinom, mu = mu_nbinom) * length(recomb_counts)
# 进行卡方检验
chisq_test_nbinom <- chisq.test(obs_freq, p = exp_freq_nbinom, rescale.p = TRUE)
print(chisq_test_nbinom)
rpois()



hist(recomb_counts, breaks = max(recomb_counts) - min(recomb_counts), probability = TRUE, 
     main = "Histogram of Recombination Counts", xlab = "Recombination Count")

# 泊松分布拟合

lambda_poisson <- mean(recomb_counts)
x <- 0:max(recomb_counts)
pois_fit <- dpois(x, lambda_poisson)
lines(x, pois_fit, type = "b", col = "red", lwd = 2)
