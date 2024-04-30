library(dplyr)
library(tidyverse)
data <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/noisymputer result/parent1_proportion.csv")
data <- data[1:26,]
rownames(data) <- data[,1]
plot_df <- data %>%
  pivot_longer(cols = -X, names_to = 'Identity', values_to = 'Proportion')
###TEST###
plot_df2 <- plot_df %>% filter(Identity == "X1b" | Identity == "X1f")
plot_df2$identity[plot_df2$identity == "X1b"] <- "Before imputation"
plot_df2$identity[plot_df2$identity == "X1f"] <- "Imputation_after"
p1 <- ggplot(data = plot_df2,aes(x = group , y =Proportion  ,fill = identity ))
p1 + geom_bar(position = "dodge", stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text(aes(label = Proportion) ,position = position_dodge(width = 1),size = 2.5,vjust = -1)

### LOOP
for (i in 1:12) {
  print(i)
  plot_df2 <-  plot_df %>% filter(Identity == paste0("X",i,"b") | Identity == paste0("X",i,"f"))
  plot_df2$Identity[plot_df2$Identity == paste0("X",i,"b")] <- "Before imputation"
  plot_df2$Identity[plot_df2$Identity == paste0("X",i,"f")] <- "Imputation_after"
  colnames(plot_df2) <- c("Individual","Identity","Proportion")
  p1 <- ggplot(data = plot_df2,aes(x =Individual , y =Proportion  ,fill = Identity ))+
    geom_bar(position = "dodge", stat = "identity")+
    ggtitle(paste0("Parent1(NP) proportion before/after imputation(chr",i,")"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_text(aes(label = Proportion) ,position = position_dodge(width = 2),size = 2,vjust = -1)
  pdf(file = paste0("C:/Users/lingo1st/Dropbox/林冠瑜/noisymputer result/parent1_proportion_chr",i,".pdf"))
  print(p1)
  dev.off()
}
i=1
