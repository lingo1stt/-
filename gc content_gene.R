library("seqinr")
install.packages("ggsignif")                     
library("ggsignif")   
library(tidyverse)
###data
whole_gene_data <- read.fasta("C:/Users/lingo1st/OneDrive/桌面/IRGSP-1.0_gene_2024-01-11.fasta/IRGSP-1.0_gene_2024-01-11.fasta",
                              forceDNAtolower = T)
all_genes_p <-read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gc content/all_gene_poisson.csv")
all_genes_p <- all_genes_p$x
all_genes_l <- read.csv("C:/Users/lingo1st/Dropbox/林冠瑜/gc content/all_gene_local.csv" )
all_genes_l <- all_genes_l$x
##
head(names(whole_gene_data))
p <- sub(pattern = "g",replacement = "t",all_genes)
grep(p[1],names(whole_gene_data))
names(whole_gene_data)[242]
GC3(whole_gene_data[[names(whole_gene_data)[242]]])

###先將各熱點區間內基因在whole_gene_data中的位置記錄下來

##local (0.490)
t_local <- unique(sub(pattern = "g",replacement = "t",all_genes_l))
gene_local_pos <- c()
for (i in 1:length(t_local)) {
  gene_local_pos <- c(gene_local_pos,grep(t_local[i],names(whole_gene_data)))
}
gene_local_gc <- c()
for (i in 1:length(gene_local_pos)) {
  gene_local_gc <- c(gene_local_gc,mean(GC1( whole_gene_data[[gene_local_pos[[i]]]]), GC2(whole_gene_data[[gene_local_pos[[i]]]]), GC3(whole_gene_data[[gene_local_pos[[i]]]])))
}

gene_local_gc
mean(gene_local_gc)
#####poisson (mean:0.494)
t_poisson <- unique(sub(pattern = "g",replacement = "t",all_genes_p))
gene_poisson_pos <- c()
for (i in 1:length(t_poisson)) {
  gene_poisson_pos <- c(gene_poisson_pos,grep(t_poisson[i],names(whole_gene_data)))
}
gene_poisson_gc <- c()
for (i in 1:length(gene_poisson_pos)) {
  gene_poisson_gc <- c(gene_poisson_gc,mean(GC1( whole_gene_data[[gene_poisson_pos[[i]]]]), GC2(whole_gene_data[[gene_poisson_pos[[i]]]]), GC3(whole_gene_data[[gene_poisson_pos[[i]]]])))
}
mean(gene_poisson_gc)
###其他基因之gc content (0.485)
gene_other_pos <- c(1:45019)
gene_other_pos <- gene_other_pos[-c(gene_local_pos,gene_poisson_pos)]
gene_other_gc <- c()
for (i in 1:length(gene_other_pos)) {
  gene_other_gc <- c(gene_other_gc,mean(GC1( whole_gene_data[[gene_other_pos[[i]]]]), GC2(whole_gene_data[[gene_other_pos[[i]]]]), GC3(whole_gene_data[[gene_other_pos[[i]]]])))
}

###重疊之基因(0.502)
intersect(all_genes_p,all_genes_l)
t_overlap <- unique(sub(pattern = "g",replacement = "t",intersect(all_genes_p,all_genes_l)))
gene_overlap_pos <- c()
for (i in 1:length(t_overlap)) {
  gene_overlap_pos <- c(gene_overlap_pos,grep(t_overlap[i],names(whole_gene_data)))
}
gene_overlap_gc <- c()
for (i in 1:length(gene_overlap_pos)) {
  gene_overlap_gc <- c(gene_overlap_gc,mean(GC1( whole_gene_data[[gene_overlap_pos[[i]]]]), GC2(whole_gene_data[[gene_overlap_pos[[i]]]]), GC3(whole_gene_data[[gene_overlap_pos[[i]]]])))
}

mean(gene_overlap_gc)
gc_aov_df <- data.frame(gc = c(gene_local_gc,gene_poisson_gc,gene_other_gc,gene_overlap_gc),Type = rep(c("L","P","Non-hotspot","Overlap"),times = c(145,1063,43877,66)) )

model <- aov(formula = gc~Type,data = gc_aov_df)
summary(model)
TukeyHSD(model,ordered = T)
write.csv(gc_aov_df,file = "C:/Users/lingo1st/Dropbox/林冠瑜/gc content/all gc content value.csv")
write.csv(all_genes,file ="C:/Users/lingo1st/Dropbox/林冠瑜/gc content/all_gene_poisson.csv" )
write.csv(all_genes_l,file ="C:/Users/lingo1st/Dropbox/林冠瑜/gc content/all_gene_local.csv" )
###plot

p1 <- ggplot(data = gc_aov_df,aes(x = Type,y = gc,fill = Type))
p1 + geom_boxplot(outliers = F)+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") + 
  scale_x_discrete(limits = c("L","P","Overlap","Non-hotspot"))+ylab("GC content")+
  geom_signif(comparisons = list(c("Overlap", "L"),
                                 c("Overlap","Non-hotspot")),
              map_signif_level = TRUE,
              y_position = c(0.75,0.85))
