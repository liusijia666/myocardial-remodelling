rm(list = ls())
setwd("CIBERSORT/")
source("CIBERSORT.R")
load("../merge chip/merge.Rdata")

write.table(outTab,file = "outTab.txt",
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")
#https://www.nature.com/articles/nmeth.3337#MOESM207
#下载Supplementry table 1（LM22）
#https://rdrr.io/github/singha53/amritr/src/R/supportFunc_cibersort.R
#下载CIBERSORT.R

# Define LM22 file

TME.results = CIBERSORT("LM22.txt", "outTab.txt", perm = 1000, QN = TRUE)

# output CIBERSORT results
write.table(TME.results, "TME.results.output.txt", 
            sep = "\t", row.names = T, col.names = T, quote = F)

## 可视化
# 22种免疫细胞在所有样本种的比例
# boxplot
library(ggpubr)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(dplyr)

re <- TME.results[,-(23:25)]
dat <- re %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

## 可视化
#无序箱线图
p <- ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) +
  geom_boxplot(outlier.shape = 21,color = "black") +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = mypalette(22))
p
ggsave(p,filename = "CIBERSORT.pdf",width = 12,
       height = 10)


#按照顺序排列箱线图
a = dat %>% 
  group_by(Cell_type) %>% 
  summarise(m = median(Proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(Cell_type)

dat$Cell_type = factor(dat$Cell_type,levels = a)

p <- ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = mypalette(22))
p
ggsave(p,filename = "CIBERSORT2.pdf",width = 12,
       height = 10)


# 那就做个比较
a <- c(rep("HCM",8),rep("normal",16),rep("HF",5))
b <- rep(a,22)
dat$Group <- factor(b,levels = c("normal","HCM","HF"))

library(ggpubr)
q <- ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,3,4)])+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")
q
ggsave(q,filename = "CIBERSORT_sample.pdf",width = 15,
       height = 10)
