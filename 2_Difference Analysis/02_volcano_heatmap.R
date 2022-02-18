rm(list = ls()) 
load(file = "Difference_Analysis.Rdata")

#1.火山图----
library(dplyr)
library(ggplot2)
dat  = deg

if(T) {
  x1 = dat[order(dat$logFC,decreasing = T),] %>% 
    filter(change == "up") %>% 
    head(2)
  x2 = dat[order(dat$logFC,decreasing = F),] %>% 
    filter(change == "down") %>% 
    head(2)
  for_label = rbind(x1,x2)
}

p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p
p <- p +
  #geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )
p

ggsave(p,filename = "volcan_plot.pdf",width = 15,height = 10,limitsize = FALSE)

#2.差异基因热图----
ex2 <- outTab
#ex2 <- ex2[,-13]
group_list <- factor(c(rep("HCM",8),rep("HF",5)))
deg1 <- data.frame(symbol=deg$symbol,logfc=abs(deg$logFC)) %>% 
  arrange(.by_group = T,desc(.$logfc))
deg2 <- deg1[1:30,]$symbol
group_list = group_list[order(group_list)]

n <- ex2[deg2,]
dim(n)

#作热图
library(pheatmap)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(n) 
heatmap_plot <- pheatmap(n,show_colnames =F,
                         show_rownames = T,
                         cluster_rows  = T,
                         cluster_cols = F,
                         scale = 'row',
                         annotation_col=annotation_col)
heatmap_plot

ggsave(heatmap_plot,filename = "heatmap.pdf",width = 15,height = 10,limitsize = FALSE)
