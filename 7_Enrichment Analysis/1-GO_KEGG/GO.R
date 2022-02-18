rm(list = ls())

library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)

load('../../2_Difference Analysis/Difference_Analysis.Rdata')
common_gene <- read.table("../../4_Common genes/common_gene.txt",   
                          header = F,sep = "\t",quote = "")
deg_gene <- ex_deg[common_gene$V1,] 
deg_gene <- deg_gene[,c(7,8)]

#(1)输入数据

deg_gene <- bitr(deg_gene$symbol,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = "org.Hs.eg.db")
ex_deg <- ex_deg[deg_gene$SYMBOL,]
deg_gene$change <- ex_deg$change

gene_up = deg_gene[deg_gene$change == 'up','ENTREZID'] 
gene_down = deg_gene[deg_gene$change == 'down','ENTREZID'] 

if(T){
  #生物过程
  go_up <- enrichGO(gene = gene_up,
                    OrgDb= org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01,
                    readable = TRUE)
  
  go_down <- enrichGO(gene = gene_down,
                      OrgDb= org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      readable = TRUE)
  
}


#(3)可视化
#气泡图
up <- go_up@result
up <- up[order(up[,3],decreasing=T),]
up <- subset(up, up$pvalue < 0.05)
# summary(up$pvalue)
# summary(up$p.adjust)
up$p <- -log10(up$p.adjust)

down <- go_down@result
down <- down[order(down[,3],decreasing=T),]
down <- subset(down, down$pvalue < 0.05)
# summary(up$pvalue)
# summary(up$p.adjust)
down$p <- log10(down$p.adjust)

up$change <- "up"
down$change <- "down"
GO <- rbind(up,down)

write.table(GO,file = "GO.txt",   
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t") 
