rm(list = ls())

library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)

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

#（2）对上调/下调/所有差异基因进行富集分析
#注意这里又有个F
if(T){
  kegg_up <- enrichKEGG(gene = gene_up,
                       organism = 'hsa',
                       pvalueCutoff = 0.05)
  kegg_down <- enrichKEGG(gene = gene_down,
                   organism = 'hsa',
                   pvalueCutoff = 0.05)
}

up <- data.frame(kegg_up@result)
up <- up[up$pvalue<0.05,]
up <- arrange(up,desc(Count))
#up <- up[1:10,]


down <- data.frame(kegg_down@result)
down <- down[down$pvalue<0.05,]
down <- arrange(down,desc(Count))
#down <- down[1:10,]


up$change <- "up"
down$change <- "down"
KEGG <- rbind(up,down)


write.table(KEGG,file = "KEGG.txt",   
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t") 
