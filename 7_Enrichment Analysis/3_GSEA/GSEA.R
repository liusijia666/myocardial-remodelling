###GSEA
rm(list = ls())
library(GSEABase)
library(fgsea)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(msigdbr)  #提供MSigdb数据库基因集
load("../../2_Difference Analysis/Difference_Analysis.Rdata")

##准备基因集矩阵  基因集富集分析

common_gene <- read.table("../../4_Common genes/common_gene.txt",   
                     header = F,sep = "\t",quote = "")
ex_deg <- ex_deg[common_gene$V1,]
gene_exp <- ex_deg[,c(1,7)]
gene_exp <- gene_exp[order(gene_exp$logFC,decreasing = T),]
GSEA_exp <-  gene_exp$logFC
names(GSEA_exp) = gene_exp$symbol
GSEA_exp[1:4]
a <- gmtPathways("c5.go.bp.v7.4.symbols.gmt")
##GSEA分析
fgseaRes <- fgsea(pathways = a, 
                  stats = GSEA_exp,
                  minSize=1,
                  maxSize=500,
                  nperm=10000)  #越大越稳定

class(fgseaRes)
fgseaRes <- fgseaRes[order(fgseaRes$pval),]
sum(fgseaRes[, fgseaRes$pval < 0.05])
aaa <- subset(fgseaRes, fgseaRes$pval < 0.05)

topPathwaysUp <- aaa[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- aaa[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(a[topPathways], GSEA_exp, aaa, 
              gseaParam = 0.5)
