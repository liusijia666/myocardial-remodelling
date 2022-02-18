load("../2_Difference Analysis/Difference_Analysis.Rdata")
#差异分析的结果与WGCNA取交集
WGCNA <- read.table('../1-Coexpression Analysis/WGCNA.txt',    
                    header = T,
                    sep = "\t",quote = "")
common_gene <- intersect(rownames(ex_deg),WGCNA$gene.inModule.)
write.table(common_gene,file = "common_gene.txt",
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t")
