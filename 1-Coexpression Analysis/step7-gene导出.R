##Step8: 模块的导出
rm(list = ls())
library(plyr)
load('step1.rdata')
load('step3.rdata')
load('step6.rdata')

# Select gene
# 正相关性最高
module = "pink";
gene = colnames(datExpr) 
inModule = (moduleColors==module)
modgene1 =data.frame(gene[inModule]) 
# 正相关性最高
# module = "black";
# gene = colnames(datExpr)
# inModule = (moduleColors==module)
# modgene2 = data.frame(gene[inModule])
# 负相关性最高
module = "salmon";
gene = colnames(datExpr)
inModule = (moduleColors==module)
modgene3 = data.frame(gene[inModule])

modgene = rbind(modgene1,modgene3)
#modgene = modgene1
write.table(modgene,file = "WGCNA.txt",
            row.names = F,
            col.names = T,
            quote = F,
            sep = "\t")

