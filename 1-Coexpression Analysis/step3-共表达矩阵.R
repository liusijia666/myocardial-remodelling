##step3:构建共表达矩阵

#计算基因间的邻接性，根据邻接性计算基因间的相似性，然后推出基因间
#的相异性系数，并据此得到基因间的系统聚类树。然后按照混合动态剪切树
#的标准，设置每个基因模块最少的基因数目为30。
# 计算资源允许的情况下最好放在一个block里面。
# power: 上一步计算的软阈值
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
rm(list = ls())
library(WGCNA)
load(file = "step2.rdata")
nGenes = ncol(datExpr)
net = blockwiseModules(datExpr, power = 12,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,corType = "pearson",
                       saveTOMFileBase = "myocardial remodelling TOM",
                       verbose = 3,maxBlockSize = nGenes)
table(net$colors)
save(net,file = "step3.rdata")
# 根据模块中基因数目的多少，一般降序排列，1-最大模块数。