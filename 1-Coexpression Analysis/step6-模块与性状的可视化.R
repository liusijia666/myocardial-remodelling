##step5:模块与性状的关系
#通过模块与各种表型的相关系数，挑选模块进行下游分析
rm(list = ls())
library(WGCNA)
load('step1.rdata')
load('step3.rdata')

#明确样本数和基因数
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

design = model.matrix(~0 + group_list)
colnames(design) = levels(group_list)
#基因与基因之间的相关性系数矩阵(0/1矩阵)，用阈值判断基因相关与否
moduleColors <- labels2colors(net$colors)#数字标签的向量转换为与标签对应的颜色向量。
# 用颜色标签重新计算MEs
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
#模块特征基因（module eigengene）被定义为给定模块的第一个主成分，
#是该模块最具代表性的表达模式。
MEs = orderMEs(MEs0) #不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# 显示相关性和它们的p值
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8,14.5,1.5,1.5))#保证图形导出完整
# 在热图图中显示相关值
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = c("control","HCM","HF"),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


save(MEs,nGenes,design,nSamples,moduleColors,file = 'step6.rdata')
