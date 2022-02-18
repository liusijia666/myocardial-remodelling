##Step5:基因和筛选模型的可视化（TOM图）
#直接load计算好的TOM结果,否则需要再计算一遍，比较耗费时间
rm(list = ls())
load('step2.rdata')
load('step3.rdata')

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)

load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
plotTOM = dissTOM^7
diag(plotTOM) = NA
geneTree = net$dendrograms[[1]]

TOMplot(plotTOM,geneTree, moduleColors[net$blockGenes[[1]]],
        main = "Network heatmap plot, all genes",
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))


# TOM图是指基因和筛选模型的可视化表示。
# WGCNA认为基因之间的简单的相关性不足以计算共表达，所以它利用邻近矩
# 阵，又计算了一个新的邻近矩阵。一般来说，TOM就是WGCNA分析的最终结果，
#后续的只是对TOM的下游注释。


