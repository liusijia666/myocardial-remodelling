##step4:模块可视化(层级聚类树展示各个模块)
rm(list = ls())
load('step1.rdata')
load('step3.rdata')

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)#将标绘数字标签转换为颜色
table(moduleColors)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#将所有的基因分配到相应的模块
#不同的颜色来代表那些所有的模块，其中灰色默认是无法归类于任
#何模块的那些基因，如果灰色模块里面的基因太多，那么对表达矩阵
#挑选基因的步骤可能就不成功。

