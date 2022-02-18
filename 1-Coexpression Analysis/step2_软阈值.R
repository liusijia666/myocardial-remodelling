##step2:计算软阈值

rm(list = ls())
library("WGCNA")
load('step1.rdata')

#通过绘制样本聚类查看分组信息和有无异常样本，在后续分析中是否
#需要剔除
sampleTree = hclust(dist(datExpr),method = "average")
plot(sampleTree,main = "Sample clustering to detect outliers",sub = "",
     xlab = "")

#确定软阈值(在得到两个基因的相关系数值后，利用软阈值来决定
#两个基因在构建网络时是否连边)
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#多软阈值功率无标度拓扑分析，选择合适的软阈值功率用于网络建设。
par(mfrow = c(1,2))#一页两图
cex1 = 0.75
plot(sft$fitIndices$Power,
     -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices$Power,
     -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     labels = powers,
     cex = cex1,
     col = "red")
abline(h = 0.85, col = "RED")
#横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
#网络越符合无标度特征 (non-scale)
# Mean connectivity
plot(sft$fitIndices$Power,
     sft$fitIndices$mean.k.,
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices$Power,
     sft$fitIndices$mean.k.,
     labels = powers,
     cex = cex1,
     col = "red")

best_beta = sft$powerEstimate
best_beta
# 软阈值为12
#若R^2达不到0.8，或平均连接度没有小于100时，前面数据未剔除成功

save(datExpr,sampleTree,sft,file = "step2.rdata")

