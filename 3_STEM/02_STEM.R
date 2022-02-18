rm(list = ls())
library(Mfuzz)

#该示例中，行为基因或蛋白名称，列为时间样本（按时间顺序提前排列好，若存在生物学重复需提前取均值）
protein <- read.delim('stem_input.txt', row.names = 1, check.names = FALSE)
protein <- as.matrix(protein)

## Mfuzz聚类时要求是一个ExpressionSet类型的对象，
#所以需要先用表达量构建这样一个对象。
mfuzz_class <- new('ExpressionSet',exprs = protein)

#剔除异常值
mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
#NA替换缺失值
mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
#根据标准差去除样本间差异太小的基因
mfuzz_class <- filter.std(mfuzz_class, min.std = 0)

### 2. 标准化：聚类时需要用一个数值来表征不同基因间的距离，Mfuzz中采用的是欧式距离，
# 由于普通欧式距离的定义没有考虑不同维度间量纲的不同，所以需要先进行标准化
mfuzz_class <- standardise(mfuzz_class)

#Mfuzz 基于 fuzzy c-means 的算法进行聚类
#需手动定义目标聚类群的个数
#需要设定随机数种子，以避免再次运行时获得不同的结果
set.seed(123)
cluster_num <- 4
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))

#作图
#time.labels 参数设置时间轴，需要和原基因表达数据集中的列对应
#颜色、线宽、坐标轴、字体等细节也可以添加其他参数调整，此处略，详见函数帮助
mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, 
            mfrow=c(2,2),
            time.labels = colnames(protein),
            centre=T,
            centre.col="red")

#最后，提取所有蛋白所属的聚类群，并和它们的原始表达值整合在一起
protein_cluster <- mfuzz_cluster$cluster
protein_cluster <- cbind(protein[names(protein_cluster), ], protein_cluster)
head(protein_cluster)
# protein_cluster <- data.frame(protein_cluster)
# protein_cluster$aa <- rownames(protein_cluster)
# CD86、IL13RA1、MYB在cluster4,JAK2、PLCL1在cluster3
write.table(protein_cluster, 'protein_cluster.txt', sep = '\t', col.names = NA, quote = FALSE)