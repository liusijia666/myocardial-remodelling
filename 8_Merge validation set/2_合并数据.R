##心肌肥厚合并数据
rm(list = ls())
load("step0output.Rdata")
load("step1output.Rdata")
library(stringr)
library(dplyr)
ex1 <- data.frame(ex1)
ex1$symbol <- rownames(ex1)

a <- ex1$symbol       #韦恩图
b <- ex2$symbol
write.table(a,'GSE36961_genes.txt',
            row.names = F,
            col.names = F,
            sep = '\t',quote = F)
write.table(b,'GSE2656_genes.txt',
            row.names = F,
            col.names = F,
            sep = '\t',quote = F)

merge_eset=inner_join(ex1,ex2,by="symbol")
rownames(merge_eset) <- merge_eset$symbol
merge_eset <- merge_eset[,-which(colnames(merge_eset)=="symbol")]
dim(merge_eset)

##去除批次效应
library(sva)
library(stringr)
data <- as.matrix(merge_eset)

batchType=c(rep(1,145),rep(2,12))
modType=c(rep("Control",39),rep("HCM",106),
          rep("Control",4),rep("HF",8)
)
mod=model.matrix(~as.factor(modType))
#mod=as.factor(modType)
outTab2=data.frame(ComBat(data,batchType,mod,par.prior = TRUE))
dim(outTab2)
boxplot(outTab2)

a <- c(40:145,150:157)
outTab2 <- outTab2[,a]
save(outTab2,file = "merge.Rdata")
