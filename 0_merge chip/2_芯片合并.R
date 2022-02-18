##心衰合并数据
rm(list = ls())
load("step0.Rdata")
load("step1.Rdata")
library(dplyr)
library(sva)
library(ggplot2)

a <- ex0$symbol        #韦恩图
b <- ex1$symbol
write.table(a,'GSE32453_genes.txt',
            row.names = F,
            col.names = F,
            sep = '\t',quote = F)
write.table(b,'GSE1145_genes.txt',
            row.names = F,
            col.names = F,
            sep = '\t',quote = F)

merge_eset=inner_join(ex0,ex1,by="symbol")
rownames(merge_eset) <- merge_eset$symbol
merge_eset <- merge_eset[,-1]
dim(merge_eset)

##去除批次效应

data <- as.matrix(merge_eset)

batchType=c(rep(1,13),rep(2,16))
modType=c(rep("HCM",8),rep("control",5),
          rep("control",11),rep("HF",5)
)
mod=model.matrix(~as.factor(modType))
#mod=as.factor(modType)
outTab=data.frame(ComBat(data,batchType,mod,par.prior = TRUE))
dim(outTab)
boxplot(outTab)

pd0=pd0[,c(1,2)]
pd1=pd1[,c(1,2)]
outPd=rbind(pd0,pd1)

p = identical(rownames(outPd),colnames(outTab));p
if(!p) outTab = outTab[,match(rownames(outPd),colnames(outTab))]

save(outPd,outTab,file = "merge.Rdata")
