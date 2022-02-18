rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
library(tidyr)
gse = "GSE32453"    #肥厚型梗死性心肌病

eSet <- getGEO(gse, 
               destdir = '.',
               getGPL = F)
#(1)提取表达矩阵exp
exp <- exprs(eSet[[1]])
exp[1:4,1:4]
exp = log2(exp+1)
exp = na.omit(exp)
table(!complete.cases(exp))

#(2)提取临床信息
pd0 <- pData(eSet[[1]])
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd0),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd0),colnames(exp))]

#(4)提取芯片平台编号
gpl <- eSet[[1]]@annotation

#1.group_list(实验分组)和ids(芯片注释)

library(stringr)
library(dplyr)
#2.探针注释----
temp = read.table("GPL6104-11576.txt",
                  header = T,skip = 26,
                  sep = "\t",quote = "",fill = T)

rowname <- data.frame(ID=rownames(exp))
deg <- left_join(rowname,temp,by='ID')
ex0 <- mutate(data.frame(symbol=deg$Symbol,exp)) %>% 
  filter(!duplicated(.$symbol))
rownames(ex0) <- ex0[,1]

save(pd0,ex0,file = "step0.Rdata")   #肥厚型梗死性心肌病_GSE32453
