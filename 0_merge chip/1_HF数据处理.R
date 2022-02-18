#数据下载
rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
library(stringr)
load("step1output.Rdata")      #已经下载好的GSE1145信息
#加载已经处理好的GSE1145表达矩阵ex2和临床信息pd2

# gse = "GSE1145"    #心衰
# eSet <- getGEO(gse,
#                destdir = '.',
#                getGPL = F)
#(1)提取表达矩阵exp
ex2 <- exprs(eSet[[1]])
ex2[1:4,1:4]
ex2 = log2(ex2+1)

table(!complete.cases(ex2))
#ex2 = na.omit(ex2)
#(2)提取临床信息
pd2 <- pData(eSet[[1]])

#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd2),colnames(ex2));p
if(!p) exp = exp[,match(rownames(pd2),colnames(exp))]

# 挑选出来自HCM的HF样本和正常样本
gr = pd2$title
k1 = str_detect(gr,"PGA_PA-N");table(k1)
k2 = str_detect(gr,"PGA-Hs_H");table(k2)
ex2 = ex2[,rownames(pd2)[k1|k2]]
pd2 = pd2[k1|k2,]
pd2$title <- c(rep("normal",11),rep("HF",5)) 

#(4)提取芯片平台编号
gpl <- eSet[[1]]@annotation

#1.group_list(实验分组)和ids(芯片注释)

library(stringr)
library(dplyr)
#2.探针注释----
temp = read.table("GPL570-55999.txt",
                  header = T,skip = 16,
                  sep = "\t",quote = "",fill = T)

rowname <- data.frame(ID=rownames(ex2))
deg <- left_join(rowname,temp,by='ID')
ex2 <- mutate(data.frame(symbol=deg$Gene.Symbol,ex2)) %>% 
  filter(!is.na(.$symbol)) %>% 
  filter(!(.$symbol) == "",) %>%   
  filter(!duplicated(.$symbol)) 
rownames(ex2) <- ex2[,1]

pd1 <- pd2
ex1 <- ex2
save(pd1,ex1,file = "step1.Rdata")   #心衰
