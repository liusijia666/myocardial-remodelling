#数据下载
rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072*2)
gse = "GSE36961"#心肌肥厚
eSet <- getGEO(gse, 
               destdir = '.',
               getGPL = F)
#(1)提取表达矩阵exp
ex1 <- exprs(eSet[[1]])
ex1[1:4,1:4]
table(!complete.cases(ex1))

#(2)提取临床信息
pd1 <- pData(eSet[[1]])
#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd1),colnames(ex1));p
if(!p) ex1 = ex1[,match(rownames(pd1),colnames(ex1))]

#(4)提取芯片平台编号
gpl <- eSet[[1]]@annotation

save(pd1,ex1,file = "step0output.Rdata")    #心肌肥厚


#数据下载
#----
#数据下载
rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
library(tidyr)
#(1)读入表达矩阵exp
exp <- read.table("GSE2656_series_matrix.txt", #心衰
                  header = T, sep = "", quote = "\"'",
                  fill = T,skip = 65,row.names = 1)
exp[1:4,1:4]
exp = na.omit(exp)
table(!complete.cases(exp))
exp <- exp[,c(12,13,14,15,33,34,35,36,37,38,39,40)]    
#前四个是来源于HCM的HF患者，后八个是正常心脏样本
#group <- c(rep("HF",4),rep("normal",8))
#2.探针注释----
temp = read.table("GSE2656.txt",
                  header = T,skip = 6,
                  sep = "\t",quote = "",fill = T)
## 利用UCSC进行GB_ACC注释

library(dplyr)
anno=read.table("refGene.txt",sep="\t")
#提取第2列和第13列，分别为GB_ACC和gene symbol
anno=anno[,c(2,13)]
#给这两列命名
names(anno)=c("GB_ACC","symbol")
#去重复
anno=unique(anno)
#将GB_ACC设置成行名
rownames(anno)=anno$GB_ACC

## 添加symbol列
#通过GB_ACC获取对应的gene symbol
symbol=as.character(anno[temp$GB_ACC,"symbol"])
#将NA转换成空
#symbol[is.na(symbol)]=""
#将gene symbol这一列加入到原来的探针注释文件中
temp$symbol=symbol

rowname <- data.frame(ID=rownames(exp))
temp$ID <- as.character(temp$ID)
deg <- left_join(rowname,temp,by='ID')
ex1 <- mutate(data.frame(symbol=deg$symbol,exp)) %>% 
  filter(!duplicated(.$symbol))
ex2 = na.omit(ex1)
table(!complete.cases(ex2))

save(ex2,file = "step1output.Rdata")
