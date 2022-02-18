rm(list = ls()) 
load("../0_merge chip/merge.Rdata")

#差异分析，用limma包来做
#需要表达矩阵和group_list，不需要改

library(limma)
library(stringr)
library(dplyr)
#挑选出肥厚型心肌病与心衰
gr = outPd$title
k1 = str_detect(gr,"HCM");table(k1)
k2 = str_detect(gr,"HF");table(k2)
outPd = outPd[k1|k2,]

outTab = outTab[,rownames(outPd)]

group_list <- factor(c(rep("HCM",8),rep("HF",5)))
design=model.matrix(~group_list)
fit=lmFit(outTab,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)

#为deg数据框添加几列

deg <- mutate(deg,symbol=rownames(deg))
head(deg)

#2.加symbol列，火山图要用

#按照symbol列去重复

#3.加change列,标记上下调基因
logFC_t=0.5
P.Value_t = 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
change = ifelse(k1,
                "down",
                ifelse(k2,
                       "up",
                       "stable"))
table(change)
deg <- mutate(deg,change)
ex_deg <- deg[k1|k2,]
save(outTab,group_list,deg,logFC_t,P.Value_t,ex_deg,file = "Difference_Analysis.Rdata")
