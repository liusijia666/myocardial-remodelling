#--------WGCNA---------------------

###Step1-数据准备
rm(list = ls())
options(stringsAsFactors = F)
library(WGCNA)
load("../0_merge chip/merge.Rdata")

exprSet <- outTab
datExpr = t(exprSet)
datExpr[1:4,1:4]

a <- goodSamplesGenes(datExpr,verbose = 3)
a$allOK

if(!a$allOK){
  if(sum(!a$goodGenes)>0)
    printFlush(paste("Removing genes:",paste(names(datExpr[!a$goodGenes],collapse = ","))))
  if(sum(!a$goodSamples)>0)
    printFlush(paste("Removing samples:",paste(rownames(datExpr)[!a$goodSamples],collapse = ",")))
  datExpr = datExpr[a$goodSamples,a$goodGenes]
}

group_list <- factor(c(rep("HCM",8),rep("control",5),
                       rep("control",11),rep("HF",5)))

save(datExpr,group_list,file = 'step1.rdata')
