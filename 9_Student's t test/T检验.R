##批量T检验，确定Hub genes
rm(list = ls())
library(dplyr)
load("../8_Merge validation set/merge.Rdata")

#获取Hub genes
gene <- read.table("../6-AUC/hub_gene.txt")
hub_gene <- data.frame(gene=intersect(rownames(outTab2),gene$V1))
hub_gene <- t(outTab2[hub_gene$gene,])
hub_gene <- mutate(data.frame(sample=rownames(hub_gene),hub_gene))

#table(outPd2$title)
type=c(rep("HCM",106),rep("HF",8))

dir.create( 'Tjianyan')
setwd("Tjianyan/")
for (i in 2:ncol(hub_gene)){
  
  dat=as.numeric(unlist(lapply(hub_gene[i],function(x) strsplit(as.character(x),"%"))))
  
  filename=paste(names(hub_gene)[i],".png",sep="")
  
  png(file=filename, width = 1200, height = 1000)
  
  b=t.test(dat~type)
  
  txt=paste("p-value=",round(b$p.value[1],digits=20),sep="")
  
  plot(factor(type,levels=c("HCM","HF")),dat,ylab="percent(%)",main=names(hub_gene)[i],sub=txt,cex.axis=1.5,cex.sub=2,cex.main=2,cex.lab=1.5)
  
  dev.off()
}
  
