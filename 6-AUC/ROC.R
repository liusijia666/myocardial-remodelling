rm(list= ls())
library(ggplot2)
library(pROC)
library(stringr) 
library(dplyr)
load("../2_Difference Analysis/Difference_Analysis.Rdata")

ppi_result <- read.csv("../5-PPI/degree_top100.csv",
                       header = TRUE, sep = ",",skip = 1)
#outTab$aa <- rownames(outTab)
group_list <- factor(c(rep("HCM",8),rep("HF",5)))
ppi_result <- outTab[ppi_result$Name,]
#rownames(ppi_result) <- ppi_result$symbol
ppi_result <- t(ppi_result)
type <- if_else(group_list=='HCM',0,1)
ppi_result <- data.frame(ppi_result,type=type)

d <- data.frame()
e <- data.frame()
for (i in 1:(ncol(ppi_result)-1)) {
  b <- roc(ppi_result$type,ppi_result[,i],ci=T)
  c <- matrix(b[["ci"]],ncol = 3,byrow = T)
  e <- rbind(e,b$auc)
  d <- rbind(d,c)
}
d <- d[,-2]
colnames(d) <- c('lower','upper')
colnames(e) <- "auc"
AUC <- data.frame(symbol=colnames(ppi_result[1:(ncol(ppi_result)-1)]),
                auc=e,
                d=d)

AUC <- AUC[order(-AUC[,2]),]
g <- subset(AUC, AUC$auc > 0.99)
Hub_gene <- g$symbol
write.table(Hub_gene,'hub_gene.txt',
            row.names = F,
            col.names = F,
            sep = '\t',quote = F)
