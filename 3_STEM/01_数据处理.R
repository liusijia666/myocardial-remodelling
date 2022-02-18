###Stem输入文件
#Expression Daa Info
rm(list = ls())
load(file = "../0_merge chip/merge.Rdata")
deg_gene <- outTab
#由于STEM只能处理相同个数的生物重复学样本，因此直接将重复样本
#按照STEM的算法取中位数
HCM <- apply(deg_gene[,1:8],1,median)
control <- apply(deg_gene[,9:24],1,median)
heart_failure <- apply(deg_gene[,25:29],1,median)
stem_input <- data.frame(cbind(control,HCM,heart_failure))
gene_name <- rownames(stem_input)
stem_input <- cbind(gene_name,stem_input)
colnames(stem_input) <- c('gene_name','control','HCM', 'heart_failure')
write.table(stem_input,file = "stem_input.txt",
            quote = F,
            sep = "\t",
            col.names = T,
            row.names = F)


##可视化
##绘制出4个Hub genes 基因表达趋势箱线图
Hub_genes <- c("IL13RA1","CD86","JAK2","MYB","PLCE1")
#protein <- data.frame(protein)
dat <- deg_gene[which(rownames(deg_gene) %in% Hub_genes),]
dat <- data.frame(t(dat))
group <- factor(c(rep("HCM",8),rep("control",16),rep("heart_failure",5)))
dat$group <- group

library(ggplot2)
p = list()
for(i in 1:(ncol(dat)-1)){
  p[[i]] = ggplot(data = dat,aes_string(x = "group",y=colnames(dat)[i]))+
    geom_boxplot(aes(color = group))+
    geom_jitter(aes(color = group))+
    stat_summary(fun=median, geom="line", aes(group=1))+
    theme_bw()
}
library(patchwork)
wrap_plots(p,nrow = 1,guides = "collect")
ggsave(filename = "box_plots.pdf",width = 15,height = 5)
