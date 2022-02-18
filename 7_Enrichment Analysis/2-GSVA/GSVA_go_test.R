###GSVA
rm(list = ls())

library(msigdbr)
library(GSVA)
library(limma)
library(ggplot2)
library(GSEABase)
library(stringr)
library(dplyr)

load("../../2_Difference Analysis/Difference_Analysis.Rdata")
##准备表达矩阵

commom_gene <- read.table("../../4_Common genes/common_gene.txt",   
                          header = F,sep = "\t",quote = "")
gene_exp <- outTab[commom_gene$V1,]
##准备基因集矩阵
#GSVA富集数据库来源
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

#---------------------------------------------
#up pathway
go <- read.table("../1-GO_KEGG/GO.txt",   
                 header = T,sep = "\t",quote = "")

go_name <- paste('GOBP',go[,2]) %>%
  str_replace_all(' ','_') %>% 
  str_replace_all('-','_') %>% 
  toupper()    #将通路名字与提供基因集的数据库格式调整成一样

b1=msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)

go_list=list()
for(i in go_name){
  go_df=b1[b1$gs_name==i,]
  go_list[i]=list(go_df$gene_symbol)
}

##GSVA构建新的矩阵

go_Es = gsva(as.matrix(gene_exp),go_list,min.sz=1,
                kcdf="Gaussian",method="gsva")

##limma差异分析
#p = identical(rownames(outPd),colnames(go_Es));p
group_list <- factor(c(rep("HCM",8),rep("HF",5)))
design=model.matrix(~group_list)
fit=lmFit(go_Es,design)
fit=eBayes(fit)
GO=topTable(fit,coef=2,number = Inf)

#选取GSVA排名前10的上调通路、下调通路
GO <- GO[order(GO$logFC,decreasing = T),]
GO_down <- tail(GO,10)
GO_up <- head(GO,10)

#----------------------------------------------
#up pathway
#查看10条通路在GO富集中的结果
go$pathway <- go_name
GO_up$pathway <- rownames(GO_up)
GO_up <- GO_up[,c(1,7)]
GO_Bubble <- inner_join(go,GO_up,by="pathway")
GO_Bubble <- GO_Bubble[,-13]
#GO_Bubble <- filter(up,up$pathway==GO_up$pathway)??

#GO富集气泡图可视化
p <- ggplot(GO_Bubble,aes(x=GeneRatio,y=Description))+
  geom_point(aes(size=Count,color=-1*log10(pvalue)))+
  scale_colour_gradient(low="blue",high="red")+
  labs(
    color=expression(-log[10](P.value)),
    size="Genenumber",
    x="",
    y="Pathway name",
    title="upregulated in pathway")+
  theme_bw()+
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks=element_blank(),
    legend.position = "left",
    plot.title = element_text(colour = "red"))

p

#ggsave(filename = "GO_up_Bubble.pdf",width = 15,height = 5)

##10条通路GSVA柱状图
q <- ggplot(GO_up,aes(x = reorder(GO_up$pathway,GO_up$logFC),y = logFC,fill = pathway))+
  geom_bar(stat = "identity",colour="white",width = 0.8,position = position_dodge(0.7))+
  xlab("")+
  ylab("logFC of GSVA score,heart fail versus HCM")+
  coord_flip()+
  theme_bw()+
  #theme(plot.margin = unit(c(1,11,1,10),"lines"),complete = TRUE)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ #去掉网格
  #theme(axis.text.y = element_blank())+  #删除刻度标签
  scale_x_discrete(position = "top")+
  theme(axis.line.x = element_line(colour = "black"),
        axis.ticks.y = element_blank(),    #去掉刻度
        panel.border = element_blank(),   #去掉边框
        legend.position = 'none')+
  #scale_y_discrete(position = "right")+
  geom_hline(aes(yintercept=0),linetype="dashed",size=1)
q

#ggsave(filename = "GO_up_bar.pdf",width = 15,height = 5)

library(patchwork)
a <- p+q

ggsave(filename = "GO_up.pdf",width = 15,height = 5)

#---------------------------------------------
##down pathway

GO_down$pathway <- rownames(GO_down)
GO_down <- GO_down[order(GO_down[,1],decreasing=T),]
GO_down <- GO_down[,c(1,7)]

#查看10条通路在GO富集中的结果
#down$pathway <- go_name
GO_Bubble <- inner_join(go,GO_down,by="pathway")
GO_Bubble <- GO_Bubble[,-13]
#GO_Bubble <- filter(up,up$pathway==GO_up$pathway)??

#GO富集气泡图可视化
aa <- ggplot(GO_Bubble,aes(x=GeneRatio,y=Description))+
  geom_point(aes(size=Count,color=-1*log10(pvalue)))+
  scale_colour_gradient(low="blue",high="red")+
  labs(
    color=expression(-log[10](P.value)),
    size="Genenumber",
    x="",
    y="Pathway name",
    title="downregulated in pathway")+
  theme_bw()+
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks=element_blank(),
    #legend.position = "left",
    plot.title = element_text(colour = "blue",vjust = "-115",hjust = "0.5")
  )

aa

#ggsave(filename = "GO_up_Bubble.pdf",width = 15,height = 5)

##10条通路GSVA柱状图
bb <- ggplot(GO_down,aes(x = reorder(pathway,logFC),y = logFC,fill = pathway))+
  geom_bar(stat = "identity",colour="white",width = 0.8,position = position_dodge(0.7))+
  xlab("")+
  ylab("logFC of GSVA score,heart fail versus HCM")+
  coord_flip()+
  theme_bw()+
  #theme(plot.margin = unit(c(1,11,1,10),"lines"),complete = TRUE)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ #去掉网格
  #theme(axis.text.y = element_blank())+  #删除刻度标签
  scale_x_discrete(position = "top")+
  theme(axis.line.x = element_line(colour = "black"),
        axis.ticks.y = element_blank(),    #去掉刻度
        panel.border = element_blank(),   #去掉边框
        legend.position = 'none')+
  #scale_y_discrete(position = "right")+
  scale_x_discrete(position = "bottom")+
  geom_hline(aes(yintercept=0),linetype="dashed",size=1)
bb

#ggsave(filename = "GO_up_bar.pdf",width = 15,height = 5)

library(patchwork)
cc <- bb+aa

ggsave(filename = "GO_down.pdf",width = 15,height = 5)
