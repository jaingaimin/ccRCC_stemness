
setwd("~/vip39/TCGA/cancer_stemness/")

load("./cancerstem_sub.rds")

.libPaths()
.libPaths(c(
  #"~/vip39/R/x86_64-pc-linux-gnu-library/4.2/",
  "/refdir/Rlib_4.3",
  "~/vip39/R/x86_64-pc-linux-gnu-library/4.0/",
  "~/vip39/R/x86_64-pc-linux-gnu-library/4.1/",
  "/home/data/t030432/R/x86_64-pc-linux-gnu-library/4.3/"))
###免疫分析
##cibersort immune score
##自定义基因集做GSVA 绘制热图 带p值

options(stringsAsFactors = F)
library(GSVA)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(limma)
library(data.table)
pathway <- fread("./Stemness_gene_new.csv",header = T,data.table = F)
pathway <- as.data.frame(pathway)

pathway_list <- lapply(pathway, function(x) {
  unique(na.omit(x)) 
})

pathway_list <- vector("list",length(pathway))

for (i in seq_along(pathway)) {
  pathway_list[[i]] <- unique(na.omit(pathway[,i]))
}

names(pathway_list) <- colnames(pathway)

for (i in 1:36) {
  pathway_list[[i]] <- pathway_list[[i]][pathway_list[[i]]!=""]
}

setwd("~/TCGA/cancer_stemness/GSVA/")
head(KIRC_mRNA_fpkm)
gsva_matrix_BD <- gsva(as.matrix(KIRC_mRNA_fpkm), pathway_list,method='gsva',
                       kcdf='Gaussian',abs.ranking=TRUE)

gsva_matrix_japan <- gsva(as.matrix(japan_mRNA_fpkm), pathway_list,method='gsva',
                       kcdf='Gaussian',abs.ranking=TRUE)

save(gsva_matrix_BD,gsva_matrix_japan,file = "~/TCGA/cancer_stemness/GSVA/TCGA_JAPAN_GSVA.rds")
load("./GSVA/TCGA_JAPAN_GSVA.rds")
#write.csv(gsva_matrix_BD,file = "gsva_matrix_BD.csv")
save(sub,file = "./stemness_sub.rds")

load("./stemness_sub.rds")
####单因素分析预后分析
library(tidyverse,lib.loc = "/home/data/refdir/Rlib")
library(ggplot2)
library(ggstatsplot)
library(survival)
library(stringr)
library(viridis)
library(scales)
library(survminer)
dim(sub)
Coxoutput <- data.frame(OS=sub$EVENT,
                        OS.time=sub$OS
)
rownames(Coxoutput) <- sub$Sample
Coxoutput2 <- as.data.frame(t(KIRC_mRNA_fpkm[unique(RNA_gene),sub$Sample]))
Coxoutput <- cbind(Coxoutput,Coxoutput2)

realdata <- Coxoutput
realdata[1:3,1:6]
dir.create("COX")
setwd("./COX/")
Coxoutput=data.frame()
for(i in colnames(realdata[,3:ncol(realdata)])){
  cox <- coxph(Surv(OS.time, OS) ~ realdata[,i], data = realdata)
  coxSummary = summary(cox)
  Coxoutput=rbind(Coxoutput,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                                  z=coxSummary$coefficients[,"z"],
                                  pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                  lower=coxSummary$conf.int[,3],
                                  upper=coxSummary$conf.int[,4]))
}
for(i in c(2:6)){
  Coxoutput[,i] <- as.numeric(as.vector(Coxoutput[,i]))
}
Coxoutput <- arrange(Coxoutput,pvalue)  %>% #按照p值排序
  filter(pvalue < 0.01) 
dim(Coxoutput)
#[1] 121   6
#保存到文件
write.csv(Coxoutput,'cox_output.csv', row.names = F)


Coxoutput <- read.csv("cox_output.csv")
head(Coxoutput)

plotCoxoutput <- filter(Coxoutput,HR <=0.92 | HR>= 1.15)  #选择合适的top genes

#采用那篇nature的配色
ggplot(data=plotCoxoutput,aes(x=HR,y=gene,color=pvalue))+
  geom_errorbarh(aes(xmax=upper,xmin=lower),color="black",height=0,size=1.2)+
  geom_point(aes(x=HR,y=gene),size=3.5,shape=18)+ #画成菱形
  geom_vline(xintercept = 1,linetype='dashed',size=1.2)+
  scale_x_continuous(breaks = c(0.75,1,1.30))+
  coord_trans(x='log2')+ 
  ylab("Gene")+  #标签
  xlab("Hazard ratios of RNAmodify_gene in KIRC")+ 
  labs(color="P value",title ="Univariate Cox regression analysis" )+
  scale_color_viridis()+  #nature配色
  theme_ggstatsplot()+  #好看的主题，同原文一致
  theme(panel.grid =element_blank()) #去除网格线
ggsave('plot1.pdf',width = 8,height = 15)

#让误差线跟着变色
ggplot(data=plotCoxoutput,aes(x=HR,y=gene,color=pvalue))+
  geom_errorbarh(aes(xmax=upper,xmin=lower,color = pvalue),height=0,size=1.2)+
  geom_point(aes(x=HR,y=gene),size=3.5,shape=18)+   #画成菱形
  geom_vline(xintercept = 1,linetype='dashed',size=1.2)+
  scale_x_continuous(breaks = c(0.75,1,1.30))+
  coord_trans(x='log2')+ 
  ylab("Gene")+  #标签
  xlab("Hazard ratios of RNAmodify_gene in KIRC")+ 
  labs(color="P value",title ="" )+
  scale_color_viridis()+ 
  theme_ggstatsplot()+  #好看的主题，同原文一致
  theme(panel.grid =element_blank()) #去除网格线
ggsave('plot2.pdf',width = 8,height = 15)

#自己DIY主题和颜色
ggplot(data=plotCoxoutput,aes(x=HR,y=gene,color=pvalue))+
  geom_errorbarh(aes(xmax=upper,xmin=lower),color='black',height=0,size=1.2)+
  geom_point(aes(x=HR,y=gene),size=3.5,shape=18)+   #画成菱形
  geom_vline(xintercept = 1,linetype='dashed',size=1.2)+
  scale_x_continuous(breaks = c(0.75,1,1.30))+
  coord_trans(x='log2')+ 
  ylab("Gene")+  #标签
  xlab("Hazard ratios of RNAmodify_gene in KIRC")+ 
  labs(color="P value",title ="" )+
  scale_color_gradient2(low = muted("skyblue"),mid ="white",high =muted('pink'),midpoint = 0.025)+ #樱花色
  theme_bw(base_size = 12)+   #字体大小，格式
  theme(panel.grid =element_blank(),  #去除网格线，开始DIY主题
        axis.text.x = element_text(face="bold", color="black", size=9),    #各个字体大小
        axis.text.y = element_text(face="bold",  color="black", size=9),
        axis.title.x = element_text(face="bold", color="black", size=11),
        axis.title.y = element_text(face="bold",color="black", size=11),
        legend.text= element_text(face="bold", color="black", size=9),
        legend.title = element_text(face="bold", color="black", size=11),
        panel.border = element_rect(colour = 'black',size=1.4))   #边框粗细
ggsave('plot3.pdf',width = 8,height = 15)

save(RNA_gene,file = "~/TCGA/shuainao_model/senescencegene.rds")

##
##直接分型研究

library(ConsensusClusterPlus)
##提取出T1到T3
##这里面只取出sub生存大于30天的
expr <- as.data.frame(gsva_matrix_BD)
dim(expr)
expr <- expr[,sub$Sample]
#shuailaogene_need <- Coxoutput$gene
#save(shuailaogene_need,file = "~/TCGA/shuainao_model/shuailaogene.rds")
#[1]  Coxoutput
setwd("../")
dir.create('ConsensusCluster/')
results = ConsensusClusterPlus(as.matrix(expr),
                               maxK=9,
                               reps=100,
                               pItem=0.8,
                               pFeature=1,
                               title='ConsensusCluster/',
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123456,
                               plot="png")
Kvec = 2:9
x1 = 0.1; x2 = 0.9 
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") 
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
optK
#[1] 2
PAC <- as.data.frame(PAC)
PAC$K <- 2:9
library(ggplot2)
ggplot(PAC,aes(factor(K),PAC,group=1))+
  geom_line()+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_point(size=4,shape=21,color='darkred',fill='skyblue')+
  ylab('Proportion of ambiguous clustering')+
  xlab('Cluster number K')
library(export)
library(eoffice)
eoffice::topptx(file='ConsensusCluster/PAC.pptx',width=4.5,height=4.5)
topptx(file='ConsensusCluster/PAC.pptx',width=4.5,height=4.5)
eoffice::tofigure(file='ConsensusCluster/PAC.pdf')
ggsave(filename = "ConsensusCluster/PAC.pdf",width = 5,height = 5)

icl <- calcICL(results,title = 'ConsensusCluster/',plot = 'png')
clusterNum=2    
cluster=results[[clusterNum]][["consensusClass"]]

sub <- data.frame(Sample=names(cluster),Cluster=cluster)
sub$Cluster <- paste0('C',sub$Cluster)
table(sub$Cluster)
# C1  C2 
# 231 285
head(sub)

rownames(sub) <- sub$Sample
sub$EVENT <- as.vector(allclin[rownames(sub),"OS"])
sub$OS <- as.vector(allclin[rownames(sub),"OS.time"])

sub <- sub[sub$OS>30,]
library(survival)
library(survminer)
library(monocle)
library(DDRTree)
library(ggsci)
class(sub)
meta <- sub[,1:4]
head(sub)
colnames(meta) <- c("ID","cluster","event","time")
meta <- meta[,-1]
meta$time <- meta$time/365
head(meta)
sfit <- survfit(Surv(time, event)~cluster, data=meta)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
sfit <- survfit(Surv(OS,EVENT) ~ Cluster,data = sub)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
sfit <- survfit(Surv(OS,EVENT)~Cluster,data = sub)
sfit <- survfit(Surv(OS.time,EVENT)~Cluster,data = sub)

sub$PFI <- allclin[rownames(sub),"PFI"]
sub$PFI.time <- allclin[rownames(sub),"PFI.time"]
sfit <- survfit(Surv(PFI.time, PFI)~Cluster, data=sub)

mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black"))
ggsurvplot(sfit,
           data = meta,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("C1","C2"), 
           legend.title="cluster",
           xlab="Time (years)",
           ylab='Overall survival',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)
ggsave(filename = "OS_twogroup.pdf",height = 6,width = 6)

ggsurvplot(sfit,
           data = meta,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("C1","C2"), 
           legend.title="cluster",
           xlab="Time (days)",
           #ylab='Overall survival',
           ylab='Progression Free Survival',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)


###PFI
meta$PFI <- allclin[rownames(meta),"PFI"]
meta$PFI.time <- allclin[rownames(meta),"PFI.time"]/365
sfit <- survfit(Surv(PFI.time, PFI)~cluster, data=meta)
ggsurvplot(sfit,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("C1","C2"), 
           legend.title="cluster",
           xlab="Time (years)",
           ylab='Progression Free Survival',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)
ggsave(filename = "PFI_twogroup.pdf",height = 6,width = 6)
T1_meta <- meta%>%
  filter(pstage=="T1")
sfit <- survfit(Surv(PFI.time, PFI)~cluster, data=T1_meta)
sfit <- survfit(Surv(time, event)~cluster, data=T1_meta)
ggsurvplot(sfit,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("C1","C2"), 
           legend.title="cluster",
           xlab="Time (years)",
           ylab='Progression Free Survival',
           risk.table=TRUE,
           break.time.by = 2,
           risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)


T2_meta <- meta%>%
  filter(pstage=="T2")
sfit <- survfit(Surv(PFI.time, PFI)~cluster, data=T2_meta)
sfit <- survfit(Surv(time, event)~cluster, data=T2_meta)
ggsurvplot(sfit,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("C1","C2"), 
           legend.title="cluster",
           xlab="Time (years)",
           ylab='Progression Free Survival',
           risk.table=TRUE,
           break.time.by = 2,
           risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)
T3_meta <- meta%>%
  filter(pstage=="T3")
sfit <- survfit(Surv(PFI.time, PFI)~cluster, data=T3_meta)
sfit <- survfit(Surv(time, event)~cluster, data=T3_meta)
ggsurvplot(sfit,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("C1","C2"), 
           legend.title="cluster",
           xlab="Time (years)",
           ylab='Progression Free Survival',
           risk.table=TRUE,
           break.time.by = 2,
           risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)

T4_meta <- meta%>%
  filter(pstage=="T4")
sfit <- survfit(Surv(PFI.time, PFI)~cluster, data=T4_meta)
sfit <- survfit(Surv(time, event)~cluster, data=T4_meta)
ggsurvplot(sfit,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("C1","C2"), 
           legend.title="cluster",
           xlab="Time (years)",
           ylab='Progression Free Survival',
           risk.table=TRUE,
           break.time.by = 2,
           risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)
####绘制热图####
#加上正常组织三组绘图
dir.create("./complexheatmap")
setwd("/home/data/vip39/TCGA/cancer_stemness/complexheatmap/")
library(ComplexHeatmap)
library(stringr)
library(pheatmap)
library(gplots)
library(grid)
library(scico)
library(circlize)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
###准备两个数据 mygene_data 和 Subtype 分组信息

normal_id <- colnames(gsva_matrix_BD)[str_sub(colnames(gsva_matrix_BD),14,15)=="11"]
c(normal_id,sub$Sample)
Subtype <- data.frame(Subtype=c(rep("N",length(normal_id)),sub$Cluster),
                      id=c(normal_id,rownames(sub)))
rownames(Subtype) <- Subtype$id
Subtype <- Subtype[-2]
mygene_data <- gsva_matrix_BD[,rownames(Subtype)]
#Subtype <- Subtype[com_sam,,drop = F]
head(Subtype)
table(Subtype$Subtype)
# C1  C2   N 
# 231 285  72  
## 用前面的自定义函数计算组间统计差异
comprTab <- cross_subtype_compr(expr = mygene_data, # 或log2(mygene_data + 1)，如果用参数检验，请注意对数转化；若非参均可
                                subt = Subtype,
                                #two_sam_compr_method = "wilcox", # 两组"t.test", "wilcox"
                                multi_sam_compr_method = "kruskal", # 多组"anova", "kruskal"
                                res.path = ".")


# 用全部基因来画
n.show_top_gene <- nrow(mygene_data)

# 按分组排序
subt.order <- Subtype[order(Subtype$Subtype),,drop = F]
indata <- mygene_data[comprTab$gene[1:n.show_top_gene],rownames(subt.order)]



# 数据标准化和边界设置
plotdata <- t(scale(t(indata)))
plotdata[plotdata > 2] <- 2
plotdata[plotdata < -2] <- -2

# 调整行名
blank <- "    " # 行名和p值之间的间隔
p.value <- comprTab$adjusted.p.value[1:n.show_top_gene]
sig.label <- ifelse(p.value < 0.001,"****",
                    ifelse(p.value < 0.005,"***",
                           ifelse(p.value < 0.01,"**",
                                  ifelse(p.value < 0.05,"*",""))))
p.label <- formatC(p.value, # 将p值变成保留两位小数的科学计数法
                   format = "e",
                   digits = 2)

add.label <- str_pad(paste0(rownames(plotdata),sig.label), # 固定行名宽度并再右侧补齐" "
                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                     side = "right")

annCol <- subt.order # 获得排序后的亚型注释信息，这里只有一个变量需要注释
colnames(annCol)[1] <- paste(str_pad(colnames(annCol)[1], # 注释列名补上"P-value"，宽度和刚才一致
                                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                                     side = "right"),
                             "P-value",
                             sep = blank)

annColors <- list(c( "C1"="#1F78B4FF", "C2"="#E31A1CFF","N"="#33A02CFF")) # 如果有多个变量要注释颜色请补充c()
names(annColors) <- colnames(annCol)[1] # 如果有多个变量要注释颜色请补充每张list的name

# 绘制热图
pheatmap(cellheight = 15, cellwidth = 1,
         mat = plotdata, # 输入数据
         scale = "none", # 不标准化因为数据已经被标准化
         annotation_col = annCol, # 列注释信息
         annotation_colors = annColors, # 列注释对应的颜色
         cluster_cols = F, # 列不聚类
         cluster_rows = F, # 行不聚类
         show_colnames = F, # 不显示列名
         show_rownames = T, # 显示基因名
         annotation_legend = F, # 不显示图例
         #labels_row = paste(add.label, p.label, sep=blank), # 自定义样本名义blank作间隔
         #fontfamily = "mono", # 关键，使用fixed font而不是proportional font
         gaps_col = c(231,516), # 根据自己的数据设置空白间隔的位置
         filename = "heatmapPvalue.pdf")
dev.off()
dev.new()
library(paletteer)
pheatmap(cellheight = 15, cellwidth = 1,
         mat = plotdata, # 输入数据
         scale = "none", # 不标准化因为数据已经被标准化
         annotation_col = annCol, # 列注释信息
         annotation_colors = annColors, # 列注释对应的颜色
         cluster_cols = F, # 列不聚类
         cluster_rows = F, # 行不聚类
         show_colnames = F, # 不显示列名
         show_rownames = T, # 显示基因名
         annotation_legend = F, # 不显示图例
         gaps_col = c(231,516),
         color = paletteer_c("scico::berlin", n = 100),
         labels_row = paste(add.label, p.label, sep=blank), # 自定义样本名义blank作间隔
         #fontfamily = "Times New Roman"
         #filename = "heatmapPvalue.pdf",
         height = 5,
         width = 15)

pheatmap(cellheight = 10, cellwidth = 1,
         mat = plotdata, # 输入数据
         scale = "none", # 不标准化因为数据已经被标准化
         annotation_col = annCol, # 列注释信息
         annotation_colors = annColors, # 列注释对应的颜色
         cluster_cols = F, # 列不聚类
         cluster_rows = F, # 行不聚类
         show_colnames = F, # 不显示列名
         show_rownames = T, # 显示基因名
         annotation_legend = F, # 不显示图例
         gaps_col = c(231,516),
         color = paletteer_c("scico::berlin", n = 100),
         labels_row = paste(add.label, p.label, sep=blank),
)
table(Subtype$Subtype)

###差异分析####
setwd("/home/data/vip39/TCGA/cancer_stemness/")
colnames(KIRC_mRNA_count) <- str_sub(colnames(KIRC_mRNA_count),1,15)
exp <-KIRC_mRNA_count[,sub$Sample]

group_list <- sub$Cluster
#C1预后差,放后面
group_list = factor(group_list,levels = c("C2","C1"))
library(DESeq2)
colData <- data.frame(row.names =colnames(exp), 
                      condition=group_list)
exp <- round(exp)
range(exp)
exp = exp[apply(exp, 1, function(x) sum(x > 1) > 5), ]
if(!file.exists(paste0("CS2_CS1","dd.Rdata"))){
  dds <- DESeqDataSetFromMatrix(
    countData = exp,
    colData = colData,
    design = ~ condition)
  dds <- DESeq(dds)
  save(dds,file = paste0("CS1_CS2","_dd.Rdata"))
}
load(paste0("CS1_CS2","_dd.Rdata"))

# 两两比较
res <- results(dds, contrast = c("condition",rev(levels(group_list))))
resOrdered <- res[order(res$pvalue),] # 按照P值排序
DEG <- as.data.frame(resOrdered)
head(DEG)
table(is.na(DEG))

#添加change列标记基因上调下调
logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )

logFC_cutoff <- 2
DEG$change = as.factor(
  ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
         ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
)
head(DEG)
table(DEG$change)
# DOWN   NOT    UP 
# 110 18295   763
library(EnhancedVolcano,lib.loc = "/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.0/")
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'CS1 versus CS2',
                pCutoff = 0.05,
                FCcutoff = logFC_cutoff,
                pointSize = 2.0,
                labSize = 3.0,
                col=c('#b5b5b5','#b5b5b5','#4D4398','#F18D00'),#调整颜色
                colAlpha=0.6)

##火山图美化
EnhancedVolcano(res,
                x='log2FoldChange',
                y='pvalue',
                #lab=data$id,
                title = 'CS2 versus CS1',
                pCutoff=10e-1/20,#y轴阈值线(水平)
                FCcutoff=1.5,#x轴阈值线（垂直）
                pointSize=3,#散点大小
                labSize=4,#标签大小
                xlim=c(-6, 6),#限制X轴范围
                ylim=c(0,300),#限制Y轴范围
                col=c('#b5b5b5','#b5b5b5','#4D4398','#F18D00'),#调整颜色
                colAlpha=0.6,#调整透明度
                #selectLab=c(downvals,upvals),#使用selectLab参数选定所关注的标签
                xlab=bquote(~Log[2]~'fold change'))#将内容传递给xlab
####↓新加入↓####
#labCol='black',#标签颜色
#labFace='bold',#标签字体
#boxedLabels=TRUE,#是否在框中绘制标签
#drawConnectors=TRUE,#是否通过连线将标签连接到对应的点上
#widthConnectors=0.8,#连线的宽度
#endsConnectors="last",#连线绘制箭头的方向，可选first、both、last
#colConnectors='black')#连线的颜色




topptx(file="CS2_CS1_EnhancedVolcano.pptx")

######富集分析#####
library(clusterProfiler)
library(AnnotationHub)
library(AnnotationDbi)
#BiocManager::install("GOplot")
library(GOplot)
library(ggplot2,lib.loc = "/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.0/")
colnames(DEG)
library(ReactomePA,lib.loc = "/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.0/")
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(ggsci)
library(paletteer)
DEG$SYMBOL <- rownames(DEG)
DEG2 <-DEG%>%
  arrange(log2FoldChange)
top150down <- DEG2$SYMBOL[1:150]
top150up <- DEG2$SYMBOL[(nrow(DEG2)-149):nrow(DEG2)]
#write.table(top150up,"up_CS2_map.txt",col.names = F,row.names = F,sep = "\t",quote = F)

#write.table(top150down,"down_CS2_map.txt",col.names = F,row.names = F,sep = "\t",quote = F)

genelist_input <- DEG[,c(2,8)]
genename <- as.character(genelist_input[,2])
class(genename)
gene_map<-bitr(genelist_input$SYMBOL, #转换的列是df数据框中的SYMBOL列
               fromType = "SYMBOL",#需要转换ID类型
               toType = "ENTREZID",#转换成的ID类型
               OrgDb = "org.Hs.eg.db")#对应的物种，小鼠的是org.Mm.eg.db
gene_map <- select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
#gene_map <- select(org.Hs.eg.db, keys=genelist_input[,2], keytype="SYMBOL", columns=c("ENTREZID"))
head(gene_map)
colnames(gene_map)[1]<-"Gene"
colnames(genelist_input) <- c("logFC","Gene")
class(genelist_input$Gene)
aaa<-inner_join(gene_map,genelist_input,by = "Gene")
head(aaa)
aaa<-aaa[,-1]
aaa<-na.omit(aaa)
aaa$logFC<-sort(aaa$logFC,decreasing = T)
geneList = aaa[,2]
names(geneList) = as.character(aaa[,1])
geneList
#GSEA分析——GO
Go_gseresult <- gseGO(geneList, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont="all",  minGSSize = 10, maxGSSize = 500, pvalueCutoff=1)
#GSEA分析——KEGG
KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.5)
#GSEA分析——Reactome
Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#波浪图
ridgeplot(KEGG_gseresult,10) #输出前十个结果
ridgeplot(Go_gseresult,10) #输出前十个结果

gseaplot2(Go_Reactomeresult, geneSetID = 1:3, pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")


gseaplot2(Go_Reactomeresult, geneSetID = 1:5, pvalue_table = TRUE,
          color = paletteer_c("scico::berlin", n = 5), ES_geom = "dot")
topptx(filename = "GSEA_up.pptx")
# go <- enrichGO(gene = geneList, OrgDb = "org.Hs.eg.db", ont="all")
# barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")

gseaplot2(Go_Reactomeresult, geneSetID = c("R-HSA-381753","R-HSA-418555","R-HSA-1461957",
                                           "R-HSA-8851680","R-HSA-176974"), pvalue_table = TRUE,
          color = paletteer_c("scico::berlin", n = 5), ES_geom = "dot")


####
gsym.id <- DEG[,c("SYMBOL","log2FoldChange","change")]
head(gsym.id)
gsym.id <- gsym.id[gsym.id$change!="NOT",]

gsym <- bitr(gsym.id$SYMBOL, #基因名
             fromType = "SYMBOL", #从gene symbol
             toType = "ENTREZID", #提取ENTREZ ID
             OrgDb = "org.Hs.eg.db") #相应物种的包，人是org.Hs.eg.db


head(gsym)
gsym.id <- inner_join(gsym.id,gsym,by="SYMBOL")
head(gsym.id)
dim(gsym.id)
id.fc <- gsym.id[,c(2,4)]
#save(id.fc,file="very_easy_input_H.Rdata")

ego <- enrichGO(gene = id.fc$ENTREZID,
                #小鼠用这行
                OrgDb = org.Hs.eg.db,
                #人类用这行
                #OrgDb = org.Hs.eg.db,
                #非模式生物用这行，例如玉米
                #OrgDb = maize.db,
                ont = "BP", #或MF或CC
                pAdjustMethod = "BH",
                #pvalueCutoff  = 0.001,
                qvalueCutoff  = 0.05) 
dim(ego)
#write.csv(ego,"enrichGO_output.csv",quote = F)
#把ENTREZ ID转为gene symbol
egox <- setReadable(ego, 'org.Hs.eg.db', #物种
                    'ENTREZID')

#把基因倍数信息转成画图所需的格式
geneList <- id.fc$log2FoldChange
names(geneList) <- id.fc$ENTREZID

#画⭕️图
#BiocManager::install("ggnewscale")
library(ggnewscale,lib.loc = "/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.0/")
cnetplot(egox, 
         foldChange = geneList, 
         #foldChange = NULL, #不展示倍数
         circular = TRUE,
         #node_label = FALSE, #如果太多，就不要显示基因名了
         showCategory = 5, #显示富集的term数量，默认5
         colorEdge = TRUE)

#保存到pdf文件
#ggsave("clusterProfiler_circle.pdf", width = 8, height = 5)

#不画成⭕️，process分开，效果更好呢
cnetplot(egox, 
         foldChange = geneList, 
         #foldChange = NULL, #不展示倍数
         #circular = TRUE,
         #node_label = FALSE, #不显示基因名
         showCategory = 4, #显示的富集term数量，默认5
         colorEdge = TRUE)

#ggsave("clusterProfiler_not_circle.pdf", width = 8, height = 5)

#读入富集分析结果
#ego <- read.csv("enrichGO_output.csv", header = T)
ego[1,]
go <- data.frame(Category = "BP",
                 ID = ego$ID,
                 Term = ego$Description, 
                 Genes = gsub("/", ", ", ego$geneID), 
                 adj_pval = ego$p.adjust)

#基因变化倍数

head(id.fc)
genelist <- data.frame(ID = id.fc$ENTREZID, logFC = id.fc$log2FoldChange)             

#把富集分析和倍数整合在一起
library(GOplot,lib.loc = "/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.0/")
circ <- circle_dat(go, genelist)
head(circ)

#根据ENTREZ ID提取gene symbol
id.gsym <- bitr(circ$genes, #基因名
                fromType = "ENTREZID", #从ENTREZ ID
                toType = "SYMBOL", #提取gene symbol
                OrgDb = "org.Hs.eg.db") #相应物种的包

#把circ里的ENTREZ ID换成gene symbol
rownames(id.gsym) <- id.gsym$ENTREZID
circ.gsym <- circ
circ.gsym$genes <- id.gsym[circ$genes,]$SYMBOL
head(circ.gsym)
#参数设置
n = 5 #圈图需要选定term，这里画前面5个

chord <- chord_dat(circ, genelist, go$Term[1:n])
head(chord)

#根据ENTREZ ID提取gene symbol
id.gsym <- bitr(row.names(chord), #基因名
                fromType = "ENTREZID", #从ENTREZ ID
                toType = "SYMBOL", #提取gene symbol
                OrgDb = "org.Hs.eg.db") #相应物种的包

#把chord的列名ENTREZ ID换成gene symbol
rownames(id.gsym) <- id.gsym$ENTREZID
head(id.gsym)
chord.gsym <- chord
row.names(chord.gsym) <- id.gsym[row.names(chord),]$SYMBOL
head(chord.gsym)

GOChord(chord.gsym, 
        space = 0.02, #基因方块间隙
        gene.order = 'logFC', 
        lfc.col = c('darkgoldenrod1', 'black', 'cyan1'), #自定义变化倍数的颜色
        gene.space = 0.25, #基因名跟⭕️的相对距离
        gene.size = 3, #基因名字体大小 
        border.size = 0.1, #中间曲线的黑色边的粗细
        process.label = 8) #term字体大小

ggsave("GOChordBP.pdf", width = 12, height = 14)
ggsave("GOChordCC.pdf", width = 12, height = 14)
ggsave("GOChordMF.pdf", width = 12, height = 14)
graph2ppt(file="GOChord_CS2_CS1.pptx")

#定义足够多的颜色，后面从这里选颜色
mycol <- c("#223D6C","#D20A13","#FFD121","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

GOCluster(circ.gsym, go$Term[1:n], 
          #clust.by = 'logFC', #用变化倍数聚类
          clust.by = 'term', #用富集的term聚类
          lfc.col = c('darkgoldenrod1', 'black', 'cyan1'), #自定义变化倍数的颜色
          lfc.space = 0.1, #倍数跟树间的空隙大小
          lfc.width = 1, #变化倍数的⭕️宽度
          term.col = mycol[1:n], #自定义term的颜色
          term.space = 0.2, #倍数跟term间的空隙大小
          term.width = 0.8) #富集term的⭕️宽度
ggsave("GOCluster.pdf", width = 12, height = 14)
graph2ppt(file="GOCluster_CS2_CS1_BP.pptx")
topptx(filename = "GOCluster_CS2_CS1_BP.pptx")
if(!require(paletteer))install.packages("paletteer")
if(!require(scico))install.packages('scico')
if(!require(nord))install.packages('nord')
library(paletteer)
##CC
ego_CC <- enrichGO(gene = id.fc$ENTREZID,
                   #小鼠用这行
                   OrgDb = org.Hs.eg.db,
                   #人类用这行
                   #OrgDb = org.Hs.eg.db,
                   #非模式生物用这行，例如玉米
                   #OrgDb = maize.db,
                   ont = "CC", #或MF或CC
                   pAdjustMethod = "BH",
                   #pvalueCutoff  = 0.001,
                   qvalueCutoff  = 0.05) 
dim(ego_CC)
#write.csv(ego_CC,"enrichGOCC_output.csv",quote = F)
ego_CC[1,]
goCC <- data.frame(Category = "CC",
                   ID = ego_CC$ID,
                   Term = ego_CC$Description, 
                   Genes = gsub("/", ", ", ego_CC$geneID), 
                   adj_pval = ego_CC$p.adjust)

#基因变化倍数

head(id.fc)
genelist <- data.frame(ID = id.fc$ENTREZID, logFC = id.fc$log2FoldChange)             

#把富集分析和倍数整合在一起
circ <- circle_dat(goCC, genelist)
head(circ)

#根据ENTREZ ID提取gene symbol
id.gsym <- bitr(circ$genes, #基因名
                fromType = "ENTREZID", #从ENTREZ ID
                toType = "SYMBOL", #提取gene symbol
                OrgDb = "org.Hs.eg.db") #相应物种的包

#把circ里的ENTREZ ID换成gene symbol
rownames(id.gsym) <- id.gsym$ENTREZID
circ.gsym <- circ
circ.gsym$genes <- id.gsym[circ$genes,]$SYMBOL
head(circ.gsym)

#参数设置
n = 5 #圈图需要选定term，这里画前面5个

chord <- chord_dat(circ, genelist, goCC$Term[1:n])
head(chord)

#根据ENTREZ ID提取gene symbol
id.gsym <- bitr(row.names(chord), #基因名
                fromType = "ENTREZID", #从ENTREZ ID
                toType = "SYMBOL", #提取gene symbol
                OrgDb = "org.Hs.eg.db") #相应物种的包

#把chord的列名ENTREZ ID换成gene symbol
rownames(id.gsym) <- id.gsym$ENTREZID
head(id.gsym)
chord.gsym <- chord
row.names(chord.gsym) <- id.gsym[row.names(chord),]$SYMBOL
head(chord.gsym)

#定义足够多的颜色，后面从这里选颜色
mycol <- c("#223D6C","#D20A13","#FFD121","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

GOCluster(circ.gsym, goCC$Term[1:n], 
          #clust.by = 'logFC', #用变化倍数聚类
          clust.by = 'term', #用富集的term聚类
          lfc.col = c('darkgoldenrod1', 'black', 'cyan1'), #自定义变化倍数的颜色
          lfc.space = 0.1, #倍数跟树间的空隙大小
          lfc.width = 1, #变化倍数的⭕️宽度
          term.col = paletteer_c("scico::berlin", n = 5), #自定义term的颜色
          term.space = 0.2, #倍数跟term间的空隙大小
          term.width = 0.8) #富集term的⭕️宽度

ggsave("GOClusterCC.pdf", width = 12, height = 14)
graph2ppt(file="GOCluster_CS2_CS1_CC.pptx")
topptx(filename = "GOCluster_CS2_CS1_CC.pptx")
##MF
ego_MF <- enrichGO(gene = id.fc$ENTREZID,
                   #小鼠用这行
                   OrgDb = org.Hs.eg.db,
                   #人类用这行
                   #OrgDb = org.Hs.eg.db,
                   #非模式生物用这行，例如玉米
                   #OrgDb = maize.db,
                   ont = "MF", #或MF或MF
                   pAdjustMethod = "BH",
                   #pvalueCutoff  = 0.001,
                   qvalueCutoff  = 0.05) 
dim(ego_MF)
#write.csv(ego_MF,"enrichGOMF_output.csv",quote = F)
ego_MF[1,]
goMF <- data.frame(Category = "MF",
                   ID = ego_MF$ID,
                   Term = ego_MF$Description, 
                   Genes = gsub("/", ", ", ego_MF$geneID), 
                   adj_pval = ego_MF$p.adjust)

#基因变化倍数

head(id.fc)
genelist <- data.frame(ID = id.fc$ENTREZID, logFC = id.fc$log2FoldChange)             

#把富集分析和倍数整合在一起
circ <- circle_dat(goMF, genelist)
head(circ)

#根据ENTREZ ID提取gene symbol
id.gsym <- bitr(circ$genes, #基因名
                fromType = "ENTREZID", #从ENTREZ ID
                toType = "SYMBOL", #提取gene symbol
                OrgDb = "org.Hs.eg.db") #相应物种的包

#把circ里的ENTREZ ID换成gene symbol
rownames(id.gsym) <- id.gsym$ENTREZID
circ.gsym <- circ
circ.gsym$genes <- id.gsym[circ$genes,]$SYMBOL
head(circ.gsym)

#参数设置
n = 5 #圈图需要选定term，这里画前面5个

chord <- chord_dat(circ, genelist, goMF$Term[1:n])
head(chord)

#根据ENTREZ ID提取gene symbol
id.gsym <- bitr(row.names(chord), #基因名
                fromType = "ENTREZID", #从ENTREZ ID
                toType = "SYMBOL", #提取gene symbol
                OrgDb = "org.Hs.eg.db") #相应物种的包

#把chord的列名ENTREZ ID换成gene symbol
rownames(id.gsym) <- id.gsym$ENTREZID
head(id.gsym)
chord.gsym <- chord
row.names(chord.gsym) <- id.gsym[row.names(chord),]$SYMBOL
head(chord.gsym)

#定义足够多的颜色，后面从这里选颜色
mycol <- c("#223D6C","#D20A13","#FFD121","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7MF767")

GOCluster(circ.gsym, goMF$Term[1:n], 
          #clust.by = 'logFC', #用变化倍数聚类
          clust.by = 'term', #用富集的term聚类
          lfc.col = c('darkgoldenrod1', 'black', 'cyan1'), #自定义变化倍数的颜色
          lfc.space = 0.1, #倍数跟树间的空隙大小
          lfc.width = 1, #变化倍数的⭕️宽度
          term.col = paletteer_d("RColorBrewer::Paired",n=5), #自定义term的颜色
          term.space = 0.2, #倍数跟term间的空隙大小
          term.width = 0.8) #富集term的⭕️宽度
ggsave("GOClusterMF.pdf", width = 12, height = 14)
graph2ppt(file="GOCluster_C2_C1_MF.pptx")
topptx(filename = "GOCluster_C2_C1_MF.pptx")






####gsva#####
library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# 读入前面生成的通路中的基因列表
(load("~/huitu/FigureYa61GSVA/hallmark.gs.RData")) #保存在当前文件夹

# 读入基因表达矩阵
gsym.expr <- KIRC_mRNA_fpkm
head(gsym.expr)

# 这一句就完成了GSVA分析
gsva_es <- gsva(as.matrix(gsym.expr), gs)

# 预览GSVA分析返回的矩阵
head(gsva_es)

# 把通路的表达量保存到文件
write.csv(gsva_es, "gsva_output.csv", quote = F)



# 分组
gsva_es <- gsva_es[,sub$Sample]
group_list <- data.frame(sample = sub$Sample, 
                         group = sub$Cluster)
head(group_list)

# 设置对比
design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_es)
design

# 构建差异比较矩阵
contrast.matrix <- makeContrasts(C1-C2, levels = design)

# 差异分析，b vs. a
fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)

#把通路的limma分析结果保存到文件
#write.csv(x, "gsva_limma.csv", quote = F)

#输出t值，用做FigureYa39bar的输入数据
pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t)
#write.csv(df, "easy_input2_for39bar.csv", quote = F, row.names = F)


head(df)

#按照score的值分组
cutoff <- 1
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))

#按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('palegreen3', 'snow3', 'dodgerblue4'), guide = FALSE) + 
  
  #画2条虚线
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, #画虚线
             size = 0.3) + #线的粗细
  
  #写label
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = "inward" ) +  #字的对齐方式
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1, label=ID, color = group),
            size = 3, hjust = "outward") +  
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score, CS1 versus CS2")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴
library(ggplot2)

library(ggtheme,lib.loc = "/home/data/refdir/Rlib")
#install.packages("ggtheme")
library(ggprism)
ggplot(data = sortdf,aes(x = ID,y = score,fill = group)) +
  geom_col()+
  coord_flip() +
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, CS2 versus CS1') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave("gsva.pdf", width = 6, height = 8)
graph2ppt(file="gsva.pdf.pptx")



head(sortdf)
ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('palegreen3', 'snow3', 'dodgerblue4'), guide = FALSE) + 
  
  #画2条虚线
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, #画虚线
             size = 0.3) + #线的粗细
  
  #写label
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = "outward" ) +  #字的对齐方式
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1, label=ID, color = group),
            size = 3, hjust = "inward") +  
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score, C1 versus C2")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴
ggsave("gsva.pdf", width = 6, height = 8)



library(ggplot2)
library(BiocManager)
#BiocManager::install("enrichR")
library(enrichR,lib.loc = "/home/data/vip39/R/x86_64-pc-linux-gnu-library/4.0/")
listEnrichrSites()
#链接到human数据上
setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
#调取需要的GO或者KEGG
dbs <- c( "GO_Biological_Process_2018",
          "GO_Cellular_Component_2018",
          "GO_Molecular_Function_2018",
          "KEGG_2016")
##读取我们的差异基因文件，因为可对这些差异基因做GO和KEGG分析
# library(readxl)
# #提前设置好logFC2_FR35.xlsx所在的路径
# setwd("~/Documents/xmu/ZhaoBin")
# logGene <- read_excel("logFC2_FR35.xlsx", sheet = 1, col_names = T)
# logGene1<-logGene$Gene
list_up <-DEG[DEG$change=="UP","SYMBOL"]
list_down<-DEG[DEG$change=="DOWN","SYMBOL"]
eup <- enrichr(list_up, dbs)
edown <- enrichr(list_down, dbs)
###提取KEGG信息
up <- eup$KEGG_2016
down <- edown$KEGG_2016
#添加一个type列变量
up$type <-"up"
down$type <-"down"
#提取前15个KEGG
up <- up[c(1:15),]
up <- up[order(up$Combined.Score), ]  # sort
down <- down [c(1:15),]
down <- down[order(down$Combined.Score), ]  # sort
keggmerges <- rbind(down,up)
keggmerges <- keggmerges[!duplicated(keggmerges$Term),]
keggmerges$Term <- factor(keggmerges$Term, levels=keggmerges$Term)
#绘制柱状图
ggplot(keggmerges, aes(x=Term, y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=type), width=.5,position="dodge")  +
  scale_fill_manual(name="Expression", 
                    labels = c("Down regulated", "Up regulated"), 
                    values = c("down"="#1F78B4FF", "up"="#E31A1CFF")) + 
  labs(subtitle="Combined scores from Kegg pathways", 
       title= "DEGs KEGG Pathway") + 
  coord_flip()+
  theme_bw()
ggsave("KEGG.pdf", width = 10, height = 6)
topptx(filename = "kegg.pptx")




#####两组间热图####
head(Subtype)
head(mygene_data)
Subtype2 <- data.frame(Subtype=sub$Cluster)
rownames(Subtype2) <- rownames(sub)
mygene_data2 <- mygene_data[,rownames(sub)]
## 用前面的自定义函数计算组间统计差异
comprTab <- cross_subtype_compr(expr = mygene_data, # 或log2(mygene_data + 1)，如果用参数检验，请注意对数转化；若非参均可
                                subt = Subtype,
                                #two_sam_compr_method = "wilcox", # 两组"t.test", "wilcox"
                                multi_sam_compr_method = "kruskal", # 多组"anova", "kruskal"
                                res.path = ".")
# 用全部基因来画
n.show_top_gene <- nrow(mygene_data)
# 或者取top 20个基因来画
#n.show_top_gene <- 20 

# 按分组排序
subt.order <- Subtype[order(Subtype$Subtype),,drop = F]
indata <- mygene_data[comprTab$gene[1:n.show_top_gene],rownames(subt.order)]
```

# 开始画图

```{r}
# 数据标准化和边界设置
plotdata <- t(scale(t(indata)))
plotdata[plotdata > 2] <- 2
plotdata[plotdata < -2] <- -2

# 调整行名
blank <- "    " # 行名和p值之间的间隔
p.value <- comprTab$adjusted.p.value[1:n.show_top_gene]
sig.label <- ifelse(p.value < 0.001,"****",
                    ifelse(p.value < 0.005,"***",
                           ifelse(p.value < 0.01,"**",
                                  ifelse(p.value < 0.05,"*",""))))
p.label <- formatC(p.value, # 将p值变成保留两位小数的科学计数法
                   format = "e",
                   digits = 2)

add.label <- str_pad(paste0(rownames(plotdata),sig.label), # 固定行名宽度并再右侧补齐" "
                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                     side = "right")

annCol <- subt.order # 获得排序后的亚型注释信息，这里只有一个变量需要注释
colnames(annCol)[1] <- paste(str_pad(colnames(annCol)[1], # 注释列名补上"P-value"，宽度和刚才一致
                                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                                     side = "right"),
                             "P-value",
                             sep = blank)

annColors <- list(c("cancer"="lightblue","normal"="pink")) # 如果有多个变量要注释颜色请补充c()
names(annColors) <- colnames(annCol)[1] # 如果有多个变量要注释颜色请补充每张list的name

# 绘制热图
p <- pheatmap(cellheight = 15, cellwidth = 1,
              mat = plotdata, # 输入数据
              scale = "none", # 不标准化因为数据已经被标准化
              annotation_col = annCol, # 列注释信息
              annotation_colors = annColors, # 列注释对应的颜色
              cluster_cols = F, # 列不聚类
              cluster_rows = T, # 行不聚类
              show_colnames = F, # 不显示列名
              show_rownames = T, # 显示基因名
              #annotation_legend = F, # 不显示图例
              labels_row = paste(add.label, p.label, sep=blank) # 自定义样本名义blank作间隔
              #fontfamily = "mono", # 关键，使用fixed font而不是proportional font
              #gaps_col = c(531), # 根据自己的数据设置空白间隔的位置
              # filename = "heatmapPvalue.pdf"
)


#####两组间PCA图####
dim(sub)
head(sub)
library(ggplot2)
library(plyr)
library(ggord)
source('~/database/geom_ord_ellipse.R') #该文件位于当前文件夹

#用`prcomp`进行PCA分析
mygene_data <- gsva_matrix_BD[,sub$Sample]
dim(mygene_data)
mygene_data2 <- mygene_data[,rownames(sub)]
Subtype2 <- data.frame(Subtype=sub$Cluster)
rownames(Subtype2) <- rownames(sub)
pca.results <- prcomp(t(mygene_data2), center = TRUE, scale. = FALSE)

#定义足够多的颜色，用于展示分组
mycol <- c("#223D6C","#D20A13","#088247","#FFD121","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

ggord(pca.results, grp_in = Subtype2$Subtype, repel=TRUE,
      ellipse = FALSE, #不显示置信区间背景色
      size = 2, #样本的点大小
      alpha=0.5, #设置点为半透明，出现叠加的效果
      #如果用自定义的颜色，就运行下面这行
      cols = mycol[1:length(unique(Subtype2$Subtype))],
      arrow = NULL,txt = NULL) + #不画箭头和箭头上的文字
  theme(panel.grid =element_blank()) + #去除网格线
  
  #用yyplot添加置信区间圆圈
  geom_ord_ellipse(ellipse_pro = .95, #设置置信区间
                   size=1.5, #线的粗细
                   lty=1 ) #实线

#保存到pdf文件
ggsave("PCA_classic.pdf", width = 6, height = 6)
topptx(filename = "PCA.pptx",width = 6,height = 6)




####免疫部分#####
##添加纵坐标
library(RColorBrewer)
library(circlize)
library(gplots)
library(viridis)
library(oompaBase)
library(data.table)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

hmdat <- read.csv("~/vip39/huitu/FigureYa230immunelandscape/easy_input.csv",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
type <- read.csv("~/vip39/huitu/FigureYa230immunelandscape/easy_input_type.csv", row.names = 1)
head(type)

table(rownames(sub)%in%rownames(hmdat))

comsam <- rownames(sub)
hmdat <- hmdat[comsam,]
risk <- sub[comsam,,drop = F]
dim(hmdat)


# 拆分不同算法结果，获得类的名字
#immMethod <- sapply(strsplit(colnames(hmdat),"_",fixed = T),"[",2) #用easy_input.csv列名里的算法信息
immMethod <- type$Methods # 用easy_input_type.csv的算法那一列




library(ComplexHeatmap) 
# 最新版ComplexHeatmap好像有一些bug，使用里面的pheatmap函数会报错
# 因此我们从脚本直接加载pheatmap函数
source("~/vip39/huitu/FigureYa230immunelandscape/pheatmap_translate.R") # 位于当前文件夹，出自ComplexHeatmap_2.7.9.tar.gz

ht_opt$message = FALSE


methods.col <- brewer.pal(n = length(unique(immMethod)),name = "Paired")

# 创建注释
annCol <- data.frame(
  Subtype = risk$Cluster,
  row.names = rownames(risk),
  stringsAsFactors = F)
annRow <- data.frame(row.names = colnames(hmdat),
                     Methods = factor(immMethod,levels = unique(immMethod)),
                     stringsAsFactors = F)

# 为各注释信息设置颜色
# annColors <- list(Methods = c("TIMER" = methods.col[1], #行注释的颜色
#                               "CIBERSORT" = methods.col[2],
#                               "CIBERSORT-ABS" = methods.col[3],
#                               "QUANTISEQ" = methods.col[4],
#                               "MCPCOUNTER" = methods.col[5],
#                               "XCELL" = methods.col[6],
#                               "EPIC" = methods.col[7]),
#                   # 下面是列注释的颜色，可依此设置更多注释的颜色
#                   #"RiskScore" = greenred(64), 
#                   "RiskType" = c("high" = "red","low" = "blue"))

annColors <- list(
  "Subtype" = c("C1"="#1F78B4FF", "C2"="#E31A1CFF"))

# 数据标准化
indata <- t(hmdat)
indata <- indata[,colSums(indata) > 0] # 确保没有富集全为0的细胞
plotdata <- standarize.fun(indata,halfwidth = 2)

# 样本按risk score排序
samorder <- c(rownames(sub)[sub$Cluster=="C1"],rownames(sub)[sub$Cluster=="C2"])

# 拆分各算法的结果
plotdata1 <- plotdata[rownames(annRow[which(annRow$Methods == "TIMER"),,drop = F]),]
plotdata2 <- plotdata[rownames(annRow[which(annRow$Methods == "CIBERSORT"),,drop = F]),]
plotdata3 <- plotdata[rownames(annRow[which(annRow$Methods == "CIBERSORT-ABS"),,drop = F]),]
plotdata4 <- plotdata[rownames(annRow[which(annRow$Methods == "QUANTISEQ"),,drop = F]),]
plotdata5 <- plotdata[rownames(annRow[which(annRow$Methods == "MCPCOUNTER"),,drop = F]),]
plotdata6 <- plotdata[rownames(annRow[which(annRow$Methods == "XCELL"),,drop = F]),]
plotdata7 <- plotdata[rownames(annRow[which(annRow$Methods == "EPIC"),,drop = F]),]

# 分别画7次热图（参数基本同pheatmap里的pheatmap）
hm1 <- pheatmap(mat = as.matrix(plotdata1[,samorder]),
                border_color = NA,
                #color = bluered(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                annotation_col = annCol[samorder,,drop = F],
                annotation_colors = annColors,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$Subtype)[1],
                name = "TIMER") # 为子热图的图例命名

hm2 <- pheatmap(mat = as.matrix(plotdata2[,samorder]),
                border_color = NA,
                color = greenred(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$Subtype)[1],
                name = "CIBERSORT")

hm3 <- pheatmap(mat = as.matrix(plotdata3[,samorder]),
                border_color = NA,
                color = blueyellow(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$Subtype)[1],
                name = "CIBERSORT-ABS")

hm4 <- pheatmap(mat = as.matrix(plotdata4[,samorder]),
                border_color = NA,
                color = bluered(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$Subtype)[1],
                name = "QUANTISEQ")

hm5 <- pheatmap(mat = as.matrix(plotdata5[,samorder]),
                border_color = NA,
                color = inferno(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$Subtype)[1],
                name = "MCPCOUNTER")

hm6 <- pheatmap(mat = as.matrix(plotdata6[,samorder]),
                border_color = NA,
                color = viridis(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$Subtype)[1],
                name = "XCELL")

hm7 <- pheatmap(mat = as.matrix(plotdata7[,samorder]),
                border_color = NA,
                color = magma(64), 
                cluster_rows = F,
                cluster_cols = F, 
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$Subtype)[1],
                name = "EPIC")

pdf("immune heatmap by ComplexHeatmap met.pdf", width = 10,height = 20) # 保存前请注意RGUI里不能有任何显示的图像，否则不会pdf打不开
draw(hm1 %v% hm2 %v% hm3 %v% hm4 %v% hm5 %v% hm6 %v% hm7, # 垂直连接子热图
     heatmap_legend_side = "bottom", # 热图注释放底部
     annotation_legend_side = "bottom") # 顶部注释放底部
invisible(dev.off())

#
save(sub,file="cancerstem_sub.rds")
load("./cancerstem_sub.rds")



##单细胞部分 DAVID数据  不同分期上皮细胞#####


###药敏####


###突变谱#####



###模型可作可不做####


###单基因泛癌分析####


###肾癌治疗队列####


#####干性分数  MSI TMB estimate 等指标######
RNA_stemness <- fread("/home/data/vip39/data/pancancer/StemnessScores_RNAexp_20170127.2.tsv.gz",data.table = F)
RNA_stemness <- as.data.frame(t(RNA_stemness))
head(RNA_stemness)
colnames(RNA_stemness) <- RNA_stemness[1,]
RNA_stemness <- RNA_stemness[-1,]
head(RNA_stemness)
rownames(RNA_stemness) <- str_replace_all(rownames(RNA_stemness),"\\.","\\-")
table(rownames(sub)%in%rownames(RNA_stemness))

RNA_stemness$Sample <- rownames(RNA_stemness)

subdata <- merge(sub,RNA_stemness,by="Sample")
head(subdata)
subdata$Cluster <- factor(subdata$Cluster,levels = c("C1","C2"))

DNA_stemness <- fread("/home/data/vip39/data/pancancer/StemnessScores_DNAmeth_20170210.tsv.gz",data.table = F)
DNA_stemness <- as.data.frame(t(DNA_stemness))
head(DNA_stemness)
colnames(DNA_stemness) <- DNA_stemness[1,]
DNA_stemness <- DNA_stemness[-1,]
head(DNA_stemness)
rownames(DNA_stemness) <- str_replace_all(rownames(DNA_stemness),"\\.","\\-")
table(rownames(sub)%in%rownames(DNA_stemness))

DNA_stemness$Sample <- rownames(DNA_stemness)
subdata2 <- merge(sub,DNA_stemness,by="Sample")
head(subdata2)

library(ggpubr)
# data(ToothGrowth)
# class(ToothGrowth)
# ggboxplot(ToothGrowth,x="supp",
#           y="len",color = "supp",
#           palette = "jco",add = "jitter")
# class(subdata)

my_comparisons <- list(c("C1", "C2"), c("C2", "C3"), c("C1", "C3"))
subdata$RNAss <- as.numeric(subdata$RNAss)
p <- ggboxplot(subdata, x="Cluster",
               y = "RNAss", color = "Cluster",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
ggviolin(subdata, x = "Cluster", y = "RNAss",
         fill = "Cluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means(label = "p.signif")
#+stat_compare_means(label.y = 0.5)

subdata$RNAss <- as.numeric(subdata$RNAss)
ggviolin(subdata, x = "Cluster", y = "RNAss",
         fill = "Cluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means(label = "p.signif")
ggsave(filename = "RNAss.pdf",height = 5,width = 5)

subdata$EREG.EXPss <- as.numeric(subdata$EREG.EXPss)
ggviolin(subdata, x = "Cluster", y = "EREG.EXPss",
         fill = "Cluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means(label = "p.signif")
ggsave(filename = "EREG.EXPss.pdf",height = 5,width = 5)


subdata2$DNAss <- as.numeric(subdata2$DNAss)
subdata2$Cluster <- factor(subdata2$Cluster,levels = c("C1","C2","C3"))
ggviolin(subdata2, x = "Cluster", y = "DNAss",
         fill = "Cluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means( label = "p.signif")
ggsave(filename = "DNAss.pdf",height = 5,width = 5)
# p <- ggboxplot(subdata, x="Cluster",
#                y = "EREG.EXPss", color = "Cluster",
#                palette = "jco", add = "jitter")
# # 添加p值
# p + stat_compare_means()
# ggviolin(subdata2, x = "Cluster", y = "DNAss",
#          fill = "Cluster", palette = c("#48D1CC", "#E9967A"),alpha=.8,color = "white",add = "boxplot",
#          add.params = list(fill=c("#AFEEEE","#FFE4E1"))) + 
#   stat_compare_means()
subdata2$DMPss<- as.numeric(subdata2$DMPss)
ggviolin(subdata2, x = "Cluster", y = "DMPss",
         fill = "Cluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means(label = "p.signif")
ggsave(filename = "DMPss.pdf",height = 5,width = 5)

subdata2$ENHss<- as.numeric(subdata2$ENHss)
ggviolin(subdata2, x = "Cluster", y = "ENHss",
         fill = "Cluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means(label = "p.signif")
ggsave(filename = "ENHss.pdf",height = 5,width = 5)

HRD_score <- fread("/home/data/vip39/data/pancancer/TCGA.HRD_withSampleID.txt.gz",data.table = F)
HRD_score <- as.data.frame(t(HRD_score))
head(HRD_score)
colnames(HRD_score) <- HRD_score[1,]
HRD_score <- HRD_score[-1,]
head(HRD_score)
#rownames(DNA_stemness) <- str_replace_all(rownames(DNA_stemness),"\\.","\\-")
table(rownames(sub)%in%rownames(HRD_score))

HRD_score$Sample <- rownames(HRD_score)
subdata4 <- merge(sub,HRD_score,by="Sample")
head(subdata4)

subdata4$HRD <- as.numeric(subdata4$HRD)
subdata4$Cluster <- factor(subdata4$Cluster,levels = c("C1","C2","C3"))
subdata4$Cluster <- factor(subdata4$Cluster,levels = c("C1","C2"))
ggviolin(subdata4, x = "Cluster", y = "HRD",
         fill = "Cluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means(label = "p.signif")
ggsave(filename = "HRD.pdf",height = 5,width = 5)

####时序分析#####
## 加载需要用的包和检查数据
library(monocle)
library(DDRTree)
library(ggsci)
library(Seurat)
expr <- KIRC_mRNA_fpkm[,sub$Sample]


gene_annotation <- data.frame(geneID=rownames(expr),
                              gene_short_name=rownames(expr),
                              row.names = rownames(expr))

if ( ! identical(colnames(expr), sub$Sample)){
  data.frame(x=colnames(expr)[1:10],
             y=sub$Sample[1:10])
} else {
  cat("run successfully")
}
if (all(sub$Sample%in%colnames(expr))){
  expr <- expr[,match(sub$Sample,colnames(expr))]
} else{
  warning("unequal row number")
}

## 构建newCellDataSet
#expr <- 2^expr-1 
## 将log2 (FPKM+1) 转化为FPKM，这是为了使用作者推荐的方法DDTress
pd <- new("AnnotatedDataFrame", data=as.data.frame(sub))
fd <- new("AnnotatedDataFrame", data=as.data.frame(gene_annotation))
CDS <- newCellDataSet(as.matrix(expr), 
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily = tobit())

#CDS <- estimateSizeFactors(CDS)
#CDS <- estimateDispersions(CDS, cores=4, relative_expr = TRUE)
## 保留亚型间最显著变异的前500个基因
CDS <- detectGenes(cds = CDS, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(CDS),num_cells_expressed >= 0.1))

clustering_DEG_genes <- differentialGeneTest(CDS[expressed_genes,],
                                             fullModelFormulaStr = '~Cluster',
                                             cores=40)
ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:100]

CDS <- setOrderingFilter(CDS, ordering_genes = ordering_genes)
CDS <- reduceDimension(CDS, max_components = 4, method="DDRTress")
CDS <- orderCells(CDS)

plot_cell_trajectory(CDS, color_by = "State", cell_size = 2)+
  scale_colour_manual(values=pal_npg()(15))+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=11),
        legend.title = element_text(size=13),
        legend.text = element_text(size=13))


plot_cell_trajectory(CDS, color_by = "Cluster", cell_size = 2,
                     show_backbone = T,show_branch_points = F)+
  scale_colour_manual(values=pal_npg()(2))+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=11),
        legend.title = element_text(size=13),
        legend.text = element_text(size=13))

## C1的亚组
CDS1 <- CDS
CDS1$Cluster[CDS1$Cluster=='C1' & CDS1$State %in% 1]<- '1A'
CDS1$Cluster[CDS1$Cluster=='C1' & CDS1$State %in% c(2,3,4,5)] <- '1B'
#CDS1$Cluster[CDS1$Cluster=='C1' & CDS1$State %in% 3] <- '1C'
#CDS1$Cluster[CDS1$Cluster=='C1' & CDS2$State %in% c(0,2,3,4,6,8:15)] <- '1D'
CDS1$Cluster[CDS1$Cluster %in% 'C2'] <- 'Other'

pp <- plot_cell_trajectory(CDS1, color_by = "Cluster", cell_size = 2,
                           show_backbone = T,show_branch_points = F)+
  scale_colour_manual(values=c(pal_nejm()(4),'grey60'))+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=11),
        legend.title = element_text(size=13),
        legend.text = element_text(size=13))
pp

## C1亚组的生存曲线
tmp <- pp$data[,c(1,6)]
tmp <- merge(sub,tmp,by=1)
head(tmp)
library(survival)
library(survminer)
fit <- survfit(Surv(OS,EVENT) ~ Cluster.y,data = tmp)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black"))
ggsurvplot(fit,
           palette= c(pal_nejm()(4),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("1A","1B","Other"), 
           legend.title="C1",
           xlab="Time (years)",
           ylab='Overall survival')


## C2的亚组
## 这里为了方便解读，我仍然将C2分为两个亚组
CDS2 <- CDS
CDS2$Cluster[CDS2$Cluster=='C2' & CDS2$State %in% 3]<- '2A'
CDS2$Cluster[CDS2$Cluster=='C2' & CDS2$State %in% c(1,2,4,5)] <- '2B'
#CDS2$Cluster[CDS2$Cluster=='C2' & CDS2$State %in% 12] <- '2C'
#CDS2$Cluster[CDS2$Cluster=='C2' & CDS2$State %in% c(0,1:8,13,14,15)] <- '2D'
CDS2$Cluster[CDS2$Cluster %in% 'C1'] <- 'Other'

pp <- plot_cell_trajectory(CDS2, color_by = "Cluster", cell_size = 2,
                           show_backbone = T,show_branch_points = F)+
  scale_colour_manual(values=c(pal_nejm()(5),'grey60'))+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=11),
        legend.title = element_text(size=13),
        legend.text = element_text(size=13))
pp

## C2亚组的生存曲线
tmp2 <- pp$data[,c(1,6)]
tmp2 <- merge(sub,tmp2,by=1)

fit2 <- survfit(Surv(OS,EVENT) ~ Cluster.y,data = tmp2)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black"))
ggsurvplot(fit2,
           palette= c(pal_nejm()(6),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("2A","2B","Other"), 
           legend.title="C2",
           xlab="Time (years)",
           ylab='Overall survival')


plot_cell_trajectory(CDS,color_by = "Pseudotime")
ggsave("Pseudotime.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))


##使用seurat选择的高变基因
library(SeuratObject)
library(Seurat)
var.genes <- VariableFeatures(CDS)
mycds <- setOrderingFilter(CDS, var.genes)

disp_table <- dispersionTable(CDS)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
diff_test <- differentialGeneTest(mycds[disp.genes,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 1e-04))
p2 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=5,
                             show_rownames=T, return_heatmap=T)



p_TIME = plot_pseudotime_heatmap(CDS[ordering_genes[1:10],], num_clusters=3,
                                 show_rownames=T, return_heatmap=T)
ggsave("../时序分析结果/pseudotime_heatmap2.pdf", plot = p_TIME, width = 5, height = 10)

p_TIME = plot_pseudotime_heatmap(CDS[RNA_gene,], num_clusters=3,
                                 show_rownames=T, return_heatmap=T)
ggsave("../时序分析结果/pseudotime_heatmap3.pdf", plot = p_TIME, width = 5, height = 5)
####转录因子分析####
#Renal dysplasia  ARX,GATA3,GLI3,HOXD13,PAX2,SALL1,SIX1,SIX5,TBX1,TBX18,TP63,ZNF148,
#Renal hypoplasia GLI3,GTF2I,GTF2IRD1,HMGA2,HNF1B,MLXIPL,NFIA,PAX2,RREB1,SALL1,SALL4,SIX6,TBX1,TBX18,
#Renal hypoplasia/aplasia  ARID2,ARX,DEAF1,HNF1B,LHX1,MEOX1,PAX1,SALL4,SIX1,SIX5,SOX11,TP63,
#Renal agenesis FEZF1,HESX1,HNF1B,HNF4A,HOXD13,PAX2,SALL4,SIX1,SOX10,TFAP2A,TP63,ZIC3,
#Renal malrotation PAX2,SALL4,SIX1,
#Renal insufficiency FOS,GATA3,GCM2,GTF2I,GTF2IRD1,HNF1B,HOXA13,IKZF1,IRF5,KMT2A,LHX1,LMX1B,MAFB,MLXIPL,PAX2,PAX6,PPARG,SALL1,SIX1,SIX5,STAT4,TBX18,TCF4,WT1,ZNF423,ZNF592,
#Renal sarcoma TBX18,
#Crossed fused renal ectopia SALL4,WT1,
#Multiple renal cysts GTF2I,GTF2IRD1,RREB1,SALL1,TBX1,
#Renal cell carcinoma HNF1A,HNF1B,HNF4A,TFE3,TP53
#Renal duplication FOXC2,GTF2I,GTF2IRD1,
#Abnormality of the renal collecting system SIX1,
#Renal tubular dysfunction PDX1,PLAGL1,STAT3,ZFP57,
#Renal cyst GLI3,GLIS3,HNF1B,PRDM16,SKI,TFAP2A,ZNF148,
#Cystic renal dysplasia TBX18,
#Renal artery stenosis MLXIPL,STAT1,
#Papillary renal cell carcinoma FOXE1
#Unilateral renal agenesis GATA3,HNF1B,TBX1,
#Renal Fanconi syndrome HNF1B,HNF4A,
#Abnormal renal physiology HMGA2,
#Renal steatosis SIX1,
#Renal salt wasting NR0B1,


##EPAS1,ZEB2, 

library(RTN)
library(snow)
library(ComplexHeatmap)
library(mclust)
library(ClassDiscovery)
library(RColorBrewer)
library(gplots)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

pheno <- data.frame(CMOIC=sub$Cluster)
rownames(pheno) <- sub$Sample
tcgaBLCA <- KIRC_mRNA_fpkm[,rownames(pheno)]
tfs <- read.table("./easyinput_renal_regulon.txt",header = T)
tfs


# 取共有的基因名
regulatoryElements <- intersect(tfs$regulon, rownames(tcgaBLCA))

# 运行TNI构建程序
# we used the R package “RTN” to reconstruct transcriptional regulatory networks (regulons)
rtni_tcgaBLCA <- tni.constructor(expData = as.matrix(log2(tcgaBLCA + 1)), # 样图计算时候没有取对数, 
                                 regulatoryElements = regulatoryElements)

# 通过置换以及bootstrap计算reference regulatory network.
# mutual information analysis and Spearman rank-order correlation deduced the possible associations between a regulator and all potential target from the transcriptome expression profile, and permutation analysis was utilized to erase associations with an FDR > 0.00001. Bootstrapping strategy removed unstable associations through one thousand times of resampling with consensus bootstrap greater than 95%. 
# 这里量力而行设置多核，或者直接单核运算
options(cluster=snow::makeCluster(spec = 30, "SOCK")) # 打开4核并行计算（不确定是不是4核，不过我windows只用4，服务器我开12）
rtni_tcgaBLCA <- tni.permutation(rtni_tcgaBLCA, pValueCutoff = 1e-5, nPermutations = 1000)
rtni_tcgaBLCA <- tni.bootstrap(rtni_tcgaBLCA, nBootstraps = 1000)
stopCluster(getOption("cluster")) # 关闭并行计算

# 计算DPI-filtered regulatory network
# Data processing inequality filtering eliminated the weakest associations in triangles of two regulators and common targets
rtni_tcgaBLCA <- tni.dpi.filter(rtni_tcgaBLCA, eps = 0, sizeThreshold = TRUE, minRegulonSize = 15)

# 保存TNI对象以便后续分析
#save(rtni_tcgaBLCA, file="rtni_tcgaBLCA.RData")

# load("rtni_tcgaBLCA.RData")
# 计算每个样本的regulon活性
# Individual regulon activity was estimated by two-sided GSEA
rtnigsea_tcgaBLCA <- tni.gsea2(rtni_tcgaBLCA, regulatoryElements = regulatoryElements)
MIBC_regact <- tni.get(rtnigsea_tcgaBLCA, what = "regulonActivity")

# 保存活性对象
#save(MIBC_regact,file = "MIBC_regact.RData") 

# 加载活性对象
#(load("MIBC_regact.RData"))

# 设置颜色
clust.col <- c("#DD492E","#40548A","#32A087","#EC7D21")
blue <- "#5bc0eb"
gold <- "#ECE700"

plotdata <- standarize.fun(t(MIBC_regact$differential),halfwidth = 1.5) # 标准化regulon的活性
annCol.tcga <- pheno[order(pheno$CMOIC),,drop = F] # 构建样本注释信息，并对亚型进行排序
annColors.tcga <- list()
annColors.tcga[["CMOIC"]] <- c("C1" = clust.col[1],
                               "C2" = clust.col[2]
)
hcg <- hclust(distanceMatrix(as.matrix(MIBC_regact$differential[rownames(annCol.tcga),]), "euclidean"), "ward.D")
hm <- pheatmap(plotdata[hcg$order,rownames(annCol.tcga)],
               border_color = NA, # 热图单元格无边框
               color = colorpanel(64,low=blue,mid = "black",high=gold),
               cluster_rows = F, # 行不聚类
               cluster_cols = F, # 列聚类
               show_rownames = T, # 显示行名
               show_colnames = F, # 不显示列名
               gaps_col = c(231), # 亚型分割
               cellwidth = 0.8, # 固定单元格宽度
               cellheight = 10, # 固定单元格高度
               name = "MIBC Regulon", # 图例名字
               annotation_col = annCol.tcga[,"CMOIC",drop = F], # 样本注释
               annotation_colors = annColors.tcga["CMOIC"]) # 样本注释的对应颜色



pdf("regulon heatmap.pdf", width = 8,height = 6)
draw(hm) # 输出热图
invisible(dev.off())

