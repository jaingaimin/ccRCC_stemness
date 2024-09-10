




#hmdat <- as.data.frame(t(KIRC_mRNA_fpkm[rownames(type),sub$Sample]))

type <- read.csv("/home/data/vip39/data/immunegene5type.csv",header = T)
head(type)
table(type$Gene%in%rownames(KIRC_mRNA_fpkm))
table(rownames(sub)%in%rownames(hmdat))
type <- type[type$Gene%in%rownames(KIRC_mRNA_fpkm),]
type <- type[!duplicated(type$Gene),]
rownames(type) <- type$Gene
type <- type[-1]
head(type)
hmdat <- as.data.frame(t(KIRC_mRNA_fpkm[rownames(type),sub$Sample]))
head(hmdat)

comsam <- rownames(sub)
hmdat <- hmdat[comsam,]
risk <- sub[comsam,,drop = F]
dim(hmdat)
head(hmdat[,1:5])

# 拆分不同算法结果，获得类的名字
#immMethod <- sapply(strsplit(colnames(hmdat),"_",fixed = T),"[",2) #用easy_input.csv列名里的算法信息
immMethod <- type$Methods # 用easy_input_type.csv的算法那一列




library(ComplexHeatmap) 
# 最新版ComplexHeatmap好像有一些bug，使用里面的pheatmap函数会报错
# 因此我们从脚本直接加载pheatmap函数
source("~/huitu/FigureYa230immunelandscape/pheatmap_translate.R") # 位于当前文件夹，出自ComplexHeatmap_2.7.9.tar.gz

ht_opt$message = FALSE

# 创建注释
annCol <- data.frame(
  Subtype = risk$Cluster,
  row.names = rownames(risk),
  stringsAsFactors = F)
annRow <- data.frame(row.names = colnames(hmdat),
                     Methods = factor(immMethod,levels = unique(immMethod)),
                     stringsAsFactors = F)

annColors <- list(
  "Subtype" = c("C1"="#1F78B4FF", "C2"="#E31A1CFF"))

# 数据标准化
indata <- t(hmdat)
indata <- indata[,colSums(indata) > 0] # 确保没有富集全为0的细胞
plotdata <- standarize.fun(indata,halfwidth = 2)

# 样本按risk score排序
samorder <- c(rownames(sub)[sub$Cluster=="C1"],rownames(sub)[sub$Cluster=="C2"])

# 拆分各算法的结果
plotdata1 <- plotdata[rownames(annRow[which(annRow$Methods == "chemokine"),,drop = F]),]
plotdata2 <- plotdata[rownames(annRow[which(annRow$Methods == "chemokine_receptor"),,drop = F]),]
plotdata3 <- plotdata[rownames(annRow[which(annRow$Methods == "MHC"),,drop = F]),]
plotdata4 <- plotdata[rownames(annRow[which(annRow$Methods == "Immunoinhibitor"),,drop = F]),]
plotdata5 <- plotdata[rownames(annRow[which(annRow$Methods == "Immunostimulator"),,drop = F]),]


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
                name = "chemokine") # 为子热图的图例命名

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
                name = "chemokine_receptor")

hm3 <- pheatmap(mat = as.matrix(plotdata3[,samorder]),
                border_color = NA,
                color = blueyellow(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col =table(annCol$Subtype)[1],
                name = "MHC")

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
                name = "Immunoinhibitor")

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
                name = "Immunostimulator")


pdf("immune heatmap by ComplexHeatmap gene RNAmodify.pdf", width = 10,height = 30) # 保存前请注意RGUI里不能有任何显示的图像，否则不会pdf打不开
draw(hm1 %v% hm2 %v% hm3 %v% hm4 %v% hm5, # 垂直连接子热图
     heatmap_legend_side = "bottom", # 热图注释放底部
     annotation_legend_side = "bottom") # 顶部注释放底部
invisible(dev.off())





####绘制通路柱状图
(load("~/huitu/FigureYa146TMEbox/signature.RData"))
head(signature)

gsva_es_imm <- gsva(as.matrix(gsym.expr), signature)

head(gsva_es_imm)

# 分组
gsva_es_imm <- gsva_es_imm[,sub$Sample]
gsva_es_imm <- as.data.frame(t(gsva_es_imm))
gsva_es_imm_left <- data.frame(ID=sub$Sample,
                               TMEcluster=sub$Cluster)
pdata <- cbind(gsva_es_imm_left,gsva_es_imm)

head(pdata)
pdata[,3:ncol(pdata)] <- scale(pdata[,3:ncol(pdata)],scale = T,center = T)
#数据宽转长
pdata_melt <- melt(pdata,
                   id.vars = c("ID","TMEcluster"),
                   variable.name ="Signature",
                   value.name = "Signature_score")
head(pdata_melt)


# 使用ggplo2画图
c <- ggplot(pdata_melt,
            aes(x=Signature, y=Signature_score, 
                fill = TMEcluster, #按类填充颜色
                color = TMEcluster)) + #按类给边框着色
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", #outlier点用黑色
               outlier.size = 0.65) +
  #自定义配色
  scale_fill_manual(values= c("#E31A1C","#E7B800","#2E9FDF")) +
  #ggtitle("Gene signature score", "stratified by TME-cluster") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top")

# 标注p value
c + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))

# 标注*
# `****` = 1e-04, `***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1
p <- c + stat_compare_means(label = "p.signif")
p

ggsave(p, filename = "TME-relevant-sigature-boxplot.pdf", width = 13, height = 6)



###绘制典型免疫基因柱状图
c("CCL2","CD274","CD276","CTLA4","CXCR4","IL6","LAG3","PDCD1","TGFB1")
table(c("CCL2","CD274","CD276","CTLA4","CXCR4","IL6","LAG3","PDCD1","TGFB1")%in%rownames(KIRC_mRNA_fpkm))
gsva_es_imm <- KIRC_mRNA_fpkm[c("CCL2","CD274","CD276","CTLA4","CXCR4","IL6","LAG3","PDCD1","TGFB1"),sub$Sample]
gsva_es_imm <- as.data.frame(t(gsva_es_imm))
gsva_es_imm_left <- data.frame(ID=sub$Sample,
                               TMEcluster=sub$Cluster)
pdata <- cbind(gsva_es_imm_left,gsva_es_imm)

head(pdata)
pdata[,3:ncol(pdata)] <- scale(pdata[,3:ncol(pdata)],scale = T,center = T)
#数据宽转长
pdata_melt <- melt(pdata,
                   id.vars = c("ID","TMEcluster"),
                   variable.name ="Signature",
                   value.name = "Signature_score")
head(pdata_melt)


# 使用ggplo2画图
c <- ggplot(pdata_melt,
            aes(x=Signature, y=Signature_score, 
                fill = TMEcluster, #按类填充颜色
                color = TMEcluster)) + #按类给边框着色
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", #outlier点用黑色
               outlier.size = 0.65) +
  #自定义配色
  scale_fill_manual(values= c("#E31A1C","#E7B800","#2E9FDF")) +
  #ggtitle("Gene signature score", "stratified by TME-cluster") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top")

# 标注p value
c + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))

# 标注*
# `****` = 1e-04, `***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1
p <- c + stat_compare_means(label = "p.signif")
p

ggsave(p, filename = "TME-relevant-gene-boxplot.pdf", width = 13, height = 6)


###estimate


setwd("/home/data/vip39/TCGA/jiaowang_KIRC/met_model/met20211216/estimate")
dat=log2(KIRC_mRNA_fpkm[,sub$Sample]+1)
dat[1:4,1:4]
library(estimate)
estimate <- function(dat,pro){ 
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate) 
  filterCommonGenes(input.f=input.f, 
                    output.f=output.f ,
                    id="GeneSymbol") 
  estimateScore(input.ds = output.f,
                output.ds=output.ds, 
                platform="illumina")  
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro='Met'
scores=estimate(dat,pro)

scores <- as.data.frame(scores)
rownames(scores) <- str_replace_all(rownames(scores),"\\.","\\-")
table(gsva_es_imm_left$ID==rownames(scores))


pdata <- cbind(gsva_es_imm_left,scores)

head(pdata)

ggviolin(pdata, x = "TMEcluster", y = "ESTIMATEScore",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means(label = "p.signif")
ggsave(filename = "ESTIMATEScore.pdf",height = 5,width = 5)

ggviolin(pdata, x = "TMEcluster", y = "ImmuneScore",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means(label = "p.signif")
ggsave(filename = "ImmuneScore.pdf",height = 5,width = 5)

ggviolin(pdata, x = "TMEcluster", y = "StromalScore",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means( label = "p.signif")
ggsave(filename = "StromalScore.pdf",height = 5,width = 5)

# ggviolin(pdata, x = "TMEcluster", y = "ESTIMATEScore",
#          fill = "TMEcluster", palette = c("#48D1CC", "#E9967A"),alpha=.8,color = "white",add = "boxplot",
#          add.params = list(fill=c("#AFEEEE","#FFE4E1"))) + 
#   stat_compare_means()




pdata[,3:ncol(pdata)] <- scale(pdata[,3:ncol(pdata)],scale = T,center = T)
#数据宽转长
pdata_melt <- melt(pdata,
                   id.vars = c("ID","TMEcluster"),
                   variable.name ="Signature",
                   value.name = "Signature_score")
head(pdata_melt)



# 使用ggplo2画图
c <- ggplot(pdata_melt,
            aes(x=Signature, y=Signature_score, 
                fill = TMEcluster, #按类填充颜色
                color = TMEcluster)) + #按类给边框着色
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", #outlier点用黑色
               outlier.size = 0.65) +
  #自定义配色
  scale_fill_manual(values= c("#E31A1C","#E7B800","#2E9FDF")) +
  #ggtitle("Gene signature score", "stratified by TME-cluster") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top")

# 标注p value
c + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))

# 标注*
# `****` = 1e-04, `***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1
p <- c + stat_compare_means(label = "p.signif")
p

ggsave(p, filename = "TME-relevant-estimate-boxplot.pdf", width = 6, height = 6)



###tide####
#行是基因 列是样本 fpkm不需要log2处理
setwd("~/TCGA/jiaowang_KIRC/met_model/met20211216/TIDE/")
dat_TIDE <- KIRC_mRNA_fpkm[,sub$Sample]
TIDE <- dat_TIDE
# 这里为了得到比较好的结果，采用two direction median centered
TIDE <- sweep(TIDE,2, apply(TIDE,2,median,na.rm=T))
TIDE <- sweep(TIDE,1, apply(TIDE,1,median,na.rm=T))
write.table(TIDE,"TIDE_input.self_subtract",sep = "\t",row.names = T,col.names = NA,quote = F)


TIDE.res <- read.csv("~/TCGA/jiaowang_KIRC/met_model/met20211216/TIDE/met_TIDE_result.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)


###绘制不同分组的TIDE评分图
table(sub$Sample%in%rownames(TIDE.res))
TIDE.res <- TIDE.res[sub$Sample,]

gsva_es_imm_left <- data.frame(ID=sub$Sample,
                               TMEcluster=sub$Cluster)
pdata <- cbind(gsva_es_imm_left,TIDE.res)

head(pdata)
#save(pdata,file = "pdata_TIDE.rds")
# 默认方法为 method = "kruskal.test"
library(ggpubr)
# ggviolin(pdata, x = "TMEcluster", y = "TIDE",
#          fill = "TMEcluster", palette = c("#3CB371", "#48D1CC", "#E9967A"),alpha=.8,color = "white",add = "boxplot",
#          add.params = list(fill=c("#AFEEEE","#FFE4E1"))) + 
#   stat_compare_means()

ggviolin(pdata, x = "TMEcluster", y = "TIDE",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means(label = "p.signif")
ggsave(filename = "TIDE.pdf",height = 5,width = 5)
# ggviolin(pdata, x = "TMEcluster", y = "TIDE",
#          palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",
#          add = "boxplot",add.params = list(fill=c("#AFEEEE","#FFE4E1"))) + 
#   stat_compare_means()
ggviolin(pdata, x = "TMEcluster", y = "Dysfunction",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means( label = "p.signif")
ggsave(filename = "Dysfunction.pdf",height = 5,width = 5)

# ggviolin(pdata, x = "TMEcluster", y = "Dysfunction",
#          fill = "TMEcluster", palette = c("#48D1CC", "#E9967A"),alpha=.8,color = "white",add = "boxplot",
#          add.params = list(fill=c("#AFEEEE","#FFE4E1"))) + 
#   stat_compare_means()

ggviolin(pdata, x = "TMEcluster", y = "MSI Expr Sig",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means( label = "p.signif")
ggsave(filename = "MSI Expr Sig.pdf",height = 5,width = 5)



ggviolin(pdata, x = "TMEcluster", y = "Dysfunction",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means( label = "p.signif")
ggsave(filename = "Dysfunction.pdf",height = 5,width = 5)

ggviolin(pdata, x = "TMEcluster", y = "Exclusion",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means(label = "p.signif")
ggsave(filename = "Exclusion.pdf",height = 5,width = 5)

ggviolin(pdata, x = "TMEcluster", y = "CAF",
         fill = "TMEcluster", palette = c("#3CB371", "#48D1CC"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#90EE90","#AFEEEE"))) + 
  stat_compare_means( label = "p.signif")
ggsave(filename = "CAF.pdf",height = 5,width = 5)


# ggviolin(pdata, x = "TMEcluster", y = "MSI Expr Sig",
#          fill = "TMEcluster", palette = c("#48D1CC", "#E9967A"),alpha=.8,color = "white",add = "boxplot",
#          add.params = list(fill=c("#AFEEEE","#FFE4E1"))) + 
#   stat_compare_means()

# ggviolin(pdata, x = "TMEcluster", y = "Dysfunction",
#          color = "TMEcluster", palette = "neje",add = c("boxplot","jitter")) + 
#   stat_compare_means()
# ggviolin(pdata, x = "TMEcluster", y = "Exclusion",
#          color = "TMEcluster", palette = "neje",add = c("boxplot","jitter")) + 
#   stat_compare_means()
# ggviolin(pdata, x = "TMEcluster", y = "CAF",
#          color = "TMEcluster", palette = "neje",add = c("boxplot","jitter")) + 
#   stat_compare_means()
# ggviolin(pdata, x = "TMEcluster", y = "MSI Expr Sig",
#          color = "TMEcluster", palette = "neje",add = c("boxplot","jitter")) + 
#   stat_compare_means()
# meta = data.frame(Sample = paste0("sample",1:100),
#                   Response = sample(c("SD/PD","CR/PR"),100,replace = T),
#                   risk = sample(c("High","Low"),100,replace = T))
# head(meta)
# str(meta)
# library(dplyr)
# dat = count(meta,risk,Response)
# dat = dat %>% group_by(risk) %>%
#   summarise(Response = Response,n = n/sum(n))
# dat$Response = factor(dat$Response,levels = c("SD/PD","CR/PR"))
# dat

library(dplyr)
meta = data.frame(Sample=pdata$ID,
                  Response=pdata$Responder,
                  risk=pdata$TMEcluster)
meta2 = data.frame(Sample=pdata$ID,
                   Response=pdata$`No benefits`,
                   risk=pdata$TMEcluster)
save(meta,file = "metaTIDE.rds")
str(meta)
#meta$Response <- factor(meta$Response,levels = c("False","True"))
# head(meta)
# class(meta)
# meta <- as.matrix(meta)
# dat = count(meta,risk,Response)
# dat = dat %>% group_by(risk) %>% 
#   summarise(Response = Response,n = n/sum(n))
# dat$Response = factor(dat$Response,levels = c("False","True"))
# dat
# library(ggplot2)
# ggplot(data = dat)+
#   geom_bar(aes(x = risk,y = n,
#                fill = Response),
#            stat = "identity")+
#   scale_fill_manual(values =  c("#f87669","#2fa1dd"))+
#   geom_text(aes(x = risk,y = n,
#                 label = scales::percent(n)),
#             color = "white",
#             size = 6,
#             position = position_fill(vjust = 0.5))+
#   theme_minimal()+
#   theme(legend.position = "top")
#一句话搞定
library(ggstatsplot)
ggbarstats(data = meta,x=Response,y=risk)
ggsave(filename = "Immune_response.pdf",height = 6,width = 6)
ggbarstats(data = meta2,x=Response,y=risk)
ggsave(filename = "Immune_response2.pdf",height = 6,width = 6)
###治疗反应###


###TIP####
stepscore <- fread("~/database/TIP/KIRC/ssGSEA.normalized.score.txt",header = T,data.table = F)

head(stepscore[,1:4])
stepscore <- stepscore%>%
  column_to_rownames("Steps")
colnames(stepscore) <- str_sub(colnames(stepscore),1,15)
rownames(stepscore)[1:3] <- c("Step1.Release of cancer cell antigens",
                              "Step2.Cancer antigen presentation",
                              "Step3.Priming and activation")
rownames(stepscore)[c(21,22,23)] <- c("Step5.Infiltration of cancer cells into tumors",
                                      "Step6.Recognitioin of cancer cells by T cells",
                                      "Step7.Killing of cancer cells")
#stepscore <- as.data.frame(t(stepscore))


gsva_es_imm <- stepscore[,sub$Sample]
gsva_es_imm <- as.data.frame(t(gsva_es_imm))
gsva_es_imm_left <- data.frame(ID=sub$Sample,
                               TMEcluster=sub$Cluster)
pdata <- cbind(gsva_es_imm_left,gsva_es_imm)

head(pdata)
pdata[,3:ncol(pdata)] <- scale(pdata[,3:ncol(pdata)],scale = T,center = T)
#数据宽转长
pdata_melt <- melt(pdata,
                   id.vars = c("ID","TMEcluster"),
                   variable.name ="Signature",
                   value.name = "Signature_score")
head(pdata_melt)


# 使用ggplo2画图
c <- ggplot(pdata_melt,
            aes(x=Signature, y=Signature_score, 
                fill = TMEcluster, #按类填充颜色
                color = TMEcluster)) + #按类给边框着色
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", #outlier点用黑色
               outlier.size = 0.65) +
  #自定义配色
  scale_fill_manual(values= c("#E31A1C","#E7B800","#2E9FDF")) +
  #ggtitle("Gene signature score", "stratified by TME-cluster") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top")

# 标注p value
c + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))

# 标注*
# `****` = 1e-04, `***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1
p <- c + stat_compare_means(label = "p.signif")
p

ggsave(p, filename = "TME-relevant-TIP-boxplot.pdf", width = 20, height = 6)





####immunecell###
immunecell <- fread("~/database/TIP/KIRC/CIBER_KIRC_lm14_allsample_Result.txt",header = T,data.table = F)
head(immunecell)
immunecell$Mixture <- str_sub(immunecell$Mixture,1,15)
immunecell <- immunecell[!duplicated(immunecell$Mixture),]
head(immunecell)

immunecell <- immunecell[-c(16,17,18)]
rownames(immunecell) <- immunecell$Mixture
immunecell <- immunecell[-1]
head(immunecell)
immunecell <- as.data.frame(t(immunecell))



gsva_es_imm <- immunecell[,sub$Sample]
gsva_es_imm <- as.data.frame(t(gsva_es_imm))
gsva_es_imm_left <- data.frame(ID=sub$Sample,
                               TMEcluster=sub$Cluster)
pdata <- cbind(gsva_es_imm_left,gsva_es_imm)

head(pdata)
#pdata[,3:ncol(pdata)] <- scale(pdata[,3:ncol(pdata)],scale = T,center = T)
#数据宽转长
pdata_melt <- melt(pdata,
                   id.vars = c("ID","TMEcluster"),
                   variable.name ="Signature",
                   value.name = "Signature_score")
head(pdata_melt)


# 使用ggplo2画图
c <- ggplot(pdata_melt,
            aes(x=Signature, y=Signature_score, 
                fill = TMEcluster, #按类填充颜色
                color = TMEcluster)) + #按类给边框着色
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", #outlier点用黑色
               outlier.size = 0.65) +
  #自定义配色
  scale_fill_manual(values= c("#E31A1C","#E7B800","#2E9FDF")) +
  #ggtitle("Gene signature score", "stratified by TME-cluster") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top")

# 标注p value
c + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))

# 标注*
# `****` = 1e-04, `***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1
p <- c + stat_compare_means(label = "p.signif")
p

ggsave(p, filename = "TME-relevant-TIPimmunecell-boxplot.pdf", width = 13, height = 6)




####5个基因的相关性####
cordata <- as.data.frame(t(KIRC_mRNA_fpkm[unique(RNA_gene),sub$Sample]))
save(cordata,file = "./cor_RNAmodify_data.rds")

library(ggcorrplot)
library(corrplot)
#data("mtcars")
#计算pearson系数
dp <- cor(cordata,method = "pearson")
df=dp
#计算spearman系数
dsp <- cor(cordata,method = "spearman")

#清空df下三角的数据#
df[lower.tri(df,diag = F)] <- NA
#合并两种数据进入同一个矩阵
#把spearman下三角赋值到pearson矩阵的下三角
df[lower.tri(df)] <- dsp[lower.tri(dsp)]
#默认绘图
corrplot(df,addCoef.col = "black")

#计算p值
res1 <- cor.mtest(df,conf.level=.95)

#美化，就用这个版本的
corrplot.mixed(df,lower = "pie",upper = "square",
               addgrid.col = "black",
               p.mat=res1,insig="label_sig",tl.cex = 0.5,
               sig.level=c(.001,.01,.05),pch.cex=1.5,pch.col="black",
               tl.col="black",
               upper.col = colorRampPalette(c("#1E3163","#00C1D4","#FFED99","#FF7600"))(15),
               lower.col = colorRampPalette(c("#5D8233","#39A9CB","#FB9300","#E02401"))(15))
ggsave(filename = "agio_gene_cor.pdf", width = 8, height = 8)

###5个基因和ssgsea相关性####
### 用ssGSEA来量化浸润水平
### 1.加载marker
load(file = "/home/data/vip39/data/cellMarker_ssGSEA.Rdata")
### 2.加载表达量
#load(file = "tcga_panmRNA_expr.Rdata")
expr <- KIRC_mRNA_fpkm[,sub$Sample]
test <- expr[1:10,1:10]
expr <- as.matrix(expr)
library(GSVA)
### 挺耗时间的，调用了12个线程，17:12开始, 17:37结束
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")


tcga_gsva <- as.data.frame(t(gsva_data))

tcga_expr <- KIRC_mRNA_fpkm[,sub$Sample]


gene <- RNA_gene
immuscore <- function(gene){
  y <- as.numeric(tcga_expr[gene,])
  colnames <- colnames(tcga_gsva)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(tcga_gsva[,x]), y , method="spearman")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}

#以FOXP3为例，测试一下函数
immuscore("REN")
genelist <- RNA_gene
data <- do.call(rbind,lapply(genelist,immuscore))
head(data)
data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")
data$pstar[1:20]
ggplot(data, aes(immune_cells, gene)) + 
  geom_tile(aes(fill = cor), colour = "black",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=pstar),col ="black",size = 5)+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(angle = 45, hjust = 1),# 调整x轴文字
        axis.text.y = element_text(size = 8))+#调整y轴文字
  #调整legen
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))


###5个基因和M6A相关性####
RNA_modify <- c("TRMT61A","TRMT61B","TRMT10C","TRMT6","YTHDF3","YTHDC1","YTHDF1","YTHDF2","ALKBH1","ALKBH3",
                "NSUN7","NSUN3","TRDMT1","NSUN5","NOP2","DNMT1","NSUN6","NSUN4","DNMT3A","NSUN2","DNMT3B","TET2","ALYREF",
                "KIAA1429","ZC3H13","METTL14","CBLL1","RBM15","WTAP","METTL3","RBM15B","ALKBH5","FTO","YTHDF3","FMR1","YTHDC1","YTHDC2","HNRNPA2B1","ELAVL1","HNRNPC","YTHDF1","YTHDF2","IGF2BP1","LRPPRC")
table(RNA_modify%in%rownames(KIRC_mRNA_fpkm))
tcga_gsva <- as.data.frame(t(KIRC_mRNA_fpkm[RNA_modify,sub$Sample]))

tcga_expr <- KIRC_mRNA_fpkm[,sub$Sample]


gene <- RNA_gene
immuscore <- function(gene){
  y <- as.numeric(tcga_expr[gene,])
  colnames <- colnames(tcga_gsva)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(tcga_gsva[,x]), y , method="spearman")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}

#以FOXP3为例，测试一下函数
immuscore("REN")
genelist <- RNA_gene
data <- do.call(rbind,lapply(genelist,immuscore))
head(data)
data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")
data$pstar[1:20]
ggplot(data, aes(immune_cells, gene)) + 
  geom_tile(aes(fill = cor), colour = "black",size=1)+
  scale_fill_gradient2(low = "#1F78B4FF",mid = "white",high = "#6A3D9AFF")+
  geom_text(aes(label=pstar),col ="black",size = 5)+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(angle = 45, hjust = 1),# 调整x轴文字
        axis.text.y = element_text(size = 8))+#调整y轴文字
  #调整legen
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))

####不同分期预后分析####
T1_sub <- sub%>%
  filter(sub$pstage=="T1")
meta$pstage <- sub$pstage
T1_sub <- meta%>%
  filter(meta$pstage=="T1")

sfit3 <- survfit(Surv(OS,EVENT)~Cluster,data = T1_sub)
sfit3 <- survfit(Surv(event,event)~cluster,data = T1_sub)
ggsurvplot(sfit3, conf.int=F, pval=TRUE)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black"))
ggsurvplot(sfit3,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("C1","C2"), 
           legend.title="cluster",
           xlab="Time (days)",
           ylab='Overall survival',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)


T2_sub <- sub%>%
  filter(sub$pstage=="T2")

sfit4 <- survfit(Surv(OS,EVENT)~Cluster,data = T2_sub)
ggsurvplot(sfit4, conf.int=F, pval=TRUE)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black"))
ggsurvplot(sfit4,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("C1","C2"), 
           legend.title="cluster",
           xlab="Time (days)",
           ylab='Overall survival',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)


T3_sub <- sub%>%
  filter(sub$pstage=="T3")

sfit5 <- survfit(Surv(OS,EVENT)~Cluster,data = T3_sub)
ggsurvplot(sfit5, conf.int=F, pval=TRUE)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black"))
ggsurvplot(sfit5,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("C1","C2"), 
           legend.title="cluster",
           xlab="Time (days)",
           ylab='Overall survival',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)

T4_sub <- sub%>%
  filter(sub$pstage=="T4")

sfit6 <- survfit(Surv(OS,EVENT)~Cluster,data = T4_sub)
ggsurvplot(sfit6, conf.int=F, pval=TRUE)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black"))
ggsurvplot(sfit6,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("C1","C2"), 
           legend.title="cluster",
           xlab="Time (days)",
           ylab='Overall survival',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)

##
class(allclin)
allclin <- as.data.frame(allclin)
rownames(allclin) <- allclin$sample

sub$PFI <- allclin[sub$Sample,"PFI"]
sub$PFItime <- allclin[sub$Sample,"PFI.time"]




###临床信息表格####
#196模仿数据绘图
pd <- read_tsv("~/database/KIRC_DATA/TCGA-KIRC.GDC_phenotype.tsv.gz")
pd <- as.data.frame(pd)
pd$idnew <- str_sub(pd$submitter_id.samples,1,15)
pd <- pd[!duplicated(pd$idnew),]
rownames(pd) <- pd$idnew
head(pd[,1:4])
load("~/data/TCGA-KIRC_clin_info.Rdata")
dat <- data.frame(Risk=ifelse(sub$Cluster=="C2","High","Low"),
                  Status=ifelse(sub$EVENT==0,"Live","Dead"),
                  M=pd[sub$Sample,"pathologic_M"],
                  N=pd[sub$Sample,"pathologic_N"],
                  T=str_sub(pd[sub$Sample,"pathologic_T"],1,2),
                  Stage=sub$pstage)
dat <- na.omit(dat)

library(dplyr)
library(ggplot2)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


# 按Risk分成High和Low，计算各列数值。
gname <- "Risk"
vname <- setdiff(colnames(dat), gname)
pie.high <- pie.low <- list()
fisher.p <- c()
for (i in vname) {
  tmp <- table(dat[,gname], dat[,i])
  p <- format(fisher.test(tmp)$p.value,digits = 2)
  names(p) <- i
  fisher.p <- c(fisher.p, p)
  
  pie.dat <- 
    tmp %>% as.data.frame() %>% group_by(Var1) %>% mutate(Pct = Freq/sum(Freq)) %>% as.data.frame()
  
  # 表格内的两行对应Risk的两类：Risk high和Risk low
  pie.high[[i]] <- pie.dat[which(pie.dat$Var1 == "High"),]
  pie.low[[i]] <- pie.dat[which(pie.dat$Var1 == "Low"),]
}


# 设置颜色
black  <- "#1E1E1B"
blue   <- "#3C4E98"
yellow <- "#E4DB36"
orange <- "#E19143"
green  <- "#57A12B"
cherry <- "#8D3A86"

# 创建颜色
status.col <- c("grey80",black)
stage.col <- alpha(blue, c(0.4, 0.6, 0.8, 1))
M.col <- c(yellow, orange)
N.col <- alpha(green, c(0.5, 0.7, 1))
T.col <- alpha(cherry, c(0.4, 0.6, 0.8, 1))

# 硬核base plot一块一块画，当然也可以把其中的pie chart提取出来后期AI或者PPT拼接也是比较方便的
pdf("pieTable.pdf",width = 7, height = 5)
showLayout <- F # 默认不在最终pdf的首页显示layout结构，不过建议初次绘制的时候改为TRUE看一下，方便理解

# 设置画面布局，相同数字代表同一区块，数字越多代表该区块所占面积越大（一共25个区域）
layout(matrix(c( 1, 1, 1,  2, 2, 2,  3, 3, 3,  4, 4, 4,  5, 5, 5,  6, 6, 6,
                 7, 7, 7,  8, 8, 8,  9, 9, 9, 10,10,10, 11,11,11, 12,12,12,
                 7, 7, 7,  8, 8, 8,  9, 9, 9, 10,10,10, 11,11,11, 12,12,12,
                 13,13,13, 14,14,14, 15,15,15, 16,16,16, 17,17,17, 18,18,18,
                 13,13,13, 14,14,14, 15,15,15, 16,16,16, 17,17,17, 18,18,18,
                 19,19,19, 20,20,20, 21,21,21, 22,22,22, 23,23,23, 24,24,24,
                 25,25,25, 25,25,25, 25,25,25, 25,25,25, 25,25,25, 25,25,25),
              byrow = T,nrow = 7))

if(showLayout) {
  layout.show(n = 25) # 直观展示画布分布
}

#-------------------------#
# 画布区域1-6：绘制图抬头 #
#-------------------------#

par(bty="n", mgp = c(0,0,0), mar = c(0,0,0,0), lwd = 2) # 基础参数，各边界距离为0
plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "ccRCC",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Status",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Stage",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "M",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "N",cex = 2, col = "white") # 显示图标题

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "T",cex = 2, col = "white") # 显示图标题

#--------------------------------------#
# 画布区域7-12：绘制High组抬头和扇形图 #
#--------------------------------------#

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "High\n(n = 100)",cex = 2, col = "white") # 显示图标题

# High group
pie(pie.high$Status$Pct, 
    col = status.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$Stage$Pct, 
    col = stage.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$M$Pct, 
    col = M.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$N$Pct, 
    col = N.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.high$`T`$Pct, 
    col = T.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
abline(v = par("usr")[2], col = "black") # 右侧封上黑线

#--------------------------------------#
# 画布区域13-18：绘制Low组抬头和扇形图 #
#--------------------------------------#

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     "Low\n(n = 200)",cex = 2, col = "white") # 显示图标题

# Low group
pie(pie.low$Status$Pct, 
    col = status.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$Stage$Pct, 
    col = stage.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$M$Pct, 
    col = M.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$N$Pct, 
    col = N.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)

pie(pie.low$`T`$Pct, 
    col = T.col, 
    border = "white",  
    radius = 1, 
    labels = NA,
    init.angle = 90)
symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
abline(v = par("usr")[2], col = "black") # 右侧封上黑线

#--------------------------------#
# 画布区域19-24：绘制空抬头和p值 #
#--------------------------------#

plot(1,1,
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black") # 背景涂黑

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["Status"]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black") # 底部封上黑线

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["Stage"]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black")

plot(1,1,col = "white",
     xlab = "",xaxt = "n",# 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["M"]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black")

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["N"]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black")

plot(1,1,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
     (par("usr")[3]+par("usr")[4])/2,
     paste0("p = ",fisher.p["T"]),cex = 1.5, col = "black") # 显示图标题
abline(h = par("usr")[3], col = "black") # 底部封上黑线
abline(v = par("usr")[2], col = "black") # 右侧封上黑线

#----------------------#
# 画布区域25：绘制图例 #
#----------------------#

plot(0,0,col = "white",
     xlab = "",xaxt = "n", # 不显示x坐标轴
     ylab = "",yaxt = "n") # 不显示y坐标轴
legend("topleft",
       legend = c("Alive","Dead",
                  "I","II","III","IV",
                  "M0","M1",
                  "N0","N1","N2",
                  "T1","T2","T3","T4"),
       fill = c(status.col,
                stage.col,
                M.col,
                N.col,
                T.col),
       border = NA, # 图例颜色没有边框
       bty = "n", # 图例没有边框
       cex = 1.2,
       #box.lwd = 3,
       x.intersp = 0.05,
       y.intersp = 1,
       text.width = 0.075, # 图例的间隔
       horiz = T) # 图例水平放置

# 关闭图像句柄
invisible(dev.off())



####每个基因对于预后的评估####
##OS###


##PFI###



###准备scissor输入数据####
bulk_dataset <- KIRC_mRNA_fpkm[,sub$Sample]
bulk_survival <- data.frame(TCGA_patient_barcode=sub$Sample,
                            OS_time=sub$OS,
                            Status=sub$EVENT,
                            Cluster=sub$Cluster)
save(bulk_dataset,bulk_survival,file = "met_for_scissor.rds")
