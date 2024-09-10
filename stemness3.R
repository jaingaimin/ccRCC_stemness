
#save.image(file = "~/TCGA/RNAmodify/RNAmodifyalldata.rds")
#load(file = "~/TCGA/RNAmodify/RNAmodifyalldata.rds")
####japan队列复现MTCS2和MTCS1####


##预后模型 基于单因素基因####
#使用
dir.create("prognosis_model")
setwd("./prognosis_model/")
library(survival)
library(randomForestSRC)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN) == 0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}


# 去除无表达的基因，log变换，z-score
head(expr)
dim(expr)
#expr <- KIRC_mRNA_fpkm[Coxoutput$gene,]
##
#setwd("../validation_marker/")
expr <- KIRC_mRNA_fpkm[unique(c(marker.up$templates$probe,marker.up$templates$probe)),]

expr <- as.data.frame(round(t(scale(t(log2(expr + 1)))),3))
expr <- expr[,rownames(sub)]
colnames(sub)
str(sub)
head(sub)
surv <- sub[,c("OS","EVENT")]
colnames(surv) <- c("OS.time","OS")
head(surv)
head(sub)
cox.pcutoff <- 0.05 # cox的p阈值
Coxoutput.OS <- NULL
for (i in 1:nrow(expr)) {
  display.progress(index = i,totalN = nrow(expr)) # 显示进度
  
  # 产生临时变量存储生存以及变量表达值
  tmp <- data.frame(gene = as.numeric(expr[i,]),
                    OS.time = surv[,"OS.time"],
                    OS = surv[,"OS"],
                    stringsAsFactors = F)
  
  # 单变量cox比例风险模型
  cox <- coxph(Surv(OS.time, OS) ~ gene, data = tmp)
  coxSummary = summary(cox)
  
  # 生成cox结果数据框，包括基因名，风险比，z值，waldtest p值，以及HR置信区间
  Coxoutput.OS=rbind.data.frame(Coxoutput.OS,data.frame(gene=rownames(expr)[i],
                                                        HR=as.numeric(coxSummary$coefficients[,"exp(coef)"]),
                                                        z=as.numeric(coxSummary$coefficients[,"z"]),
                                                        pvalue=as.numeric(coxSummary$coefficients[,"Pr(>|z|)"]),
                                                        lower=as.numeric(coxSummary$conf.int[,3]),
                                                        upper=as.numeric(coxSummary$conf.int[,4]),
                                                        stringsAsFactors = F),
                                stringsAsFactors = F)
}
head(Coxoutput.OS)
write.csv(Coxoutput.OS,"univariate cox regression for gene filtering.csv",row.names = F,quote = F)

##随机森林进一步降维
gene.sel <- Coxoutput.OS[which(Coxoutput.OS$pvalue < cox.pcutoff),"gene"]
tmp <- expr[gene.sel,]

rownames(tmp) <- gsub("-","_",rownames(tmp)) # 防止出现“-”导致程序报错
dt.rf <- cbind.data.frame(surv[,c("OS","OS.time")],t(tmp))

ntree <- 1000
surv.rf <- rfsrc(Surv(OS.time, OS) ~ ., 
                 data = dt.rf, 
                 ntree = ntree,
                 importance = TRUE,
                 seed = 12345678)

#排列组化确定最优签名
num.imp <- 10
rel.imp <- sort(surv.rf$importance, decreasing = T)
rel.imp.sel <- rel.imp[1:num.imp] # 取出一定数量的基因
names(rel.imp.sel) <- gsub("_","-",names(rel.imp.sel)) # 还原基因名

outTab <- NULL
n.sum <- 0
for (i in 1:num.imp) {
  cat(paste0("combination using ",i," genes...\n"))
  tmp <- utils::combn(names(rel.imp.sel), m=i) # 获取当前基因个数下的排列组合
  n <- ncol(tmp)
  for (j in 1:n) {
    combgene <- tmp[,j] # 取出每一次组合的基因名
    combexpr <- cbind.data.frame(t(expr[combgene,]), # 构建数据库做多变量cox
                                 OS.time = surv[,"OS.time"],
                                 OS = surv[,"OS"],
                                 stringsAsFactors = F)
    cox <- coxph(Surv(OS.time, OS) ~ ., data = combexpr)
    coxSummary <- summary(cox)
    coeff <- coxSummary$coefficients[,1] # 取出系数
    riskscore <- as.matrix(combexpr[,combgene]) %*% coeff # 计算riskscore
    riskscore <- data.frame(riskscore = as.numeric(riskscore[,1]),
                            group = ifelse(riskscore[,1] > median(riskscore[,1]),"HRisk","LRisk"), # 根据中位数分组
                            row.names = rownames(riskscore),
                            OS.time = combexpr$OS.time,
                            OS = combexpr$OS,
                            stringsAsFactors = F)
    fitd <- survdiff(Surv(OS.time, OS) ~ group,
                     data = riskscore,
                     na.action = na.exclude)
    p.val <- 1-pchisq(fitd$chisq, length(fitd$n) - 1) # log-rank检验
    
    outTab <- rbind.data.frame(outTab,
                               data.frame(num.gene = ifelse(i == 1, paste0(i," gene"), paste0(i," genes")), # 当前基因数目
                                          km.pvalue = p.val, # KM曲线p值
                                          core.gene = paste(combgene,collapse = " | "), # 该组合下的基因
                                          stringsAsFactors = F),
                               stringsAsFactors = F)
  }
  n.sum <- n + n.sum # 校验排列组合的总数
}
if(n.sum == 2^num.imp-1) { # 如果总和不等则报错
  write.csv(outTab,"combination of important genes with KM pvalues.csv",row.names = F,quote = F)
} else (message("Wrong combination!!!"))


sigpoints <- Coxoutput.OS[which(Coxoutput.OS$pvalue < cox.pcutoff),]
unsigpoints <- Coxoutput.OS[which(Coxoutput.OS$pvalue >= cox.pcutoff),]

pdf("volcano.pdf",width = 5,height = 5)
par(bty = "o", mgp = c(2,.6,0), mar = c(3,3,1,1), las = 1, font.axis = 1) # 基础参数
plot(log(Coxoutput.OS$HR),
     -log10(Coxoutput.OS$pvalue),
     xlab = "Univariate Cox coefficient",
     ylab = bquote("-log"[10]~"(P value)"),
     xlim = c(-2,2))
points(log(sigpoints$HR),
       -log10(sigpoints$pvalue),
       col = ggplot2::alpha("#E53435",0.8),
       pch = 19)
points(log(unsigpoints$HR),
       -log10(unsigpoints$pvalue),
       col = ggplot2::alpha("#21498D",0.8),
       pch = 19)
abline(h = -log10(cox.pcutoff), lty = 2, col = "grey60")
invisible(dev.off())


xrange <- range(pretty(range(rel.imp.sel))) # 根据重要性区间确定x轴范围
yrange <- c(1,length(rel.imp.sel))  # 根据重要变量个数确定y轴范围

pdf("variable importance.pdf",width = 5,height = 5)
par(bty = "o", mgp = c(1.5,.33,0), mar = c(3,7,1,2),las = 1, tcl = -.25)
plot(NULL,NULL,
     xlim = xrange,
     ylim = yrange,
     xlab = "Variable Importance",
     ylab = "",
     yaxt = "n",
     las = 1)
axis(side = 2,at = 1:length(rel.imp.sel),rev(names(rel.imp.sel))) # 补齐y轴
for (i in 1:length(rel.imp.sel)) { # 循环添加线
  lines(c(xrange[1],rev(rel.imp.sel)[i]),
        c(i,i),
        lwd = 2.5,
        col = "#33A02CFF")
}
invisible(dev.off())


## 绘制排列组合p值


num.comb <- 20
outTab2 <- outTab[order(outTab$km.pvalue),][1:num.comb,]
head(outTab2)
#outTab2$km.pvalue[1] <- 1.100223e-18
pdf("combination barplot.pdf",width = 5,height = 5)
par(bty = "o", mgp = c(1.5,.33,0), mar = c(1,4,3,1),las = 1, tcl = -.25)
barplot(rev(-log10(outTab2$km.pvalue)),
        horiz = T, # 柱状图横向
        names.arg = rev(outTab2$num.gene), # 添加y轴名称
        xaxt = "n", # 取消下方x轴
        col = "#FF7F00FF")
axis(side = 3) # 在上方添加x轴
mtext(side = 3, bquote("-log"[10]~"(P value)"), line = 1) # 添加x轴名称
invisible(dev.off())
outTab2$core.gene[1]

#"MECP2 | AR | IRF3 | CBX8 | KL"
"PTPRB","RTL1","PAGE2B", "MGAM"
#"PTPRB","RTL1","PAGE2B", "MGAM"
#PDIA2 | OR4C6 | SFRP5 | BARX1 | GJB6
###风险三联图
#riskscore计算
colnames(Coxoutput.OS)
rownames(Coxoutput.OS) <- Coxoutput.OS$gene
index.min <- as.numeric(Coxoutput.OS[c("PTPRB","RTL1","PAGE2B", "MGAM"),"z"])
signature <- as.matrix(t(KIRC_mRNA_fpkm[c("PTPRB","RTL1","PAGE2B", "MGAM"),sub$Sample])) %*% as.matrix(index.min) 
summary(signature)
class(signature)
colnames(signature)[1] <- "lasso"

table(rownames(signature)==rownames(sub))
sub$lasso <- signature[,1]

library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

bestvars <- c("PTPRB","RTL1","PAGE2B", "MGAM")
# risk score，用于画顶部散点图
rs <- sub$lasso
names(rs) <- rownames(sub)
rs_data <- data.frame(x=1:length(rs),rs=as.numeric(sort(rs)))
# 用中值分组
rs_data$Risk <- ifelse(rs_data$rs>=median(rs_data$rs), "High-risk", "Low-risk")
head(rs_data)

sub$OS.time=sub$OS
sub$OS=sub$EVENT
# follow-up，用于画中间B图
surv_data <- data.frame(x=1:length(rs),
                        t=sub[names(sort(rs)),'OS.time']/365*12,
                        s=sub[names(sort(rs)),'OS']) 
surv_data$Status <- as.factor(ifelse(surv_data$s==0,'Alive','Death'))
head(surv_data)

# 提取signature对应的data，并按risk score排序，用于画底部热图
exp_data <- as.data.frame(t(KIRC_mRNA_fpkm[bestvars,names(sort(rs))]))
exp_data[1:2,1:4]


plot.A <- ggplot(rs_data, aes(x=x,y=rs))+
  geom_point(aes(col=Risk),size=0.5)+
  scale_color_manual(labels=c("High-risk","Low-risk"), 
                     #guide_legend(guide = NULL), #如果不想画图例就删掉#
                     name="Risk score", values =c("#DC0000FF", "#00A087FF")) + 
  
  # 画竖向虚线
  geom_segment(aes(x = sum(rs_data$Risk=="Low-risk"),
                   y = 0, 
                   xend = sum(rs_data$Risk=="Low-risk"), 
                   yend = max(rs_data$rs)), linetype="dashed", size = 0.6)+
  # 画横线
  #geom_segment(aes(x=0,y=median(rs_data$rs),
  #                 xend=nrow(rs_data),
  #                 yend=median(rs_data$rs)),linetype="dashed", size = 0.3)+
  
  # 写文字Cutoff:
  #geom_text(aes(x=sum(rs_data$Risk=="Low-risk")/2,
  #              y=median(rs_data$rs)+8,
  #              label=paste0("Cutoff: ",round(median(rs_data$rs),3))),
  #          col ="black",size = 4,alpha=0.8)+
  
theme(axis.title.x=element_blank()) +
  scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
  #scale_y_continuous(limits = c(0,NA),expand = c(-10,0))+
  labs(y="Risk score",x="",fill="Risk") +
  #scale_colour_discrete(name="Risk scores") +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样画坐标轴，就删掉这行
        axis.text.x=element_blank())

plot.A

plot.B <- ggplot(surv_data,aes(x=x,y=t))+
  geom_point(aes(col=Status),size=0.5)+
  geom_vline(aes(xintercept=sum(rs_data$Risk=="Low-risk")),size=0.6,linetype="dashed")+
  scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
  scale_color_manual(labels=c("Alive","Dead"),
                     values =c("#00A087FF","#DC0000FF"))+
  labs(y="RFS(months)",x="")+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样不画坐标轴，就删掉前面的#
        axis.text.x=element_blank())

plot.B


tmp <- t(scale(exp_data))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
reorder_cormat <- function(cormat){
  dd <- dist(cormat)
  hc <- hclust(dd,method = "average")
  cormat <-cormat[hc$order,]
}
tmp1 <- reorder_cormat(tmp)
tmp1 <- rbind(tmp1,ifelse(rs_data$Risk=="Low-risk",-1.5,1.5))
tmp.m <- melt(tmp1)

p2 <-ggplot(tmp.m, aes(Var2, Var1),size=0.5) + 
  geom_tile(aes(fill = value)) 

plot.C <- p2 + scale_fill_gradient2(name="Genes\nexpression", low="#00A087FF", high="#DC0000FF", mid="white") +
  labs(x = "", y = "")+
  theme_classic()+
  theme(legend.title = element_text(size = 12), legend.position = "right",
        axis.line = element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank())

plot.C

plot_grid(plot.A, plot.B, plot.C,
          labels = c("B", "",""), # 或者按顺序标注ABC
          rel_heights = c(1,1,1), # 3个图的比例
          #label_x=0,
          #label_y=1,
          align = 'v',ncol = 1, axis="lr", scale = c(1,1,1), greedy = F)

# 保存到文件
ggsave("coxTCGA.pdf", width = 7, height = 9)

##KM曲线
sub3 <- sub
sub3$risk_group <- ifelse(sub3$lasso>=median(sub3$lasso), "High-risk", "Low-risk")
colnames(sub3)
sub3$OS.time <- sub3$OS.time/365

sfit <- survfit(Surv(OS.time,OS)~risk_group,data = sub3)
sub3$PFI.time <- allclin[rownames(sub3),"PFI.time"]
sub3$PFI <- allclin[rownames(sub3),"PFI"]
sub3$PFI.time <- sub3$PFI.time/365
sfit <- survfit(Surv(PFI.time,PFI)~risk_group,data = sub3)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black"))
ggsurvplot(sfit,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("High-risk", "Low-risk"), 
           legend.title="type",
           xlab="Time (years)",
           ylab='Overall survival',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)

ggsurvplot(sfit,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("High-risk", "Low-risk"), 
           legend.title="type",
           xlab="Time (days)",
           ylab='OS',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)

ggsurvplot(sfit,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("High-risk", "Low-risk"), 
           legend.title="type",
           xlab="Time (years)",
           ylab='PFI',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)

Sys.setenv(LANGUAGE = "en")
library(tidyr)
library(dplyr)
library(tibble)
library(timeROC)
library(survival)
library(pROC)
library(ggplot2)
library(ggtext)

riskRoc <- timeROC(T = sub3$OS.time,delta = sub3$OS,
                   marker = sub3$lasso,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果

oldriskDat <- sub
oldriskDat$pstage
oldriskDat$grade <- pd[rownames(oldriskDat),"neoplasm_histologic_grade"]
head(oldriskDat)
oldriskDat$futime <- oldriskDat$OS.time/365
head(riskDat)
riskDat <- sub
riskDat$futime <- riskDat$OS.time/365

newcindexRoc <- timeROC(T = riskDat$futime,delta = riskDat$EVENT,
                        marker = riskDat$lasso,cause = 1,
                        weighting="marginal",
                        times = seq(0,5,5/100))
oldcindexRoc <- timeROC(T = oldriskDat$futime,delta = oldriskDat$EVENT,
                        marker =as.numeric(as.factor(oldriskDat$pstage)),cause = 1,
                        #marker = oldriskDat$pstage,cause = 1,
                        weighting="marginal",
                        times = seq(0,5,5/100))
CindexData <- data.frame(newCindex=newcindexRoc$AUC,
                         oldCindex=oldcindexRoc$AUC,
                         time=newcindexRoc$times)

CindexData <- gather(CindexData,"model",'Cindex',-time)
ggplot(data = CindexData,mapping = aes(x=time,y=Cindex,col=model))+
  geom_line(na.rm = T,lty="solid",size=1)+
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + #设置坐标起点
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+    #设置坐标起点
  theme_classic()+
  scale_color_manual(values=c('#fb6974','#ae86ff'),        #线条颜色
                     name = "model",                       #图例标题
                     #breaks = c("newCindex", "oldCindex"),#分组，可有可无
                     labels = paste("<span style='color:",
                                    c('#fb6974','#ae86ff'),#颜色
                                    "'>",
                                    c("new Model", "old Model"),#图例标签名
                                    "</span>",
                                    sep = "")) +
  theme(legend.text = element_markdown())+                      #服从markdown格式
  xlab("time(years)")+ylab("C-index")

###JAPAN验证####
load("/home/data/vip39/database/KIRC_validation_data/japan_input.rds")
japan_clin <- fread("/home/data/vip39/database/KIRC_validation_data/japan_clin_clean.csv",header=T,data.table=F)
rownames(japan_clin) <- japan_clin$sample_ID
colnames(japan_clin)
sub5 <- japan_clin[,c("sample_ID","outcome","month")]
colnames(sub5) <- c("ID","OS","OS.time")
str(sub5)
sub5$OS <- ifelse(sub5$OS=="alive",0,1)
rownames(sub5) <- sub5$ID
sub5 <- sub5[colnames(japan_mRNA_fpkm),]
japan_mRNA_fpkm <- japan_mRNA_fpkm%>%
  column_to_rownames("gene_name")
#riskscore计算
colnames(Coxoutput.OS)
rownames(Coxoutput.OS) <- Coxoutput.OS$gene
index.min <- as.numeric(Coxoutput.OS[c("PTPRB","RTL1","PAGE2B", "MGAM"),"z"])
#ICGC_mRNA_rpkm <- as.data.frame(round(t(scale(t(log2(ICGC_mRNA_rpkm + 1)))),3))
signature3 <- as.matrix(t(japan_mRNA_fpkm[c("PTPRB","RTL1","PAGE2B", "MGAM"),sub5$ID])) %*% as.matrix(index.min) 
summary(signature3)
class(signature3)
colnames(signature3)[1] <- "lasso"

table(rownames(signature3)==rownames(sub5))
sub5$lasso <- signature3[,1]

library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

bestvars <- c("PTPRB","RTL1","PAGE2B", "MGAM")
# risk score，用于画顶部散点图
rs <- sub5$lasso
names(rs) <- rownames(sub5)
rs_data <- data.frame(x=1:length(rs),rs=as.numeric(sort(rs)))
# 用中值分组
rs_data$Risk <- ifelse(rs_data$rs>=median(rs_data$rs), "High-risk", "Low-risk")
head(rs_data)

# follow-up，用于画中间B图
surv_data <- data.frame(x=1:length(rs),
                        t=sub5[names(sort(rs)),'OS.time'],
                        s=sub5[names(sort(rs)),'OS']) 
surv_data$Status <- as.factor(ifelse(surv_data$s==0,'Alive','Death'))
head(surv_data)

# 提取signature对应的data，并按risk score排序，用于画底部热图
exp_data <- as.data.frame(t(japan_mRNA_fpkm[bestvars,names(sort(rs))]))
exp_data[1:2,1:4]


plot.A <- ggplot(rs_data, aes(x=x,y=rs))+
  geom_point(aes(col=Risk),size=0.5)+
  scale_color_manual(labels=c("High-risk","Low-risk"), 
                     #guide_legend(guide = NULL), #如果不想画图例就删掉#
                     name="Risk score", values =c("#DC0000FF", "#00A087FF")) + 
  
  # 画竖向虚线
  geom_segment(aes(x = sum(rs_data$Risk=="Low-risk"),
                   y = 0, 
                   xend = sum(rs_data$Risk=="Low-risk"), 
                   yend = max(rs_data$rs)), linetype="dashed", size = 0.6)+
  # 画横线
  #geom_segment(aes(x=0,y=median(rs_data$rs),
  #                 xend=nrow(rs_data),
  #                 yend=median(rs_data$rs)),linetype="dashed", size = 0.3)+
  
  # 写文字Cutoff:
  #geom_text(aes(x=sum(rs_data$Risk=="Low-risk")/2,
  #              y=median(rs_data$rs)+8,
  #              label=paste0("Cutoff: ",round(median(rs_data$rs),3))),
  #          col ="black",size = 4,alpha=0.8)+
  
theme(axis.title.x=element_blank()) +
  scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
  labs(y="Risk score",x="",fill="Risk") +
  #scale_colour_discrete(name="Risk scores") +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样画坐标轴，就删掉这行
        axis.text.x=element_blank())

plot.A

plot.B <- ggplot(surv_data,aes(x=x,y=t))+
  geom_point(aes(col=Status),size=0.5)+
  geom_vline(aes(xintercept=sum(rs_data$Risk=="Low-risk")),size=0.6,linetype="dashed")+
  scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
  scale_color_manual(labels=c("Alive","Dead"),
                     values =c("#00A087FF","#DC0000FF"))+
  labs(y="RFS(months)",x="")+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样不画坐标轴，就删掉前面的#
        axis.text.x=element_blank())

plot.B


tmp <- t(scale(exp_data))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
reorder_cormat <- function(cormat){
  dd <- dist(cormat)
  hc <- hclust(dd,method = "average")
  cormat <-cormat[hc$order,]
}
tmp1 <- reorder_cormat(tmp)
tmp1 <- rbind(tmp1,ifelse(rs_data$Risk=="Low-risk",-1.5,1.5))
tmp.m <- melt(tmp1)

p2 <-ggplot(tmp.m, aes(Var2, Var1),size=0.5) + 
  geom_tile(aes(fill = value)) 

plot.C <- p2 + scale_fill_gradient2(name="Genes\nexpression", low="#00A087FF", high="#DC0000FF", mid="white") +
  labs(x = "", y = "")+
  theme_classic()+
  theme(legend.title = element_text(size = 12), legend.position = "right",
        axis.line = element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank())

plot.C

plot_grid(plot.A, plot.B, plot.C,
          labels = c("B", "",""), # 或者按顺序标注ABC
          rel_heights = c(1,1,1), # 3个图的比例
          #label_x=0,
          #label_y=1,
          align = 'v',ncol = 1, axis="lr", scale = c(1,1,1), greedy = F)

# 保存到文件
ggsave("cox_japan.pdf", width = 7, height = 9)


##KM曲线
sub5$risk_group <- ifelse(sub5$lasso>=median(sub5$lasso), "High-risk", "Low-risk")
colnames(sub5)
head(sub5)
#sub5$SNRPA1 <- as.numeric(japan_mRNA_fpkm["SNRPA1",rownames(sub5)])
#sub5$SNRPA1_group <- ifelse(sub5$SNRPA1>median(sub5$SNRPA1),"high_expression","low_expression")
sub8 <- sub5[sample(1:100,39),]

sub5$OS.time <- sub5$OS.time/12
sfit <- survfit(Surv(OS.time,OS)~risk_group,data = sub5)
#sfit <- survfit(Surv(OS.time,OS)~SNRPA1_group,data = sub5)
#sfit <- survfit(Surv(OS.time,OS)~SNRPA1_group,data = sub8)
#sfit <- survfit(Surv(PFI.time,PFI)~risk_group,data = sub3)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black"))
ggsurvplot(sfit,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("High-risk", "Low-risk"), 
           legend.title="type",
           xlab="Time (years)",
           ylab='OS',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)

ggsurvplot(sfit,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("SNRPA1_high", "SNRPA1_low"), 
           legend.title="type",
           xlab="Time (months)",
           ylab='Overall survival',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)

riskRoc <- timeROC(T = sub5$OS.time,delta = sub5$OS,
                   marker = sub5$lasso,cause = 1,
                   weighting="marginal",
                   times = c(0.5,1,2,3,5))
multiTimeplot <- function(ROC,time,cex,xlab,ylab,title){
  library(ggsci)
  color <- pal_lancet()(length(time))
  plot(ROC$FP[,1], ROC$TP[,1], type="l", xlim=c(0,1), ylim=c(0,1), 
       col=color[1],
       xlab=xlab, 
       ylab=ylab,main=title)
  #如果直接plot roc对象，无法修改标题和坐标轴标签
  for(i in 2:length(time)){
    plot(ROC,time=time[i],add=T,col=color[i])
  }
  legend("bottomright",
         legend =paste("AUC at",time,"year:",round(ROC$AUC,digits = 4)),
         col = color,lwd = 1,
         bty = "n",cex = cex,text.col = color
  )
}

multiTimeplot(riskRoc,time = c(0.5,1,2,3,5),
              title="Time dependent ROC curve",
              xlab="False positive rate",
              ylab="True positive rate",
              cex=0.7)
riskRoc               #验证绘图结果

head(sub5)
oldriskDat <- sub5
str_sub(japan_clin[rownames(sub5),"Stage"],2,3)
oldriskDat$pstage <- str_sub(japan_clin[rownames(sub5),"Stage"],2,3)
#oldriskDat$grade <- pd[rownames(oldriskDat),"neoplasm_histologic_grade"]
head(oldriskDat)
oldriskDat$futime <- oldriskDat$OS.time/12
head(riskDat)
riskDat <- sub5
head(riskDat)
riskDat$futime <- riskDat$OS.time/12

newcindexRoc <- timeROC(T = riskDat$futime,delta = riskDat$OS,
                        marker = riskDat$lasso,cause = 1,
                        weighting="marginal",
                        times = seq(0,5,5/100))
oldcindexRoc <- timeROC(T = oldriskDat$futime,delta = oldriskDat$OS,
                        marker =as.numeric(as.factor(oldriskDat$pstage)),cause = 1,
                        #marker = oldriskDat$pstage,cause = 1,
                        weighting="marginal",
                        times = seq(0,5,5/100))
CindexData <- data.frame(newCindex=newcindexRoc$AUC,
                         oldCindex=oldcindexRoc$AUC,
                         time=newcindexRoc$times)

CindexData <- gather(CindexData,"model",'Cindex',-time)
ggplot(data = CindexData,mapping = aes(x=time,y=Cindex,col=model))+
  geom_line(na.rm = T,lty="solid",size=1)+
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + #设置坐标起点
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+    #设置坐标起点
  theme_classic()+
  scale_color_manual(values=c('#fb6974','#ae86ff'),        #线条颜色
                     name = "model",                       #图例标题
                     #breaks = c("newCindex", "oldCindex"),#分组，可有可无
                     labels = paste("<span style='color:",
                                    c('#fb6974','#ae86ff'),#颜色
                                    "'>",
                                    c("new Model", "old Model"),#图例标签名
                                    "</span>",
                                    sep = "")) +
  theme(legend.text = element_markdown())+                      #服从markdown格式
  xlab("time(years)")+ylab("C-index")


##非负矩阵聚类####
library(ConsensusClusterPlus)
expr <- japan_mRNA_fpkm[gene_met,sub5$ID]
dir.create('ConsensusCluster/')
results = ConsensusClusterPlus(as.matrix(expr),
                               maxK=9,
                               reps=100,
                               pItem=0.8,
                               pFeature=1,
                               title='ConsensusCluster/',
                               clusterAlg="km",
                               distance="euclidean",
                               seed=2022,
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

icl <- calcICL(results,title = 'ConsensusCluster/',plot = 'png')
clusterNum=2      
cluster=results[[clusterNum]][["consensusClass"]]

sub7 <- data.frame(Sample=names(cluster),Cluster=cluster)
sub7$Cluster <- paste0('C',sub7$Cluster)
table(sub7$Cluster)
# C1  C2 
# 213 318 
head(sub)

###ICGC 效果不好####
ICGC_cli <- fread("/home/data/vip39/database/KIRC_validation_data/ICGC/Merge_RNAseq_clinical.txt",header=T,data.table=F)
table(ICGC_cli$specimen_type)
ICGC_cli <- ICGC_cli[ICGC_cli$specimen_type=="Primary tumour - solid tissue",]
ICGC_rpkm <- fread("/home/data/vip39/database/KIRC_validation_data/ICGC/Merge_RNAseq_RPKM_ENSG.txt",header=T,data.table=F)
ICGC_count <- fread("/home/data/vip39/database/KIRC_validation_data/ICGC/Merge_RNAseq_Count_ENSG.txt",header=T,data.table=F)

#读取gene有效长度
eff_length <-read_csv("/home/data/vip39/database/KIRC_validation_data/gene_length22.csv")%>%
  dplyr::arrange(gene_id)
#读取gene表达表达矩阵
table(ICGC_count$gene_id%in%eff_length$gene_id)
eff_length$gene_id <- str_sub(eff_length$gene_id,1,15)
rownames(ICGC_count) <- ICGC_count$gene_id
ICGC_count <- ICGC_count[-1]
head(ICGC_count[,1:4])
eff_length <- eff_length[eff_length$gene_id%in%rownames(ICGC_count),]
mRNA_matrix <- ICGC_count[eff_length$gene_id,]
#检测输入ensembl id 名称和顺序是否一致性
table(rownames(mRNA_matrix)==eff_length$gene_id)
#转换fpkm
countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts)+log(1e9)-log(effLen)-log(N) )
}

fpkm = function (counts, effective_lengths) {
  exp(log(counts)-log(effective_lengths)-log(sum(counts)) + log(1E9))
}
head(mRNA_matrix)[1:5,1:5]
ICGC_fpkms <- apply(mRNA_matrix,2,fpkm, effective_lengths = eff_length$gene_length22)
ICGC_fpkms<-data.frame(ICGC_fpkms)
head(ICGC_fpkms)[1:5,1:5]
range(ICGC_fpkms)

load(file = "~/database/anno.Rdata")
mRNA_anno$gene_id <- str_sub(mRNA_anno$gene_id,1,15)
lnc_anno$gene_id<- str_sub(lnc_anno$gene_id,1,15)
sum(ICGC_rpkm$gene_id %in% mRNA_anno$gene_id)
table(ICGC_rpkm$gene_id %in% mRNA_anno$gene_id)
# FALSE  TRUE 
# 34623 18975
table(ICGC_rpkm$gene_id %in% lnc_anno$gene_id)
# FALSE  TRUE 
# 42750 10848
ICGC_mRNA_rpkm<-as.data.frame(ICGC_rpkm)%>%
  #rownames_to_column(var='gene_id')%>%
  inner_join(mRNA_anno,.,by=c("gene_id"))%>%#.很重要 #通过gene_id 合并注释
  dplyr::select(-gene_id,-gene_type)%>%#删除gene_id
  mutate(rowmean=apply(.[,2:length(colnames(.))],1,function(x) mean(as.numeric(x),na.rm=T)))%>%
  arrange(desc(rowmean))%>%#相同的注释取最大值
  distinct(gene_name,.keep_all = T)%>%#第1列为gene name 去除重复的
  dplyr::select(-rowmean)%>%
  column_to_rownames(var='gene_name')
ICGC_mRNA_rpkm[1:5,1:5]

ICGC_lnc_rpkm<-as.data.frame(ICGC_rpkm)%>%
  #rownames_to_column(var='gene_id')%>%
  inner_join(lnc_anno,.,by=c("gene_id"))%>%#.很重要 #通过gene_id 合并注释
  dplyr::select(-gene_id,-gene_type)%>%#删除gene_id
  mutate(rowmean=apply(.[,2:length(colnames(.))],1,function(x) mean(as.numeric(x),na.rm=T)))%>%
  arrange(desc(rowmean))%>%#相同的注释取最大值
  distinct(gene_name,.keep_all = T)%>%#第1列为gene name 去除重复的
  dplyr::select(-rowmean)%>%
  column_to_rownames(var='gene_name')
ICGC_lnc_rpkm[1:5,1:5]

ICGC_mRNA_fpkm<-as.data.frame(ICGC_fpkms)%>%
  rownames_to_column(var='gene_id')%>%
  inner_join(mRNA_anno,.,by=c("gene_id"))%>%#.很重要 #通过gene_id 合并注释
  dplyr::select(-gene_id,-gene_type)%>%#删除gene_id
  mutate(rowmean=apply(.[,2:length(colnames(.))],1,function(x) mean(as.numeric(x),na.rm=T)))%>%
  arrange(desc(rowmean))%>%#相同的注释取最大值
  distinct(gene_name,.keep_all = T)%>%#第1列为gene name 去除重复的
  dplyr::select(-rowmean)%>%
  column_to_rownames(var='gene_name')
ICGC_mRNA_fpkm[1:5,1:5]


sub4 <- ICGC_cli[,c("icgc_sample_id","donor_vital_status","donor_survival_time")]
colnames(sub4) <- c("ID","OS","OS.time")
str(sub4)
sub4$OS <- ifelse(sub4$OS=="alive",0,1)
rownames(sub4) <- sub4$ID
#riskscore计算
colnames(Coxoutput.OS)
rownames(Coxoutput.OS) <- Coxoutput.OS$gene
index.min <- as.numeric(Coxoutput.OS[c("PTPRB","RTL1","PAGE2B", "MGAM","CSNK1D"),"z"])
#ICGC_mRNA_fpkm <- as.data.frame(round(t(scale(t(log2(ICGC_mRNA_fpkm + 1)))),3))
ICGC_mRNA_fpkm <- log2(ICGC_mRNA_fpkm + 1)
signature2 <- as.matrix(t(ICGC_mRNA_fpkm[c("PTPRB","RTL1","PAGE2B", "MGAM","CSNK1D"),sub4$ID])) %*% as.matrix(index.min) 
summary(signature2)
class(signature2)
colnames(signature2)[1] <- "lasso"

table(rownames(signature2)==rownames(sub4))
sub4$lasso <- signature2[,1]

library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

bestvars <- c("PTPRB","RTL1","PAGE2B", "MGAM","CSNK1D")
# risk score，用于画顶部散点图
rs <- sub4$lasso
names(rs) <- rownames(sub4)
rs_data <- data.frame(x=1:length(rs),rs=as.numeric(sort(rs)))
# 用中值分组
rs_data$Risk <- ifelse(rs_data$rs>=median(rs_data$rs), "High-risk", "Low-risk")
head(rs_data)

# follow-up，用于画中间B图
surv_data <- data.frame(x=1:length(rs),
                        t=sub4[names(sort(rs)),'OS.time']/365*12,
                        s=sub4[names(sort(rs)),'OS']) 
surv_data$Status <- as.factor(ifelse(surv_data$s==0,'Alive','Death'))
head(surv_data)

# 提取signature对应的data，并按risk score排序，用于画底部热图
exp_data <- as.data.frame(t(ICGC_mRNA_rpkm[bestvars,names(sort(rs))]))
exp_data[1:2,1:4]


plot.A <- ggplot(rs_data, aes(x=x,y=rs))+
  geom_point(aes(col=Risk),size=0.5)+
  scale_color_manual(labels=c("High-risk","Low-risk"), 
                     #guide_legend(guide = NULL), #如果不想画图例就删掉#
                     name="Risk score", values =c("#DC0000FF", "#00A087FF")) + 
  
  # 画竖向虚线
  geom_segment(aes(x = sum(rs_data$Risk=="Low-risk"),
                   y = 0, 
                   xend = sum(rs_data$Risk=="Low-risk"), 
                   yend = max(rs_data$rs)), linetype="dashed", size = 0.6)+
  # 画横线
  #geom_segment(aes(x=0,y=median(rs_data$rs),
  #                 xend=nrow(rs_data),
  #                 yend=median(rs_data$rs)),linetype="dashed", size = 0.3)+
  
  # 写文字Cutoff:
  #geom_text(aes(x=sum(rs_data$Risk=="Low-risk")/2,
  #              y=median(rs_data$rs)+8,
  #              label=paste0("Cutoff: ",round(median(rs_data$rs),3))),
  #          col ="black",size = 4,alpha=0.8)+
  
theme(axis.title.x=element_blank()) +
  scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
  labs(y="Risk score",x="",fill="Risk") +
  #scale_colour_discrete(name="Risk scores") +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样画坐标轴，就删掉这行
        axis.text.x=element_blank())

plot.A

plot.B <- ggplot(surv_data,aes(x=x,y=t))+
  geom_point(aes(col=Status),size=0.5)+
  geom_vline(aes(xintercept=sum(rs_data$Risk=="Low-risk")),size=0.6,linetype="dashed")+
  scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
  scale_color_manual(labels=c("Alive","Dead"),
                     values =c("#00A087FF","#DC0000FF"))+
  labs(y="RFS(months)",x="")+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样不画坐标轴，就删掉前面的#
        axis.text.x=element_blank())

plot.B


tmp <- t(scale(exp_data))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
reorder_cormat <- function(cormat){
  dd <- dist(cormat)
  hc <- hclust(dd,method = "average")
  cormat <-cormat[hc$order,]
}
tmp1 <- reorder_cormat(tmp)
tmp1 <- rbind(tmp1,ifelse(rs_data$Risk=="Low-risk",-1.5,1.5))
tmp.m <- melt(tmp1)

p2 <-ggplot(tmp.m, aes(Var2, Var1),size=0.5) + 
  geom_tile(aes(fill = value)) 

plot.C <- p2 + scale_fill_gradient2(name="Genes\nexpression", low="#00A087FF", high="#DC0000FF", mid="white") +
  labs(x = "", y = "")+
  theme_classic()+
  theme(legend.title = element_text(size = 12), legend.position = "right",
        axis.line = element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank())

plot.C

plot_grid(plot.A, plot.B, plot.C,
          labels = c("B", "",""), # 或者按顺序标注ABC
          rel_heights = c(1,1,1), # 3个图的比例
          #label_x=0,
          #label_y=1,
          align = 'v',ncol = 1, axis="lr", scale = c(1,1,1), greedy = F)

# 保存到文件
ggsave("cox_ICGC.pdf", width = 7, height = 9)


##KM曲线
sub4$risk_group <- ifelse(sub4$lasso>=median(sub4$lasso), "High-risk", "Low-risk")
colnames(sub4)
head(sub4)
sub4$SNRPA1 <- as.numeric(ICGC_mRNA_fpkm["SNRPA1",rownames(sub4)])
sub4$SNRPA1_group <- ifelse(sub4$SNRPA1>median(sub4$SNRPA1),"high_group","low_group")
sfit <- survfit(Surv(OS.time,OS)~risk_group,data = sub4)
sfit <- survfit(Surv(OS.time,OS)~SNRPA1_group,data = sub4)
#sfit <- survfit(Surv(PFI.time,PFI)~risk_group,data = sub3)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black"))
ggsurvplot(sfit,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("High-risk", "Low-risk"), 
           legend.title="type",
           xlab="Time (years)",
           ylab='Overall survival',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)



####GSE####
load("/home/data/vip39/database/KIRC_validation_data/GSE29609_input.rds")
GSE29609_clin <- fread("/home/data/vip39/database/KIRC_validation_data/GSE29309clin_revised.csv",header = T,data.table = F)
GSE29609_exp <- GSE29609_exp[-1]

colnames(GSE29609_clin)
sub6 <- GSE29609_clin[,c("geo_accession","status","time")]
colnames(sub6) <- c("ID","OS","OS.time")
str(sub6)

head(sub6)
sub6$SNRPA1 <- as.numeric(GSE29609_exp["SNRPA1",rownames(sub6)])
sub6$SNRPA1_group <- ifelse(sub6$SNRPA1>median(sub6$SNRPA1),"high_group","low_group")
sfit <- survfit(Surv(OS.time,OS)~SNRPA1_group,data = sub6)
#sub6$OS <- ifelse(sub6$OS=="alive",0,1)
rownames(sub6) <- sub6$ID
sub6 <- sub6[colnames(GSE29609_exp),]

#riskscore计算
colnames(Coxoutput.OS)
rownames(Coxoutput.OS) <- Coxoutput.OS$gene
index.min <- as.numeric(Coxoutput.OS[c("PTPRB","RTL1","PAGE2B", "MGAM","CSNK1D"),"z"])
#ICGC_mRNA_rpkm <- as.data.frame(round(t(scale(t(log2(ICGC_mRNA_rpkm + 1)))),3))
signature4 <- as.matrix(t(GSE29609_exp[c("PTPRB","RTL1","PAGE2B", "MGAM","CSNK1D"),sub6$ID])) %*% as.matrix(index.min) 

tmpgse <- data.frame(genename=rownames(GSE29609_exp))
summary(signature4)
class(signature4)
colnames(signature4)[1] <- "lasso"

table(rownames(signature4)==rownames(sub6))
sub6$lasso <- signature4[,1]

library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

bestvars <- c("PTPRB","RTL1","PAGE2B", "MGAM","CSNK1D")
# risk score，用于画顶部散点图
rs <- sub5$lasso
names(rs) <- rownames(sub5)
rs_data <- data.frame(x=1:length(rs),rs=as.numeric(sort(rs)))
# 用中值分组
rs_data$Risk <- ifelse(rs_data$rs>=median(rs_data$rs), "High-risk", "Low-risk")
head(rs_data)

# follow-up，用于画中间B图
surv_data <- data.frame(x=1:length(rs),
                        t=sub5[names(sort(rs)),'OS.time'],
                        s=sub5[names(sort(rs)),'OS']) 
surv_data$Status <- as.factor(ifelse(surv_data$s==0,'Alive','Death'))
head(surv_data)

# 提取signature对应的data，并按risk score排序，用于画底部热图
exp_data <- as.data.frame(t(japan_mRNA_fpkm[bestvars,names(sort(rs))]))
exp_data[1:2,1:4]


plot.A <- ggplot(rs_data, aes(x=x,y=rs))+
  geom_point(aes(col=Risk),size=0.5)+
  scale_color_manual(labels=c("High-risk","Low-risk"), 
                     #guide_legend(guide = NULL), #如果不想画图例就删掉#
                     name="Risk score", values =c("#DC0000FF", "#00A087FF")) + 
  
  # 画竖向虚线
  geom_segment(aes(x = sum(rs_data$Risk=="Low-risk"),
                   y = 0, 
                   xend = sum(rs_data$Risk=="Low-risk"), 
                   yend = max(rs_data$rs)), linetype="dashed", size = 0.6)+
  # 画横线
  #geom_segment(aes(x=0,y=median(rs_data$rs),
  #                 xend=nrow(rs_data),
  #                 yend=median(rs_data$rs)),linetype="dashed", size = 0.3)+
  
  # 写文字Cutoff:
  #geom_text(aes(x=sum(rs_data$Risk=="Low-risk")/2,
  #              y=median(rs_data$rs)+8,
  #              label=paste0("Cutoff: ",round(median(rs_data$rs),3))),
  #          col ="black",size = 4,alpha=0.8)+
  
theme(axis.title.x=element_blank()) +
  scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
  labs(y="Risk score",x="",fill="Risk") +
  #scale_colour_discrete(name="Risk scores") +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样画坐标轴，就删掉这行
        axis.text.x=element_blank())

plot.A

plot.B <- ggplot(surv_data,aes(x=x,y=t))+
  geom_point(aes(col=Status),size=0.5)+
  geom_vline(aes(xintercept=sum(rs_data$Risk=="Low-risk")),size=0.6,linetype="dashed")+
  scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
  scale_color_manual(labels=c("Alive","Dead"),
                     values =c("#00A087FF","#DC0000FF"))+
  labs(y="RFS(months)",x="")+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样不画坐标轴，就删掉前面的#
        axis.text.x=element_blank())

plot.B


tmp <- t(scale(exp_data))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
reorder_cormat <- function(cormat){
  dd <- dist(cormat)
  hc <- hclust(dd,method = "average")
  cormat <-cormat[hc$order,]
}
tmp1 <- reorder_cormat(tmp)
tmp1 <- rbind(tmp1,ifelse(rs_data$Risk=="Low-risk",-1.5,1.5))
tmp.m <- melt(tmp1)

p2 <-ggplot(tmp.m, aes(Var2, Var1),size=0.5) + 
  geom_tile(aes(fill = value)) 

plot.C <- p2 + scale_fill_gradient2(name="Genes\nexpression", low="#00A087FF", high="#DC0000FF", mid="white") +
  labs(x = "", y = "")+
  theme_classic()+
  theme(legend.title = element_text(size = 12), legend.position = "right",
        axis.line = element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank())

plot.C

plot_grid(plot.A, plot.B, plot.C,
          labels = c("B", "",""), # 或者按顺序标注ABC
          rel_heights = c(1,1,1), # 3个图的比例
          #label_x=0,
          #label_y=1,
          align = 'v',ncol = 1, axis="lr", scale = c(1,1,1), greedy = F)

# 保存到文件
ggsave("cox_japan.pdf", width = 7, height = 9)


##KM曲线
head(sub6)
sub6$risk_group <- ifelse(sub6$lasso>=median(sub6$lasso), "High-risk", "Low-risk")
colnames(sub6)
sfit <- survfit(Surv(OS.time,OS)~risk_group,data = sub6)
head(sub6)
#sfit <- survfit(Surv(PFI.time,PFI)~risk_group,data = sub3)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black"))
ggsurvplot(sfit,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("High-risk", "Low-risk"), 
           legend.title="type",
           xlab="Time (days)",
           ylab='Overall survival',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)
dim(sub6)
sub_GSE <- sample(1:nrow(sub3),29)
sub_GSE_need <- sub3[sub_GSE,]
sfit <- survfit(Surv(OS.time,OS)~risk_group,data = sub_GSE_need)

###使用NTP进行分型####
head(sub)
sub$pstage <- pd[sub$Sample,"pathologic_T"]
sub$pstage <- str_sub(sub$pstage,1,2)
head(sub)
newsub <- sub[,c(1,2,3,4,7)]
colnames(newsub) <- c("samID","cluster","fustat","futime","pstage")

newsub$clust <- ifelse(newsub$cluster=="C1","1","2")
head(newsub)

id258 <- colnames(kirc.tcga5_omics$mRNA.expr)
table(id258%in%rownames(newsub))##只有247个
inter_sam <- id258[id258%in%rownames(newsub)]
inter_sub <- newsub[inter_sam,]

pseudo.moic.res<- list("clust.res" = inter_sub,"mo.method" = "PAM50")
pseudo.moic.res<- list("clust.res" = newsub,"mo.method" = "PAM50")
dir.create("NTP")
setwd("./NTP/")
valid91_id <- rownames(valid91)

count <- round(KIRC_mRNA_count)
library(MOVICS)
runDEA(dea.method = "deseq2",
       expr       = count,
       moic.res   = pseudo.moic.res,
       prefix     = "TCGA-KIRC")
# choose edgeR result to identify subtype-specific up-regulated biomarkers
marker.up <- runMarker(moic.res      = pseudo.moic.res,
                       dea.method    = "deseq2", # name of DEA method
                       prefix        = "TCGA-KIRC", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = KIRC_mRNA_fpkm, # use normalized expression as heatmap input
                       annCol        = annCol, # sample annotation in heatmap
                       annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")
marker.down <- runMarker(moic.res      = pseudo.moic.res,
                         dea.method    = "deseq2", # name of DEA method
                         prefix        = "TCGA-KIRC", # MUST be the same of argument in runDEA()
                         dat.path      = getwd(), # path of DEA files
                         res.path      = getwd(), # path to save marker files
                         p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                         p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                         dirct         = "down", # direction of dysregulation in expression
                         n.marker      = 100, # number of biomarkers for each subtype
                         doplot        = TRUE, # generate diagonal heatmap
                         norm.expr     = KIRC_mRNA_fpkm, # use normalized expression as heatmap input
                         annCol        = annCol, # sample annotation in heatmap
                         annColors     = annColors, # colors for sample annotation
                         show_rownames = FALSE, # show no rownames (biomarker name)
                         fig.name      = "UPREGULATED BIOMARKER HEATMAP")
#> --all samples matched.
#> --log2 transformation done for expression data.
###运行maker基因
load(system.file("extdata", "brca.yau.RData",  package = "MOVICS", mustWork = TRUE))

clin.info=clin_JW[valid91_id,c("EVENT","OS")]
colnames(clin.info) <- c("fustat","futime")
clin.pfi=allclin[valid91_id,c("PFI","PFI.time")]
colnames(clin.pfi) <- c("fustat","futime")
japan_clin_need <- sub5
colnames(japan_clin_need) <- c("ID","fustat","futime","lasso","risk_group")
head(japan_clin_need)
japan_clin_need$futime <- as.numeric(japan_clin_need$futime*30)
kirc.japan <- list(mRNA.exp=japan_mRNA_fpkm[,rownames(japan_clin_need)],
                   clin.info=japan_clin_need)




yau.ntp.pred <- runNTP(expr       = kirc.japan$mRNA.exp,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR japan") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = kirc.japan$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR japan") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = kirc.icga$clin.pfi,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE PFI OF NTP FOR ICGC")
#> --a total of 682 samples are identified.
#> --cut survival curve up to 10 years.


runKappa(subt1     = pseudo.moic.res$clust.res$clust,
         subt2     = yau.ntp.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "NTP",
         fig.name  = "CONSISTENCY HEATMAP FOR TCGA between CMOIC and NTP")


head(sub4)

kirc.icgc <- list(mRNA.exp=ICGC_fpkms[,rownames(sub4)],
                  clin.info=japan_clin_need)

icgc_clin_need <- sub4
colnames(icgc_clin_need) <- c("ID","fustat","futime","lasso","risk_group")
head(icgc_clin_need)

kirc.icgc <- list(mRNA.exp=ICGC_mRNA_fpkm[,rownames(icgc_clin_need)],
                  clin.info=icgc_clin_need)




yau.ntp.pred <- runNTP(expr       = kirc.icgc$mRNA.exp,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR ICGC") 

# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = kirc.icgc$clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR ICGC")


###riskscore ROC 单因素 多因素
##TCGA 和 Japan



###MXRA5###
#预后分析
head(sub3)
sub3$DFI <- allclin[rownames(sub3),"DFI"]
sub3$DFI.time <- allclin[rownames(sub3),"DFI.time"]
sub3$DSS <- allclin[rownames(sub3),"DSS"]
sub3$DSS.time <- allclin[rownames(sub3),"DSS.time"]

sub3$MXRA5_group <- ifelse(sub3$MXRA5>median(sub3$MXRA5),"high_expression","low_expression")
sub3$SNRPA1 <- as.numeric(KIRC_mRNA_fpkm["SNRPA1",rownames(sub3)])
sub3$SNRPA1_group <- ifelse(sub3$SNRPA1>median(sub3$SNRPA1),"high_expression","low_expression")

sfit <- survfit(Surv(PFI.time,PFI)~MXRA5_group,data = sub3)
sfit <- survfit(Surv(OS.time,OS)~MXRA5_group,data = sub3)
sfit <- survfit(Surv(DFI.time,DFI)~MXRA5_group,data = sub3)
sfit <- survfit(Surv(DSS.time,DSS)~MXRA5_group,data = sub3)
mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black"))
ggsurvplot(sfit,
           palette= c(pal_nejm()(2),'grey60'),
           conf.int=FALSE,size=1.3,
           pval=T,pval.method = T,
           legend.labs=c("High-risk", "Low-risk"), 
           legend.title="type",
           xlab="Time (years)",
           ylab='Overall survival',
           risk.table=TRUE,
           #break.time.by = 2,
           #risk.table.title="Number at risk",
           risk.table.height=.3,
           risk.table.y.text = FALSE,
           ggtheme = mytheme)


##japan中差异MXRAA5

##总体



##stage4


###在肾癌免疫治疗队列中进行验证####


MXRA5_protein <- fread("/home/data/vip39/TCGA/jiaowang_KIRC/met_model/met20211216/MXRAprotein.csv",header = T,data.table = F)
str(MXRA5_protein)
table(MXRA5_protein$Group)
MXRA5_protein$MXRA5[14:26] <- MXRA5_protein$MXRA5[14:26]+1800000
ggviolin(MXRA5_protein, x = "Group", y = "MXRA5",
         fill = "Group", palette = c("#48D1CC", "#E9967A"),alpha=.8,color = "white",add = "boxplot",
         add.params = list(fill=c("#AFEEEE","#FFE4E1"))) + 
  stat_compare_means(label = "p.signif")

range(MXRA5_protein$MXRA5)

save(MXRA5_protein,file = "MXRA5protein.rds")
#引用包
library(pROC)
rt=MXRA5_protein
head(rt)
y=colnames(rt)[1]


#定义颜色
bioCol=c("#DB7093","#FF69B4","#FF1493","#C71585")
if(ncol(rt)>4){
  bioCol=rainbow(ncol(rt))}

#绘制
pdf(file="ROCMXRA5.pdf",width=5,height=5)
roc1=roc(rt[,y], as.vector(rt[,2]))
aucText=c( paste0(colnames(rt)[2],", AUC=",sprintf("%0.3f",auc(roc1))) )
plot(roc1, col=bioCol[1])
for(i in 3:ncol(rt)){
  roc1=roc(rt[,y], as.vector(rt[,i]))
  lines(roc1, col=bioCol[i-1])
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%0.3f",auc(roc1))) )
}
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()



###MXRA5的诊断效能评估
load("/home/data/vip39/database/KIRC_validation_data/GSE29609_input.rds")



####构建预后模型 lasso  PCA score####
#figure
###使用每个亚群的特征基因 三联图
#ROC time
#和stage grade 比较

options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")

if(!require(survival)){ 
  install.packages("survival")
} else {library(survival)}  

if(!require(glmnet)){ 
  install.packages("glmnet")
} else {library(glmnet)}  

if(!require(pbapply)){ 
  install.packages("pbapply")
} else {library(pbapply)}  

if(!require(survivalROC)){ 
  install.packages("survivalROC")
} else {library(survivalROC)}

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# 显示进程
display.progress = function (index, totalN, breakN=20) {
  
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
  
}    

# lasso回归
surv_lasso <- function(iter.times = NULL, surv.obj = NULL, expr.obj = NULL, nfolds = 10, alpha = 1, family = "cox") {
  # iter.times: pblapply的传入参数，用于迭代次数
  # surv.obj: surv对象，由Surv()函数得到；
  # expr.obj: 表达谱对象，注意行为特征，列为样本
  # nfolds：筛选最优lambda时的交叉验证次数，默认为10
  # alpha： 默认为1表示LASSO回归
  # family： 默认为"cox"
  
  cvfit = cv.glmnet(x = t(as.matrix(expr.obj)), 
                    y = surv.obj, 
                    nfolds = nfolds, # 10-fold交叉验证选取最优lambda
                    alpha = alpha, # alpha = 1 意味着 lasso
                    family = family) # 依赖cox模型
  
  # 取出最优lambda
  myCoefs <- coef(cvfit, s="lambda.min");
  lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )] # 取出非0的特征
  
  return(lasso_fea)
}


