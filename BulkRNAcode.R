###################################Bulk-RNA part##############################
library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
library(GSEABase)
library(GSVA)
rt=read.table("xxxx.txt",header=T,check.names=F,sep="\t")
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)

max(rt)
if(max(rt)>30) rt=log2(rt+1)    

data=rt
afcon=4
conData=data[,as.vector(colnames(data)[1:afcon])]
aftreat=afcon+1
treatData=data[,as.vector(colnames(data)[aftreat:ncol(data)])]
rt=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

Diff=topTable(fit2,adjust='fdr',number=length(rownames(data)))
DIFFOUT=rbind(id=colnames(Diff),Diff)
write.table(DIFFOUT,file="DIFF_alllog.xls",sep="\t",quote=F,col.names=F)
diffSig=Diff[with(Diff, (abs(logFC) < -1 & P.Value < 0.05 )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="DIFF_af down.xls",sep="\t",quote=F,col.names=F)
diffSig=Diff[with(Diff, (abs(logFC) > 1 & P.Value < 0.05 )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="DIFF_af up.xls",sep="\t",quote=F,col.names=F)



adjP <- 0.05
aflogFC <- 1
Significant <- ifelse((Diff$P.Value < adjP & abs(Diff$logFC) > aflogFC), ifelse(Diff$logFC > aflogFC, "Up", "Down"), "Not")

# Start plotting
p <- ggplot(Diff, aes(logFC, -log10(P.Value))) +
  geom_point(aes(col = Significant), size = 3) +
  scale_color_manual(values = c(pal_npg()(2)[2], "#838B8B", pal_npg()(1))) +
  labs(title = " ") +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) +
  geom_hline(aes(yintercept = -log10(adjP)), colour = "gray", linetype = "twodash", linewidth = 1) +
  geom_vline(aes(xintercept = aflogFC), colour = "gray", linetype = "twodash", linewidth = 1) +
  geom_vline(aes(xintercept = -aflogFC), colour = "gray", linetype = "twodash", linewidth = 1)

# Adding label
point.Pvalue <- 0.01
point.logFc <- 300

Diff$symbol <- rownames(Diff)

# Filter for points to label
for_label <- Diff %>% 
  filter(abs(logFC) > point.logFc & P.Value < point.Pvalue)

# Filter for SLC7A11
slc7a11_label <- Diff %>%
  filter(symbol %in% c("SLC7A11"))

# Plot with additional labeling
pdf("DIFF_vol8086.pdf", width = 6, height = 5)
p <- p + theme_bw() + theme(panel.grid = element_blank())

p + geom_point(size = 1.5, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color = "black",
    label.size = 0.1
  ) +
  geom_point(data = slc7a11_label, aes(logFC, -log10(P.Value)), color = "green", size = 3) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = slc7a11_label,
    color = "red",
    label.size = 0.1
  )
dev.off()


source('CIBERSORT.R')
result1 <- CIBERSORT('LM22V2.txt','allcountstxt', perm = 1000, QN = T) 
library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)
my36colors <-c('#E95C59','#E5D2DD', '#5F3D69', '#F1BB72', '#F3B1A0', '#57C3F3', '#476D87', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#B53E2B','#53A85F', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', 
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)
rt=read.table(("CIBERSORT-Results.txt"),sep="\t",header=T,row.names = 1)


cibersort=as.data.frame(rt)
immune=cibersort
data=as.matrix(immune[,0:(ncol(immune)-3)])
colnames(immune)=gsub("CIBERSORT-Results"," ",colnames(immune))
colnames(data)=gsub("CIBERSORT-Results"," ",colnames(data))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data)=gsub("_"," ",colnames(data))
cluster=read.table("sample.txt", header=F, sep="\t", check.names=F, row.names=1)
colnames(cluster)="Type"
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$Type),]
gaps=c(1, as.vector(cumsum(table(data$Type))))
xlabels=levels(factor(data$Type))
data$Type=factor(data$Type, levels=c(cluster[1,1],cluster[nrow(cluster),1]))
data$GSM=rownames(data)
data=melt(data,id.vars=c("Type","GSM"))
colnames(data)=c("Type","GSM","Celltype", "Freq")
Cellratio=data
colourCount = length(unique(Cellratio$Celltype))
ggplot(Cellratio) + 
  geom_bar(aes(x =GSM, y= Freq, fill = Celltype),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Cell cycle phase',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"),
        legend.position = "right" )+   
  scale_fill_manual(values=my36colors)+
  theme_bw()+
  xlab(NULL)+
  theme(axis.text.x  = element_blank())+
  guides(
    fill = guide_legend(
      ncol = 2, 
      byrow = TRUE, 
      label.theme = element_text(size = 15),  # 
    ))+
  facet_grid(. ~ Type,scales="free")+
  theme(strip.text.x = element_text(size = 30,colour = "black"))+       
  theme(strip.background.x = element_rect(fill = c("#C1E6F3"), colour = "black"))+
  theme(
    legend.key.size = unit(33, "pt"),   # 
    legend.title = element_text(size = 20),  #
    legend.box.background = element_rect(color = "black", size = 3)  # 
    
  )
ggsave("immune.ration.pdf",width = 15,height = 8)        




###KM analysis
library(GEOquery)
library(limma)
library(affy)
library(survival)
library(survminer)

a1 <- read.table("xxxx.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
a2 <- a1[,c('futime',"fustat","SLC7A11")]
names(a2)[names(a2) == "futime"] <- "FUTIME"
names(a2)[names(a2) == "fustat"] <- "FUSTAT"


a2$FUTIME <-a2$FUTIME / 30
merged_data = a2
merged_data$SLC7A11 <- as.numeric(as.character(merged_data$SLC7A11))
merged_data <- merged_data[merged_data$FUTIME >= 1, ]

cutoff<-surv_cutpoint(merged_data, 
                      time="FUTIME",
                      event="FUSTAT",
                      variables="SLC7A11");summary(cutoff) #

###中位值
median_SLC7A11 <- median(merged_data$SLC7A11, na.rm = TRUE)
merged_data$group <- ifelse(merged_data$SLC7A11 >= median_SLC7A11, "High", "Low")
str(merged_data)
head(merged_data)
fit <- survfit(Surv(FUTIME, FUSTAT) ~ group, data = merged_data)

ggsurvplot(fit, 
           data = merged_data,     
           pval=TRUE,         
           pval.method=TRUE,  
           palette = c( "red", "#00A087FF"),
           risk.table = T, #显示风险表
           conf.int = F, 
           title = "ICGC Cohort",
           font.main = c(16, "bold", "darkblue"),
           font.x = c(16, "plain", "black"),
           font.y = c(16, "plain", "black"),
           font.tickslab = c(16, "plain", "black"),
           font.legend = c(16, "plain", "black"))   

####COX


cox_model <- coxph(Surv(FUTIME, FUSTAT) ~ SLC7A11, data = merged_data)
ggforest(cox_model, data = merged_data)




