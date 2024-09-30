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





###################################Singlecell part##############################
library(Seurat)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(tidyr)
library(viridis)
library(SeuratWrappers)
library(DoubletFinder)
library(Matrix)

my36colors <-c('#E95C59','#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)
my2colors <- c("#2ca02c","#1f77b4")
GSE140228.data <- readMM("./data/matrix.mtx") 
dim(GSE140228.data)#54574 66187
gene.names = read.delim('./data/features.tsv.gz',header = F, stringsAsFactors = FALSE)
dim(gene.names)#54574 7
gene.names[1:3,1:3]
barcode.names = read.delim("./data/barcodes.tsv.gz",header = FALSE,stringsAsFactors = FALSE)
dim(barcode.names) #66187     1
head(barcode.names)
colnames(GSE140228.data) = barcode.names$V1
rownames(GSE140228.data) = gene.names$V2 
GSE140228.data[1:4,1:4]
dim(GSE140228.data)
table(duplicated(GSE140228.data@Dimnames[[1]]))##
GSE140228.data=GSE140228.data[!duplicated(GSE140228.data@Dimnames[[1]]) ,]##
GSE140228_HCC <- CreateSeuratObject(GSE140228.data,min.cells=3,min.features = 40,project="GSE140228")
dim(GSE140228_HCC)#[1] 32285  9575

anno_data<-read.table('./data/cellinfo.tsv.gz',header = T)
#add metadata
anno_data <- anno_data[match(rownames(GSE140228_HCC@meta.data),anno_data$Barcode),]
colnames(anno_data)
GSE140228_HCC@meta.data[,4:12] <- anno_data[,1:9]
GSE140228_HCC$Sample <- paste(GSE140228_HCC$Donor,GSE140228_HCC$Tissue_sub,sep = '_')
a1 = GSE140228_HCC@meta.data[,"Donor"]
a1 = GSE140228_HCC@meta.data[,"Tissue_sub"]
class(GSE140228_HCC[['Sample']])
table(GSE140228_HCC[['Sample']])
unique(GSE140228_HCC$Sample)
class(GSE140228_HCC@meta.data[["Sample"]])
Idents(GSE140228_HCC) <- GSE140228_HCC@meta.data[["Sample"]]



###################################QC part######################################
#mitochondrial QC metrics
GSE140228_HCC[["percent.MT"]] <- PercentageFeatureSet(GSE140228_HCC, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(GSE140228_HCC@assays$RNA))
HB.genes <- rownames(GSE140228_HCC@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)] 
HB.genes
GSE140228_HCC[["percent.HB"]]<-PercentageFeatureSet(GSE140228_HCC, features=HB.genes) 
# Visualize QC metrics as a violin plot
VlnPlot(GSE140228_HCC, features = c("nFeature_RNA", "nCount_RNA", "percent.MT","percent.HB"), ncol = 4,pt.size = 0.00)
ggsave("vlnplot_before_qc.pdf",dpi = 300, width = 12, height = 6) 

##filter
allcell <- subset(GSE140228_HCC, subset = 
                    nFeature_RNA > 400 &
                    nCount_RNA >1000 &
                    nCount_RNA <20000 & 
                    percent.MT < 10 #& percent.Hb<0.1
)




Idents(allcell)<-allcell$Tissue
table(allcell$Tissue)
allcell = subset(allcell,idents = c('Normal','Tumor'))
Idents(allcell)<-allcell$Tissue
table(allcell$Tissue)

###################################Integrate data  part#########################
library(harmony)
DefaultAssay(allcell2) <- "RNA"  
allcell_harmony <- NormalizeData(allcell2) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
system.time({allcell_harmony <- RunHarmony(allcell_harmony, group.by.vars = "Sample")})
allcell_harmony <- RunUMAP(allcell_harmony, reduction = "harmony", dims = 1:30)
DimPlot(object = allcell_harmony , reduction = "umap", label = T,label.size = 6,repel = T,pt.size = 0.01)  + NoLegend()+NoAxes()#
allcell_harmony <- RunTSNE(allcell_harmony, reduction = "harmony", dims = 1:30)
DimPlot(object = allcell_harmony , reduction = "tsne", label = T,label.size = 6,repel = T,pt.size = 0.01
)  + NoLegend()+NoAxes()#分组umap


allcell_harmony <- FindNeighbors(allcell_harmony, reduction = "harmony", dims = 1:30)
allcell_harmony <- FindClusters(allcell_harmony,resolution = 1)
table(allcell_harmony$seurat_clusters)
Idents(allcell_harmony) <- allcell_harmony@meta.data$seurat_cluster
new.cluster.ids <- c(1:length(levels(allcell_harmony)))
names(x = new.cluster.ids) <- levels(x = allcell_harmony)
allcell_harmony <- RenameIdents(object = allcell_harmony, new.cluster.ids)
allcell_harmony[["Clusters"]]<-Idents(allcell_harmony)
p=DimPlot(allcell_harmony, reduction = "tsne", group.by = c("Tissue","seurat_clusters"))
p
ggsave("tsne_harmony_evidence.pdf",p,width = 12,height = 5)




################ CCI #####
allcell.list <- SplitObject(allcell2, split.by = 'Sample')
DefaultAssay(allcell2)<-"RNA"
for (i in 1:length(allcell.list)) {
  allcell.list[[i]] <- NormalizeData(allcell.list[[i]], verbose = FALSE)
  allcell.list[[i]] <- FindVariableFeatures(allcell.list[[i]], selection.method = "vst", 
                                            nfeatures = 2000, verbose = FALSE)
}

dims=1:30
allcell.anchors <- FindIntegrationAnchors(object.list = allcell.list, dims = dims)
allcell2 <- IntegrateData(anchorset = allcell.anchors, dims = dims)
allcell2 <- ScaleData(allcell2, verbose = FALSE)
allcell2 <- RunPCA(allcell2, npcs = 30, verbose = FALSE)
allcell2 <- FindNeighbors(allcell2, reduction = "pca", dims = dims)
allcell2 <- FindClusters(allcell2, resolution = 0.3)#22
head(allcell2@meta.data)

a=allcell2@meta.data[,-c(11:15)]
allcell2@meta.data=a
dim(allcell2)# 2000 29460
allcell3<-allcell_harmony@meta.data[,c(14:16)]
allcell3<-cbind(allcell3,allcell2@meta.data)
allcell2@meta.data<-allcell3

Idents(allcell2) <- allcell2@meta.data$Clusters
allcell2 <- RunTSNE(allcell2, dims =dims,perplexity=30)
#DimPlot(object = allcell2 , reduction = "tsne", label = T,label.size = 6,repel = T,pt.size = 0.01)  
#ggsave("tsne_allcell.pdf",width = 6,height = 5)
Idents(allcell_harmony)<-allcell_harmony$seurat_clusters
DimPlot(object = allcell_harmony , reduction = "tsne", label = T,label.size = 6,repel = T,pt.size = 0.01)  #分组umap
ggsave("tsne_harmony.pdf",width = 6,height = 5)
allcell = allcell2


###################################allcell part###############################

Idents(allcell_harmony)<-allcell_harmony$Clusters
Idents(allcell_harmony)<-allcell_harmony$Clusters
allcell_harmony.markers_res0.5 <- FindAllMarkers(allcell_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
allcell_harmony.markers_res0.5 <- allcell_harmony.markers_res0.5
allcell_harmony.markers_res0.5$nLog_padj <- -log10(allcell_harmony.markers_res0.5$p_val_adj) 
allcell_harmony.markers_res0.5$ratio <- allcell_harmony.markers_res0.5$pct.1/allcell_harmony.markers_res0.5$pct.2
out <- cbind(gene=allcell_harmony.markers_res0.5$gene,allcell_harmony.markers_res0.5)
write.table(out,'allcell_harmony_marker_res0.5.xls',sep = '\t',quote = F,row.names = F)
top10markers <- allcell_harmony.markers_res0.5 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10markers,"allcell_harmony_top10.csv")

Idents(allcell_harmony)<-allcell_harmony$Clusters
table(allcell_harmony$Clusters)

p=FeaturePlot(allcell_harmony,features=c('SLC7A11'), reduction = "tsne",
              pt.size = 1,raster=FALSE,order=T)
p
ggsave("SLC7A11+tsne.pdf",p,width = 6,height = 5)

######celltype name
allcell@meta.data$celltype <- as.character(allcell@meta.data$Clusters)
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(1))] <- 'CD4T' #
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(2))] <- 'CD8T' #
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(3))] <- 'NK/NKT'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(4))] <- 'CD8T'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(5))] <- 'NK/NKT'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(6))] <- 'CD8T'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(7))] <- 'Monocyte'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(8))] <- 'TAM'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(9))] <- 'CD4T'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(10))] <- 'DC'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(11))] <- 'CD8T'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(12))] <- 'B'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(13))] <- 'CD4T'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(14))] <- 'NK/NKT'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(15))] <- 'Monocyte'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(16))] <- 'B'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(17))] <- 'DC'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(18))] <- 'Unknown'
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(19))] <- 'DC'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(20))] <- 'Mast Cell'
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(21))] <- 'DC'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(22))] <- 'CD8T'#
allcell@meta.data$celltype[allcell@meta.data$Clusters %in% as.character(c(23))] <- 'Unknown'

allcell_harmony@meta.data$Clusters <- allcell@meta.data$celltype
allcell_harmony@meta.data$celltype <- allcell@meta.data$celltype

Idents(allcell_harmony) = allcell_harmony@meta.data$celltype 
unique(allcell_harmony$celltype)
table(allcell_harmony$Clusters)



###fliter unknown
Idents(allcell_harmony) <- allcell_harmony@meta.data$celltype
unique(allcell_harmony$celltype)
table(Idents(allcell_harmony))

allcell_harmony <- subset(allcell_harmony, idents =  c( "CD8T",
                                                        "CD4T",
                                                        "NK/NKT" ,
                                                        "TAM" ,
                                                        "B" ,
                                                        "DC",
                                                        'Mast Cell',
                                                        'Monocyte'
))

Idents(allcell) <- allcell@meta.data$celltype
unique(allcell$celltype)
table(Idents(allcell))

allcell <- subset(allcell, idents =  c("CD8T",
                                       "CD4T",
                                       "NK/NKT" ,
                                       "TAM" ,
                                       "B" ,
                                       "DC",
                                       'Mast Cell',
                                       'Monocyte'
))


DimPlot(
  object = allcell_harmony,
  reduction = "tsne",
  label = TRUE,
  label.size = 6,
  repel = TRUE,
  pt.size = 1,
  cols = my36colors
) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 16)  # 设置标题居中
    
  ) +
  NoLegend()+
  ggtitle('Celltype in scRNA')
ggsave("tsne_harmony.pdf",width =5.5,height = 5)


DimPlot(
  object = allcell_harmony,
  reduction = "tsne",
  label = TRUE,
  label.size = 6,
  repel = TRUE,
  pt.size = 1,
  cols = my2colors,
  group.by = "Tissue"
) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    legend.position = c(0.05, 0.05), # 
    legend.justification = c("left", "bottom"),
    legend.background = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 16)  # 
  )+ 
  ggtitle('Tissue type in scRNA')
ggsave("tsne_Tissue.pdf",width =5.5,height = 5)

DimPlot(
  object = allcell_harmony,
  reduction = "tsne",
  label = FALSE,
  label.size = 6,
  repel = TRUE,
  pt.size = 1,
  cols = my36colors,
  group.by = "Sample"
) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    legend.background = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 16)  # 
    
  )+ 
  ggtitle('Sample in scRNA')
ggsave("tsne_Sample.pdf",width =7.5,height = 5)

######marker present
library(RColorBrewer)
FeaturePlot(
  allcell_harmony,
  features = c('SLC7A11'),
  reduction = "tsne",
  pt.size = 1,
  raster = FALSE,
  order = T
) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black")+  
      coord_cartesian(xlim = c(-40, 40), ylim = c(-40, 40)) + 
      annotate("text", x = -20, y = -35, 
               label = "SLC7A11", hjust = 0, vjust = 0, 
               size = 6, colour = "black"))+
  ggtitle("SLC7A11 Expression")

ggsave("FeaturePlotSLC7A11.pdf",width = 6,height = 5)

FeaturePlot(
  allcell_harmony,
  features = c('xxx'),
  reduction = "tsne",
  pt.size = 0.5,
  raster = FALSE,
  order = T
) + 
  scale_color_gradientn(colors = brewer.pal(n = 9, name = "Reds")) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    plot.title = element_text(size = 20)  # 
  ) +   
  coord_cartesian(xlim = c(-40, 40), ylim = c(-40, 40)) + 
  ggtitle("xxx")

ggsave("FeaturePlotT.pdf",width = 6,height = 5)

allcell_harmony.markers_res0.5 <- FindAllMarkers(allcell_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
allcell_harmony.markers_res0.5 <- allcell_harmony.markers_res0.5
allcell_harmony.markers_res0.5$nLog_padj <- -log10(allcell_harmony.markers_res0.5$p_val_adj) 
allcell_harmony.markers_res0.5$ratio <- allcell_harmony.markers_res0.5$pct.1/allcell_harmony.markers_res0.5$pct.2
out <- cbind(gene=allcell_harmony.markers_res0.5$gene,allcell_harmony.markers_res0.5)
top10markers <- allcell_harmony.markers_res0.5 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10markers,"allcell_harmony_top10.csv")

genes_to_check =list(
  Mast_Cell=c('TPSAB1','TPSB2','CPA3'), 
  B =c('CD79A','MS4A1','CD19'),
  CD4T =c( 'FOXP3','TNFRSF18','TNFRSF4'), 
  TAM=c('CD68', 'C1QC','C1QB') , 
  Monocyte=c('VCAN', 'S100A8','S100A9'),
  CD8T =c( 'CD8A','CD8B','GZMK'),
  DC=c('CD1C', 'CD1E','CLEC10A'),
  NKNKT =c( 'KLRF1','GNLY','KLRC1')
  
)
p_all_markers=DotPlot(allcell_harmony, 
                      features = genes_to_check,
                      scale = T,assay='RNA')+
  theme(axis.text.x=element_text(angle=45,hjust = 1))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) 

p_all_markers
ggsave('harmonycheckmarkers.pdf',height = 6,width = 12)
FeaturePlot(allcell_harmony,features=c('SLC7A11'), reduction = "tsne",
            pt.size = 1,raster=FALSE,order=T)
ggsave("FeaturePlot_SLC7A11.pdf",width = 6,height = 5)


###################################SLC7A11+/- cells ratio  part##################
allcell_harmony$Ccl2='Negative'
k_posi<-subset(allcell_harmony,rna_SLC7A11>0)
allcell_harmony$Ccl2[match(rownames(k_posi@meta.data),rownames(allcell_harmony@meta.data))]='Positive'
VlnPlot(allcell_harmony,features ="rna_SLC7A11",group.by = 'Ccl2',pt.size = 0.001)#Positive
table(allcell_harmony$Ccl2)
ggsave("pos_SLC7A11.pdf",width = 6,height = 4)

VlnPlot(object = allcell_harmony, features = c('SLC7A11'), pt.size = 1, slot = "data")
ggsave("pos_SLC7A11.pdf",width = 6,height = 4)

####1、barplot####
colnames(allcell_harmony@meta.data)#Tissue_sub,slc7a11
tdata <- allcell_harmony@meta.data[,c('celltype','Tissue')]#
tdata1 <- table(tdata)
tdata2 <- as.data.frame.matrix(t(apply(tdata1,1,function(x){unlist(lapply(x,function(z){z/sum(x)}))})))
tdata2$group <- rownames(tdata2)
rownames(tdata2)
tdata3 <- gather(tdata2,cluster,percentage,-group)
tdata3$percentage <- tdata3$percentage*100

p <- ggplot(tdata3, aes(fill=cluster, y=group, x=percentage)) + 
  geom_bar(stat="identity", width=0.7) +
  theme_bw()+ 
  theme(axis.text = element_text(size=14),
        axis.text.x = element_blank(),  # 
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16) # 
  )+xlab("")+ylab("")+
  scale_fill_manual(values = my2colors)+
  ggtitle('Celltype Percentage')

p
ggsave("SLC7A11-tissue_ration.pdf",width = 6,height = 4)
dev.off()

colnames(allcell_harmony@meta.data)#Tissue_sub,slc7a11
tdata <- allcell_harmony@meta.data[,c('Ccl2','celltype')]#
tdata1 <- table(tdata)
tdata2 <- as.data.frame.matrix(t(apply(tdata1,1,function(x){unlist(lapply(x,function(z){z/sum(x)}))})))
tdata2$group <- rownames(tdata2)
rownames(tdata2)
tdata3 <- gather(tdata2,cluster,percentage,-group)
tdata3$percentage <- tdata3$percentage*100

p <- ggplot(tdata3, aes(fill=cluster, y=percentage, x=group)) + 
  geom_bar(stat="identity", width=0.7) +
  theme_bw()+ 
  theme(axis.text = element_text(size=14),
        axis.text.x = element_text(),  # 
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16) # 
  )+xlab("")+ylab("")+
  scale_fill_brewer(palette = "Set3")+
  ggtitle('SLC7A11+ Cell Percentage')
p
ggsave("SLC7A11-celltype_ration.pdf",width = 4,height = 4)

dev.off()
VlnPlot(allcell_harmony, features = c('SLC7A11'),cols = my36colors,pt.size = 1) 
ggsave("VinplotPlotall_SLC7A11.pdf",width = 6,height = 4)



###################################DC  part####################################

library(Seurat)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(tidyr)
library(viridis)
library(SeuratWrappers)
#library(devtools)
library(harmony)

Idents(allcell) <- allcell@meta.data$celltype
unique(allcell$celltype)
seuratDC<-subset(allcell,idents='DC')
seuratDC[["RNA"]] <- JoinLayers(seuratDC[["RNA"]])
DefaultAssay(seuratDC) <- "RNA"  
seuratDC<- NormalizeData(seuratDC, verbose = FALSE)
seuratDC<- FindVariableFeatures(seuratDC, selection.method = "vst", nfeatures = 2000, 
                                verbose = FALSE)
seuratDC<- ScaleData(seuratDC)
seuratDC<- RunPCA(seuratDC)
ElbowPlot(seuratDC)
pc.num=1:20

######harmony########
system.time({seuratDC <- RunHarmony(seuratDC, group.by.vars = "Sample")})
seuratDC <- RunTSNE(seuratDC, reduction = "harmony", dims = 1:30)
seuratDC <- FindNeighbors(seuratDC, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 0.3)
DimPlot(seuratDC, reduction = "tsne", label=T) 
p=DimPlot(seuratDC, reduction = "tsne", group.by = c("Tissue","seurat_clusters"))
p
ggsave("tsne_DC_evidence.pdf",p,width = 12,height = 5)


Idents(seuratDC) <- seuratDC@meta.data$seurat_cluster
new.cluster.ids <- c(1:length(levels(seuratDC)))
names(x = new.cluster.ids) <- levels(x = seuratDC)
seuratDC <- RenameIdents(object = seuratDC, new.cluster.ids)
seuratDC[["Clusters"]]<-Idents(seuratDC)
Idents(seuratDC) <- seuratDC$Clusters 
VlnPlot(seuratDC, features = c('CD80'),cols = my36colors,pt.size = 1) 
ggsave("DCtype_CD80.pdf",width = 6,height = 4)


Idents(seuratDC) <- seuratDC$seurat_clusters
seuratDC.markers_res0.5 <- FindAllMarkers(seuratDC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seuratDC.markers_res0.5$nLog_padj <- -log10(seuratDC.markers_res0.5$p_val_adj) 
seuratDC.markers_res0.5$ratio <- seuratDC.markers_res0.5$pct.1/seuratDC.markers_res0.5$pct.2
out <- cbind(gene=seuratDC.markers_res0.5$gene,seuratDC.markers_res0.5)
write.table(out,'seuratDCcluster.xls',sep = '\t',quote = F,row.names = F)
write.csv(seuratDC.markers_res0.5,"seuratDC_res0.3.csv")

seuratDC_res0.3=read.csv('seuratDC_res0.3.csv')
top10markers <- seuratDC_res0.3 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10markers,"seuratDC_top10.csv")
library(Seurat)

VlnPlot(seuratDC, features = c( "LILRA4", "TCF4", "IRF7", "CLEC4C",'NRP1','SPIB', #pDCs
                                "CLEC9A", "XCR1", "BATF3", "IRF8", "CADM1", #cDC1
                                "CD1C", 'FCER1A', 'CLEC10A', 'IRF4', 'SIRPA','CLEC4A','BTLA' ,#cDC2
                                "CD14", 'FCGR3A', 'FCER1A', 'CD163','MRC1','CD5', #MoDC
                                'IDO1',"CCR7",'LILRB2' , 'CD200','CD40','LAMP3' ,#MREG
                                'CD1A','CD207'        #LC
),cols = my36colors,pt.size = 0,stack=T)
ggsave("DC_MAKER.pdf",width = 12,height = 10)
###
seuratDC@meta.data$celltype <- as.character(seuratDC@meta.data$Clusters)
seuratDC@meta.data$celltype[seuratDC@meta.data$Clusters %in% as.character(c(1))] <- 'cDC2'#
seuratDC@meta.data$celltype[seuratDC@meta.data$Clusters %in% as.character(c(2))] <- 'cDC2' #
seuratDC@meta.data$celltype[seuratDC@meta.data$Clusters %in% as.character(c(3))] <- 'cDC1'#
seuratDC@meta.data$celltype[seuratDC@meta.data$Clusters %in% as.character(c(4))] <- 'cDC2'#
seuratDC@meta.data$celltype[seuratDC@meta.data$Clusters %in% as.character(c(5))] <- 'mregDC'#
seuratDC@meta.data$celltype[seuratDC@meta.data$Clusters %in% as.character(c(6))] <- 'cDC2'#
seuratDC@meta.data$celltype[seuratDC@meta.data$Clusters %in% as.character(c(7))] <- 'pDC'#
seuratDC@meta.data$celltype[seuratDC@meta.data$Clusters %in% as.character(c(8))] <- 'LCs'#
seuratDC@meta.data$celltype[seuratDC@meta.data$Clusters %in% as.character(c(9))] <- 'pDC'#
seuratDC@meta.data$celltype[seuratDC@meta.data$Clusters %in% as.character(c(10))] <- 'cDC1'#

Idents(seuratDC) <- seuratDC@meta.data$celltype
unique(seuratDC$celltype)
table(seuratDC$celltype)
table(seuratDC$seurat_clusters)

seuratDC <- subset(seuratDC, idents =  c( "cDC1" ,
                                          "cDC2" ,
                                          "mregDC",
                                          'pDC',
                                          'LCs'
))

DimPlot(
  object = seuratDC,
  reduction = "tsne",
  label = TRUE,
  label.size = 6,
  repel = TRUE,
  pt.size = 2,
  cols = my36colors
) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 16) 
  ) +
  ggtitle('Celltype in cDC')+
  NoLegend()
ggsave("tsne_seuratDC.pdf",width =5.5,height = 5)


DimPlot(
  object = seuratDC,
  reduction = "tsne",
  label = TRUE,
  label.size = 6,
  repel = TRUE,
  pt.size = 1,
  cols = my2colors,
  group.by = "Tissue"
) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    legend.position = c(0.05, 0.05), # 
    legend.justification = c("left", "bottom"),
    legend.background = element_rect(fill = "white", colour = "black")
  )+ 
  ggtitle("") 
ggsave("tsne_seuratDCTissue.pdf",width =5.5,height = 5)

DimPlot(
  object = seuratDC,
  reduction = "tsne",
  label = FALSE,
  label.size = 6,
  repel = TRUE,
  pt.size = 1,
  cols = my36colors,
  group.by = "Sample"
) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    legend.background = element_rect(fill = "white", colour = "black"),
  )+ 
  ggtitle("") #
ggsave("tsne_seuratDCSample.pdf",width =7.5,height = 5)

library(RColorBrewer)

FeaturePlot(
  seuratDC,
  features = c('SLC7A11'),
  reduction = "tsne",
  pt.size = 1,
  raster = FALSE,
  order = T
) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black"))+   
  ggtitle("SLC7A11 Expression in DC")

ggsave("FeatureseuratDCPlotSLC7A11.pdf",width = 6,height = 5)

genes_to_check =list( 
  LCs=c('CD1A','CD207', 'S100B'), 
  mregDC=c('CD274','CD200','CCR7'), 
  pDC=c('LILRA4','TCF4','CLEC4C'), 
  cDC1 =c( 'BATF3','XCR1','CLEC9A'),
  cDC2 =c('CD1C','CLEC10A','FCER1A')
)



p_all_markers <- DotPlot(seuratDC, 
                         features = genes_to_check,
                         scale = TRUE, assay = 'RNA', dot.scale = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  scale_color_gradientn(colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))


p_all_markers
ggsave('DC_checkmarkers.pdf',height = 6,width = 8)


####SLC7A11+/- DCs ratio####
seuratDC$Ccl2='Negative'
k_posi<-subset(seuratDC,rna_SLC7A11>0)
seuratDC$Ccl2[match(rownames(k_posi@meta.data),rownames(seuratDC@meta.data))]='Positive'
VlnPlot(seuratDC,features ="rna_SLC7A11",group.by = 'Ccl2',pt.size = 0.001)#Positive
table(seuratDC$Ccl2)
ggsave("DC_pos_SLC7A11.pdf",width = 6,height = 4)
VlnPlot(seuratDC,features =c('SLC7A11'),pt.size = 1)#Positive
ggsave("VinplotPlotDC_SLC7A11.pdf",width = 6,height = 4)


colnames(seuratDC@meta.data)#Tissue_sub,slc7a11
tdata <- seuratDC@meta.data[,c('Ccl2','celltype')]#提取所需数据

tdata1 <- table(tdata)
tdata2 <- as.data.frame.matrix(t(apply(tdata1,1,function(x){unlist(lapply(x,function(z){z/sum(x)}))})))
tdata2$group <- rownames(tdata2)
rownames(tdata2)
tdata3 <- gather(tdata2,cluster,percentage,-group)
tdata3$percentage <- tdata3$percentage*100

p <- ggplot(tdata3, aes(fill=cluster, y=percentage, x=group)) + 
  geom_bar(stat="identity", width=0.7) +
  theme_bw()+ 
  theme(axis.text = element_text(size=14),
        axis.text.x = element_text(),  # 设置 x 轴文本为空
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16) # 
  )+xlab("")+ylab("")+
  scale_fill_brewer(palette = "Set3")+
  ggtitle('SLC7A11+ Cell Percentage in cDC')
p
ggsave("SLC7A11-celltype_rationDC3.pdf",width = 4,height = 4)
dev.off()


tdata <- seuratDC@meta.data[,c('celltype','Tissue')]#
tdata1 <- table(tdata)
tdata2 <- as.data.frame.matrix(t(apply(tdata1,1,function(x){unlist(lapply(x,function(z){z/sum(x)}))})))
tdata2$group <- rownames(tdata2)
rownames(tdata2)
tdata3 <- gather(tdata2,cluster,percentage,-group)
tdata3$percentage <- tdata3$percentage*100

p <- ggplot(tdata3, aes(fill=cluster, y=group, x=percentage)) + 
  geom_bar(stat="identity", width=0.7) +
  theme_bw()+ 
  theme(axis.text = element_text(size=14),
        axis.text.x = element_blank(),  # 
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16) # 
  )+xlab("")+ylab("")+
  scale_fill_manual(values = my2colors)+
  ggtitle('DC Percentage')
p
ggsave("SLC7A11-tissue_rationDC.pdf",width = 6,height = 4)
dev.off()
saveRDS(seuratDC,'seuratDCname.rds')

#######DC type genes  / mregDC GO
#devtools::install_github("junjunlab/scRNAtoolVis")          
library(scRNAtoolVis)
seuratDC[["RNA"]] <- JoinLayers(seuratDC[["RNA"]])
cellselect="mregDC"
af=seuratDC
afc=af[,which(af$celltype %in% c(cellselect))]
Idents(af)=af$celltype
unique(af@meta.data$celltype)
table(Idents(af))
af.markers <- FindAllMarkers(af,only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
head(af.markers)
af.markers_MregDC <- af.markers[af.markers$cluster == 'mregDC', ]
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)
options(connectionObserver = NULL)
library(org.Hs.eg.db)
library(ggsci)#g
library(tidyverse)
library(DOSE)
ids=bitr(af.markers_MregDC$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 
af.markers_MregDC2=ids$ENTREZID

ego_BP <- enrichGO(gene = af.markers_MregDC2,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = TRUE)


ego_BP@result$Description <- str_to_title(ego_BP@result$Description)#

dotplot(ego_BP, showCategory = 10,label_format =50)+scale_color_gradient(low = "yellow", high = "red")

jjVolcano(diffData = af.markers,tile.col = corrplot::COL2('RdBu', 15)[4:10],base_size=20,min.segment.length = 0, box.padding = 0.6,
          fontface = 'italic',polar = T,col.type = "updown",log2FC.cutoff = 0.5,topGeneN = 5)
ggsave('circleDC.pdf',width = 6,height = 6)

Idents(afc)=afc$Ccl2
unique(afc@meta.data$Ccl2)
table(Idents(afc))

afc.markers <- FindAllMarkers(afc,only.pos = F,min.pct = 0.25, logfc.threshold = 0.25)
head(afc.markers)
afc.markers <- afc.markers %>% filter(rownames(afc.markers) != "SLC7A11")
jjVolcano(
  diffData = afc.markers,base_size=20)
ggsave('degmregDC.pdf',width = 6,height = 6)

###############GSEA
library(clusterProfiler)
library(org.Hs.eg.db)
positive_markers <- FindMarkers(afc,only.pos = F, ident.1='Positive',min.pct = 0.25, logfc.threshold = 0.25)

head(positive_markers)
positive_markers <- positive_markers %>% filter(rownames(positive_markers) != "SLC7A11")

positive_gene_list <- positive_markers$avg_log2FC
names(positive_gene_list) <- rownames(positive_markers) 
positive_gene_list <- sort(positive_gene_list, decreasing = TRUE)

positive_gsea_result <- gseGO(geneList = positive_gene_list,
                              OrgDb = org.Hs.eg.db,
                              ont = "BP",  # Biological Process
                              keyType = "SYMBOL",
                              minGSSize = 10,
                              maxGSSize = 500,
                              pvalueCutoff = 1,
                              verbose = T)


sortpp<-positive_gsea_result[order(positive_gsea_result@result$enrichmentScore, decreasing = T),]
head(sortpp)
sortpp <- sortpp[sortpp$pvalue < 0.01, ]


top_results1 <- head(sortpp, 10)
top_results2<- tail(sortpp,10)
top_all =rbind(top_results1,top_results2)
library(stringr)
top_all$Description <- str_to_title(top_all$Description)#
head(top_results1)

ggplot(top_results1, aes(reorder(Description, -log10(pvalue)), -log10(pvalue))) +
  geom_col(aes(fill=NES)) + coord_flip() +
  labs(x="Pathway", y="-log10(pval)") + scale_fill_gradient(low = "#639CA4FF", high = "yellow") + theme_minimal()

# 创建柱状图

ggplot(top_all, aes(x = reorder(Description, enrichmentScore), y = enrichmentScore, fill = Description)) +
  geom_bar(stat = "identity") +
  geom_text(data = top_all[which(top_all$enrichmentScore>0),],aes(x=Description, y=-0.01, label=Description),
            hjust=1, size=4)+
  geom_text(data = top_all[which(top_all$enrichmentScore<0),],aes(x=Description, y=0.01, label=Description),
            hjust=0, size=4)+
  
  coord_flip() +  # 
  xlab("") +
  ylab("enrichmentScore") +
  ggtitle("Top Enriched Gene Sets by enrichmentScore") +
  theme_minimal() +
  NoLegend()+
  scale_fill_manual(values = my36colors) +  
  theme(     
    plot.title = element_text(size = 14),  
    axis.title.x = element_text(size = 14),  
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_blank(),   
    legend.position = "none"    
    
    ggsave("GSEAposSLC7A11mregDC.pdf",width = 10,height = 6)
    library(enrichplot)
    head(positive_gsea_result$ID)
    
    
    positive_gsea_result@result$ID=positive_gsea_result@result$Description
    plot <- gseaplot2(
      positive_gsea_result,
      geneSetID = 'GO:0034605',  # 
      pvalue_table = TRUE,
      title = "SLC7A11+ MregDC",
      color = "skyblue"  # 
    )
    
    plot 
    
    ggsave("GSEAposSLC7A11heat.pdf",width = 7,height = 5)
    
    #######AddModuleScore 
    Idents(seuratDC)<-seuratDC$Tissue
    table(seuratDC$Tissue)
    seuratDC=subset(seuratDC,idents = c('Tumor'))
    table(seuratDC$Tissue)
    Idents(seuratDC)<-seuratDC$celltype
    library(ggsignif)
    library(ggpubr)
    
    table(seuratDC$celltype)
    gene <- read.table("Co-S.txt",sep='\t',header = T) 
    gene <- unique(gene)
    gene.list=as.list(gene)
    seuratDC=AddModuleScore(seuratDC,gene.list,name='Co-stimulation')
    Idents(seuratDC)<-seuratDC$celltype
    head(seuratDC)
    my_comparisons <- list( c("cDC2", "mregDC"), c("cDC1", "mregDC"), c("pDC", "mregDC"), c("LCs", "mregDC"))
    
    library(FSA)  # For Dunn's Test
    VlnPlot(seuratDC, features = "Co-stimulation1", pt.size = 0) +
      ggtitle('Co-stimulation Score') + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      stat_compare_means(method = "kruskal.test") +  # Kruskal-Wallis 检验
      ylim(-0.5, 2)
    ggsave('DCCo-stimulation1.pdf',height =4,width = 6 )
    
    data_for_dunn <- FetchData(seuratDC, vars = c("Co-stimulation1", "celltype"))
    head(data_for_dunn)
    dunn_result <- dunnTest(`Co-stimulation1` ~ celltype, data = data_for_dunn, method = "bonferroni")
    
    print(dunn_result)
    
    Idents(seuratDC)<-seuratDC$celltype
    seuratMregDC = subset(seuratDC,idents = c('mregDC'))
    table(seuratMregDC$Tissue)
    unique(seuratMregDC$Ccl2)
    table(seuratMregDC$Ccl2)
    seuratMregDC$celltype[seuratMregDC@meta.data$Ccl2 == "Positive"] <- "SLC7A11+mregDC"
    seuratMregDC$celltype[seuratMregDC@meta.data$Ccl2 == "Negative"] <- "SLC7A11-mregDC"
    Idents(seuratMregDC)<-seuratMregDC$celltype
    unique(seuratMregDC$celltype)
    table(seuratMregDC$celltype)
    cells_slc7a11_neg <- WhichCells(seuratMregDC, idents = "SLC7A11-mregDC")
    cells_slc7a11_pos <- WhichCells(seuratMregDC, idents = "SLC7A11+mregDC")
    n_pos <- 30
    set.seed(101) 
    downsampled_cells_neg <- sample(cells_slc7a11_neg, size = n_pos, replace = FALSE)
    downsampled_cells_pos <- sample(cells_slc7a11_pos, size = n_pos, replace = FALSE)
    combined_cells <- c(downsampled_cells_neg, downsampled_cells_pos)
    seuratMregDC_downsampled <- subset(seuratMregDC, cells = combined_cells)
    
    table(seuratMregDC_downsampled$celltype)
    my_comparisons <- list( c("SLC7A11-mregDC", "SLC7A11+mregDC"))
    VlnPlot(seuratMregDC_downsampled, features = "CD80",  pt.size = 0) +
      ggtitle('CD 80') + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      stat_compare_means(comparisons = my_comparisons,method = "wilcox.test") +
      ylim(-0.1, 2) 
    
    ggsave('CD80.pdf',height =4,width = 6 )
    
    VlnPlot(seuratMregDC_downsampled, features = "CD86",  pt.size = 0) +
      ggtitle('CD 86') + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      stat_compare_means(comparisons = my_comparisons,method = "wilcox.test") +
      ylim(-0.1, 2.5) 
    ggsave('CD86.pdf',height =4,width = 6 )
    
    gene <- read.table("Co-S.txt",sep='\t',header = T) 
    gene <- unique(gene)
    gene.list=as.list(gene)
    seuratMregDC_downsampled=AddModuleScore(seuratMregDC_downsampled,gene.list,name='CO')
    Idents(seuratMregDC_downsampled)<-seuratMregDC_downsampled$celltype
    
    table(seuratMregDC_downsampled$celltype)
    my_comparisons <- list( c("SLC7A11-mregDC", "SLC7A11+mregDC"))
    VlnPlot(seuratMregDC_downsampled, features = "CO1", pt.size = 0) +
      ggtitle('Co-stimulation Score') + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      stat_compare_means(comparisons = my_comparisons,method = "wilcox.test") +
      ylim(-0.6, 1) 
    ggsave('mregDCCo-stimulation.pdf',height =4,width = 6 )
    
    
    ###################################trajectory analysis Part#####################
    
    library(monocle)
    library(Seurat)
    library(igraph)
    library(CytoTRACE)
    exp1 <- as.matrix(seuratDC@assays$RNA$counts)
    exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,]
    results <- CytoTRACE(exp1,ncores = 1)
    phenot <- seuratDC$celltype
    phenot <- as.character(phenot)
    names(phenot) <- rownames(seuratDC@meta.data)
    emb <- seuratDC@reductions[["tsne"]]@cell.embeddings
    plotCytoTRACE(results,phenotype = phenot, emb = emb)
    plotCytoTRACE(results,phenotype = phenot, emb = emb)
    
    
    library(Seurat)
    #devtools::install_github('cole-trapnell-lab/monocle3')
    library(monocle3)
    library(tidyverse)
    library(patchwork)
    
    pbmc <-seuratDC
    Idents(pbmc) <- pbmc@meta.data$celltype
    
    data <- GetAssayData(pbmc, assay = 'RNA', slot = 'counts')
    cell_metadata <- pbmc@meta.data
    gene_annotation <- data.frame(gene_short_name = rownames(data))
    rownames(gene_annotation) <- rownames(data)
    cds <- new_cell_data_set(data,
                             cell_metadata = cell_metadata,
                             gene_metadata = gene_annotation)
    #NormalizeData+ScaleData+RunPCA
    cds <- preprocess_cds(cds, num_dim = 50)     
    #plot_pc_variance_explained(cds)   
    
    cds <- reduce_dimension(cds,preprocess_method = "PCA") 
    cds <- reduce_dimension(cds, reduction_method="tSNE",preprocess_method = "PCA")
    cds <- cluster_cells(cds) 
    
    cds.embed <- cds@int_colData$reducedDims$UMAP
    int.embed <- Embeddings(pbmc, reduction = "tsne")
    int.embed <- int.embed[rownames(cds.embed),]
    cds@int_colData$reducedDims$UMAP <- int.embed   
    
    plot_cells(cds, color_cells_by="seurat_clusters") 
    
    cds <- learn_graph(cds, use_partition = F,
                       close_loop = FALSE,)    
    
    
    
    head(colData(cds))
    plot_cells(cds,
               color_cells_by = "celltype",
               label_cell_groups=FALSE,
               label_leaves=TRUE,
               label_branch_points=TRUE,
               group_label_size=6,
               cell_size=1.5)                 
    root_cells <- colnames(cds)[colData(cds)$celltype == "cDC1"]
    
    cds = order_cells(cds,root_cells = root_cells )  
    plot_cells(cds,
               color_cells_by = "pseudotime",
               label_cell_groups= F,
               label_leaves=F,
               label_branch_points=F,
               graph_label_size=1.5,
               group_label_size=4,cell_size=1.5)
    
    ggsave("DCtimefromcDC1.pdf",width = 5,height = 4)
    
    plot_cells(cds,
               genes='SLC7A11',
               label_cell_groups=FALSE,cell_size = 2,
               show_trajectory_graph=FALSE)
    ggsave("pDCtimeSLC7A11.pdf",width = 5,height = 4)
    
    pseudotime_data <- data.frame(cell = colnames(cds),
                                  pseudotime = cds@principal_graph_aux[["UMAP"]]$pseudotime)
    
    metadata <- as.data.frame(cds@colData)
    merged_data <- merge(metadata, pseudotime_data, by.x = "row.names", by.y = "cell")
    head(merged_data)
    
    library(ggplot2)
    ggplot(merged_data, aes(x = pseudotime, fill = Ccl2 )) +
      geom_density(alpha = 1) +
      labs(title = "Pseudotime Distribution of DCs from Different Sources",
           x = "Pseudotime",
           y = "Density") +
      theme_minimal()
    ggsave("mregDCSLC7A11wave.pdf",width = 5,height = 4)
    
    
    
    #FeaturePlot图
    plot <- plot_cells(cds, 
                       genes = Track_genes_sig, 
                       show_trajectory_graph = FALSE, 
                       label_cell_groups = FALSE, 
                       label_leaves = FALSE, 
                       cell_size = 1)
    
    # 修改基因名称的字体大小，并添加边框
    plot + theme(axis.text.x = element_text(size = 14),  
                 axis.text.y = element_text(size = 14), 
                 strip.text = element_text(size = 16),  
                 strip.background = element_rect(color = "black", fill = NA, size = 1))  
    
    ggsave("featuretracktop9.pdf",width = 8,height = 7)
    
    
    
    plot <-plot_cells(cds, genes="SLC7A11", show_trajectory_graph=FALSE,
                      label_cell_groups=FALSE, label_leaves=FALSE,cell_size=1)
    plot + theme(axis.text.x = element_text(size = 14),  
                 axis.text.y = element_text(size = 14),  
                 strip.text = element_text(size = 16),   
                 strip.background = element_rect(color = "black", fill = NA, size = 1))  
    
    ggsave("featuretrackSLC7A11.pdf",width = 4,height = 3)
    
    
    
    ###################################cellchat Part###############################
    
    library(CellChat)
    library(tidyverse)
    library(ggalluvial)
    library(Seurat)
    library(data.table)
    library(ggsci)
    
    ##############################
    allcell$Ccl2='Negative'
    k_posi<-subset(allcell,rna_SLC7A11>0)
    allcell$Ccl2[match(rownames(k_posi@meta.data),rownames(allcell@meta.data))]='Positive'
    table(allcell$Ccl2)
    
    Idents(allcell)<-allcell$Tissue
    unique(allcell$Tissue)
    table(allcell$Tissue)
    
    allcell<-subset(allcell,idents='Tumor')
    scRNA_harmony.1<-allcell
    Idents(scRNA_harmony.1)<-scRNA_harmony.1$celltype
    unique(scRNA_harmony.1$celltype)
    scRNA_harmony.1$celltype2 <- seuratDC$celltype
    scRNA_harmony.1$celltype3 <- scRNA_harmony.1$celltype
    table(scRNA_harmony.1$celltype)
    table(scRNA_harmony.1$Ccl2)
    table(scRNA_harmony.1$celltype2)
    table(scRNA_harmony.1$celltype3)
    scRNA_harmony.1$celltype3[scRNA_harmony.1$celltype2 == "mregDC" & scRNA_harmony.1@meta.data$Ccl2 == "Positive"] <- "SLC7A11+mregDC"
    scRNA_harmony.1$celltype3[scRNA_harmony.1$celltype2 == "mregDC" & scRNA_harmony.1@meta.data$Ccl2 == "Negative"] <- "SLC7A11-mregDC"
    #scRNA_harmony.1$celltype3[scRNA_harmony.1$celltype2 == "mregDC"] <- "mregDC"
    scRNA_harmony.1$celltype3[scRNA_harmony.1$celltype2 == "cDC1"] <- "cDC1"
    scRNA_harmony.1$celltype3[scRNA_harmony.1$celltype2 == "cDC2"] <- "cDC2"
    scRNA_harmony.1$celltype3[scRNA_harmony.1$celltype2 == "pDC"] <- "pDC"
    scRNA_harmony.1$celltype3[scRNA_harmony.1$celltype2 == "LCs"] <- "LCs"
    
    Idents(scRNA_harmony.1)<-scRNA_harmony.1$celltype3
    unique(scRNA_harmony.1$celltype3)
    
    #选择亚群
    
    scRNA_harmony.1<-subset(scRNA_harmony.1,idents= c("CD8T",'TAM','Monocyte','CD4T','NK/NKT','B',"SLC7A11+mregDC","SLC7A11-mregDC","Mast Cell")) #
    
    scRNA_harmony.1$celltype <- scRNA_harmony.1$celltype3
    Idents(scRNA_harmony.1)<-scRNA_harmony.1$celltype3
    unique(scRNA_harmony.1$celltype3)
    table(scRNA_harmony.1$celltype3)
    
    
    data.input <- GetAssayData(scRNA_harmony.1, assay = "RNA", slot = "data")
    identity <- subset(scRNA_harmony.1@meta.data, select = "celltype")
    cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "celltype")
    
    CellChatDB <- CellChatDB.human
    showDatabaseCategory(CellChatDB)
    colnames(CellChatDB$interaction)
    CellChatDB$interaction[1:4,1:4]
    head(CellChatDB$cofactor)
    head(CellChatDB$complex)
    head(CellChatDB$geneInfo)
    
    
    unique(CellChatDB$interaction$annotation)
    # use Secreted Signaling for cell-cell communication analysis
    n <- c( "Secreted Signaling" ,"ECM-Receptor" ,"Cell-Cell Contact" )
    
    CellChatDB.use <- subsetDB(CellChatDB, search = n)
    cellchat@DB <- CellChatDB.use # set the used database in the object
    
    
    
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
    cellchat <- filterCommunication(cellchat, min.cells = 3)
    
    df.net <- subsetCommunication(cellchat)
    df.netp <- subsetCommunication(cellchat, slot.name = "netP")
    
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    groupSize <- as.numeric(table(cellchat@idents))
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                     weight.scale = T, label.edge= F, 
                     title.name = "Number of interactions")
    
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                     weight.scale = T, label.edge= F, 
                     title.name = "Interaction weights/strength")
    
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    netAnalysis_signalingRole_network(cellchat, signaling = "CD80", width = 8, height = 2.5, font.size = 10)
    netAnalysis_signalingRole_network(cellchat, signaling = "CD86", width = 8, height = 2.5, font.size = 10)
    ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
    ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
    ht1 + ht2
    
    