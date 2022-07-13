library(scde)
library(Seurat)
library(tidyverse)
library(svglite)
library(tidyverse)
library(cowplot)
rawcounts<-read.delim("allcountsALBEX.txt")
gn<-read.delim("genenamesALBEX.txt", header = FALSE)
rownames(rawcounts)<-gn$V1
counts<-rawcounts[-ncol(rawcounts)]
counts<-counts[5:nrow(counts),]
# cc<-clean.counts(counts)
# knn <- knn.error.models(cc, k = ncol(cc)/5, n.cores = 4, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 5, save.model.plots = TRUE, verbose = 1)  # turn on verbosity
# varinfo <- pagoda.varnorm(knn, counts = cc, trim = 3/ncol(cc), max.adj.var = 5, n.cores = 4, plot = TRUE)

set.seed(2024)
##SEURAT
counts<-CreateSeuratObject(counts)

counts[["percent.mt"]] <- PercentageFeatureSet(counts, pattern = "^MT")
dir.create("figs")
dir.create("figs/QC")
png("figs/QC/VlnPlot.png", width = 1500)
VlnPlot(counts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
dev.off()
plot1 <- FeatureScatter(counts, feature1 = "nCount_RNA", feature2 = "percent.mt")&
  geom_vline(xintercept=10^7, linetype="dashed",col="black")&
  geom_hline(yintercept=1, linetype="dashed",col="black")&
  theme(axis.text=element_text(size=20),title=element_text(size=30))&
  xlab("Number of reads in cell")&
  ylab("Percentage of reads mapped\n to mitochondrial genes")&
  labs(title="")
plot2 <- FeatureScatter(counts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")&
  geom_vline(xintercept=10^7, linetype="dashed",col="black")&
  geom_hline(yintercept=200, linetype="dashed",col="blue")&
  geom_hline(yintercept=20000, linetype="dashed",col="black")&
  theme(axis.text=element_text(size=20),title=element_text(size=30))&
    xlab("Number of reads in cell")&
    ylab("Number of unique genes\n mapped in cell")&
    labs(title="")
plot3 <- FeatureScatter(counts, feature1 = "percent.mt", feature2 = "nFeature_RNA")&
  geom_hline(yintercept=200, linetype="dashed",col="blue")&
  geom_hline(yintercept=20000, linetype="dashed",col="black")&
  geom_vline(xintercept=1, linetype="dashed",col="black")&
  theme(axis.text=element_text(size=20),title=element_text(size=30))&
    xlab("Percentage of reads mapped\n to mitochondrial genes")&
    ylab("Number of unique genes\n mapped in cell")&
    labs(title="")
png("figs/QC/Featureplots.png", width = 1500)
plot1 + NoLegend() +plot2 + NoLegend() + plot3 + NoLegend()
dev.off()
counts <- subset(counts, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & nCount_RNA< 10^7 & percent.mt < 1)
counts <- NormalizeData(counts)
counts <- FindVariableFeatures(counts, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(counts), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(counts)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
dir.create("figs/dataex/")
png("figs/QC/VFeaturePlot.png", width =1500)
plot2
dev.off()
all.genes <- rownames(counts)
counts <- ScaleData(counts, features = all.genes)
counts <- RunPCA(counts, features = VariableFeatures(object = counts))
png("figs/QC/PCAPlot.png", width = 1500)
d<-DimPlot(counts, reduction = "pca")
d&NoLegend()&
  theme(axis.text=element_text(size=20),title=element_text(size=30))&
  xlab("Principle Component 1")&
  ylab("Principle Component 2")&
  labs(title="")
dev.off()
# png("figs/QC/PCAPlot_CS.png", width = 1500)
# nonvar<-CellSelector(d)
# dev.off()
# saveRDS(nonvar, "nonvar.RDS")
nonvar<-readRDS("nonvar.RDS")
counts<-rawcounts[5:nrow(rawcounts),colnames(rawcounts) %in% nonvar]
counts<-CreateSeuratObject(counts)
counts <- NormalizeData(counts)
counts <- FindVariableFeatures(counts, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(counts), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(counts)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
dir.create("figs/dataex/")
png("figs/QC/VFeaturePlot_no_outliers.png", width =1500)
plot2
dev.off()
all.genes <- rownames(counts)
counts <- ScaleData(counts, features = all.genes)
counts <- RunPCA(counts, features = VariableFeatures(object = counts))
png("figs/QC/PCAPlot_no_outliers.png", width = 1500)
d<-DimPlot(counts, reduction = "pca")
d&NoLegend()&
  theme(axis.text=element_text(size=20),title=element_text(size=30))&
  xlab("Principle Component 1")&
  ylab("Principle Component 2")&
  labs(title="")
dev.off()


png("figs/QC/HEATMAPs.png", width = 1500)
DimHeatmap(counts, dims = 1:10, cells = 500, balanced = TRUE)&
  theme(axis.text=element_text(size=20),title=element_text(size=30))
dev.off()
counts <- JackStraw(counts, num.replicate = 100)
counts <- ScoreJackStraw(counts, dims = 1:20)
png("figs/QC/JACKSTRAW.png", width = 1500)
JackStrawPlot(counts, dims = 1:15)&
  theme(axis.text=element_text(size=20),title=element_text(size=30))
dev.off()
png("figs/QC/ELBOW.png", width = 1500)
ElbowPlot(counts)&
  theme(axis.text=element_text(size=20),title=element_text(size=30))&
  xlab("Number of Principal Components")
dev.off()
counts <- FindNeighbors(counts, dims = 1:10)
counts <- FindClusters(counts, resolution = 0.5)
head(Idents(counts), 5)
counts <- RunUMAP(counts, dims = 1:10)
metadata<-read.csv("metadata.txt")
counts@meta.data$Dev<-metadata$Developmental_stage[match(rownames(counts@meta.data),metadata$Run)]
counts@meta.data$Devnum<-as.numeric(str_extract(counts@meta.data$Dev, "(?<=day )\\w"))
png("figs/dataex/UMAPplot.png", width = 1200)
umap<-DimPlot(counts, reduction = "umap",pt.size=2)&
  theme(axis.text=element_text(size=20),title=element_text(size=30),
  legend.text=element_text(size=20))&
  labs(color="Cluster")
  devnumumap<-FeaturePlot(counts, "Devnum")&
    theme(axis.text=element_text(size=20),title=element_text(size=30),
    legend.text=element_text(size=20))&
    labs(color="Embryonic \nday",title="")
  umap+devnumumap
dev.off()

counts@meta.data$ALBEX1<-NULL
png("figs/dataex/ALBEXVlnPlot.png", width = 2000, height = 750)
  album<-FeaturePlot(counts, "ALBEX1", min.cutoff=0,max.cutoff=0.02,pt.size=3)&
    theme(axis.text=element_text(size=20),title=element_text(size=30),
    legend.text=element_text(size=20))
  album$data<-album$data[order(album$data$ALBEX1),]
  albvln<-VlnPlot(counts, features = c("ALBEX1"),pt.size=2)&
  theme(axis.text=element_text(size=20),axis.text.x=element_text(angle=0),
  title=element_text(size=30),
  legend.text=element_text(size=20))&
  xlab("Cluster")
  alb<-album+albvln
  albvln
dev.off()


png("figs/dataex/ALBEXUMAP.png", width = 2000, height = 750)
  album
dev.off()
counts@meta.data$ALBEX1<-"NO"
counts@meta.data$ALBEX1[as.numeric(counts@assays$RNA["ALBEX1",])>0]<-"YES"
metadata<-counts@meta.data
png("figs/dataex/ALBEXbarplot.png", width = 2000, height = 750)
  bp<-ggplot(metadata, aes(seurat_clusters,fill=ALBEX1))+geom_bar(position="fill")+theme(axis.text=element_text(size=20),title=element_text(size=30),
  legend.text=element_text(size=20))+xlab("Cluster")+ylab("Proportion of cells")
  bp
dev.off()
png("figs/dataex/ALBEX_inf.png", width = 2000, height = 750)
plot_grid(album&labs(title=""),
  albvln&labs(title=""),
  bp, nrow=1)
dev.off()

cluster5.markers <- FindMarkers(counts, ident.1 = 5)
cl.5.m<-subset(cluster5.markers, p_val_adj<0.05)
cl.5.m<-cl.5.m[order(cl.5.m$avg_log2FC, decreasing=TRUE),]

write.table(cl.5.m, "cluster5table.txt", sep="\t")
png("figs/dataex/APOBECumapplot.png", width = 2000, height = 1500)
FeaturePlot(counts, features = c("NODAL","NANOG","SOX2","APOBEC3A","APOBEC3B","APOBEC3H","APOBEC3C","APOBEC3D","APOBEC3F","APOBEC3G","IFITM1","IFITM3"),pt.size=3,ncol=3)&
  theme(axis.text=element_text(size=20),title=element_text(size=30))
dev.off()
counts<-RunTSNE(counts, dims = 1:10)
png("figs/dataex/TSNEplot.png", width = 1500)
DimPlot(counts, reduction = "tsne")
dev.off()

png("figs/dataex/Placentalumapplot.png", width = 2000, height = 1500)
FeaturePlot(counts, features = c("VEGFA","FLT1","PIGF"),pt.size=3)&
  theme(axis.text=element_text(size=20),title=element_text(size=30))
dev.off()

png("figs/dataex/Pluripotentumapplot.png", width = 2000, height = 1500)
FeaturePlot(counts, features = c("NODAL","NANOG","SOX2","GDF3","LEFTY1","LEFTY2"),pt.size=3)&
  theme(axis.text=element_text(size=20),title=element_text(size=30))
dev.off()


png("figs/dataex/APOBECtsneplot.png", width = 2000, height = 750)
FeaturePlot(counts, features = c("APOBEC3A","APOBEC3B","APOBEC3H","APOBEC3C","APOBEC3D","APOBEC3F","APOBEC3G"),reduction="tsne")
dev.off()

png("figs/dataex/ALBEXtsneplot.png", width = 2000, height = 750)
FeaturePlot(counts, features = c("ALBEX1", "NODAL"), reduction="tsne")
dev.off()
png("figs/dataex/ALBEXumapplot.png", width = 2000, height = 750)
DimPlot(counts, reduction = "umap")+FeaturePlot(counts, features = c("ALBEX1", "NODAL"), reduction="umap")
dev.off()

png("figs/dataex/NODALtsneplot.png", width = 2000, height = 750)
FeaturePlot(counts, features = c("ALBEX1","LEFTY1","LEFTY2","TDGF1"), reduction = "tsne")
dev.off()
png("figs/dataex/vartsneplot.png", width = 2000, height = 750)
FeaturePlot(counts, features = c("LNCPRESS1", "MEG3", "GDF3", "MEG8"), reduction = "tsne")
dev.off()
cluster.markers <- FindAllMarkers(counts)
cl.m.sig<-subset(cluster.markers, p_val_adj <0.05 & cluster ==9 )
write(cl.m.sig,"clust9genes.txt")
VlnPlot(counts,"APOBEC3B")

VlnPlot(count_ss, "ALBEX1")
emat<-data.frame(counts@assays$RNA@counts)
alb1<-emat[nrow(emat),]
corsALBEX1<-lapply(seq(nrow(emat)), function(x){
  print(x)
  g<-cor.test(as.numeric(alb1), as.numeric(emat[x,]))
  gene<-rownames(emat)[x]
  pval<-g[3]
  cor<-g[4]
  data.frame(gene,pval,cor)
})

corsALBEX1_df<-bind_rows(corsALBEX1[1:length(corsALBEX1)])
corsALBEX1_df<-corsALBEX1_df[!is.na(corsALBEX1_df$estimate),]
corsALBEX1_df<-corsALBEX1_df[order(corsALBEX1_df$estimate),]
corsALBEX1_df$gene<-fct_relevel(corsALBEX1_df$gene, corsALBEX1_df$gene)
corsALBEX1_df$sig<-corsALBEX1_df$p.value<0.05
p<-ggplot(corsALBEX1_df, aes(gene, estimate, colour = -log10(p.value))) + geom_point()+
  scale_y_continuous(breaks= seq(-1,1, by =0.05))+scale_colour_gradient(low="blue", high = "red")
ggsave("corplotALBEX1.png", p)
corsALBEX1_df$p.adj<-p.adjust(corsALBEX1_df$p.value)
FeaturePlot(counts,c("ALBEX1","APOBEC3B"), blend=TRUE,cols=c("green","red","blue"))
write(noquote(intneg), "global_corgenes_neg.txt")
nonvar<-CellSelector(DimPlot(counts, reduction = "umap"))
counts9<-counts[,colnames(counts) %in% nonvar]
FeaturePlot(counts9,c("ALBEX1","APOBEC3B"), blend=TRUE,cols=c("green","red","blue"),pt.size=3,min.cutoff=0,max.cutoff=0.02)
FeaturePlot(counts9,c("ALBEX1","APOBEC3B"), blend=TRUE,cols=c("green","red","blue"),pt.size=3,min.cutoff=0,max.cutoff=0.02)


# albexquant<-read.delim("isoforms.fpkm_tracking.1", header=FALSE)
# albexquant$albexbin[albexquant$V12>0]<-"yes"
# albexquant$albexbin[is.na(albexquant$albexbin)]<-"no"
#
#
# albexvals<-albexquant[,  c(14,15)]
# av<-data.frame(table(albexvals))
# av<-av[av$albexbin == "yes",]
# av<-na.omit(av[match(colnames(counts), av$V14),])
# ALBEX1<-av$Freq>0
# ALBEXpos<-av$V14[av$Freq>0]
# ALBEXneg<-av$V14[av$Freq==0]
# names(ALBEX1)<-av$V14
# data.frame(ALBEX1)
#
# #FeaturePlot(co, features = c('ALBEX1', "FGF4", "VSNL1", "APOBEC3C"))
#
# a3b<-albexquant[albexquant$V1 == "ALBEX1_A3B", c(12,14)]
# a3b<-na.omit(a3b[match(colnames(counts), a3b$V14),])
#
# a3a<-albexquant[albexquant$V1 == "ALBEX1_A3A", c(12,14)]
# a3a<-na.omit(a3a[match(colnames(counts), a3a$V14),])
# rownames(a3a)<-a3a$V14
# rownames(a3b)<-a3b$V14
# a3a<-log(a3a[1]+0.01)
# a3b<-log(a3b[1]+0.01)
# colnames(a3a)<-c("ALBEX1_A3A")
# colnames(a3b)<-c("ALBEX1_A3B")
albexquant<-read.delim("isoforms.fpkm_tracking.1", header=FALSE)
albexquant$albexbin[albexquant$V9>0]<-"yes"
albexquant$albexbin[is.na(albexquant$albexbin)]<-"no"


albexvals<-albexquant[,  c(14,15)]
av<-data.frame(table(albexvals))
av<-av[av$albexbin == "yes",]
av<-na.omit(av[match(colnames(counts), av$V14),])
ALBEX1<-av$Freq>0
ALBEXpos<-av$V14[av$Freq>0]
ALBEXneg<-av$V14[av$Freq==0]
names(ALBEX1)<-av$V14


a3b<-albexquant[albexquant$V1 == "ALBEX1_A3B", c(9,14)]
a3b<-na.omit(a3b[match(colnames(counts), a3b$V14),])

a3a<-albexquant[albexquant$V1 == "ALBEX1_A3A", c(9,14)]
a3a<-na.omit(a3a[match(colnames(counts), a3a$V14),])
rownames(a3a)<-a3a$V14
rownames(a3b)<-a3b$V14
APOBEC3BAS1_FPKM<-data.frame(APOBEC3BAS1_FPKM=a3a$V9+a3b$V9, row.names = rownames(a3b))
a3a<-log(a3a[1]+0.01)
a3b<-log(a3b[1]+0.01)
APOBEC3BAS1_FPKM<-log(APOBEC3BAS1_FPKM[1]+0.01)
colnames(a3a)<-c("ALBEX1_A3A")
colnames(a3b)<-c("ALBEX1_A3B")
alquant<-cbind(a3a, a3b, ALBEX1,APOBEC3BAS1_FPKM)
alquant$Run<-rownames(alquant)

co<-AddMetaData(counts, data.frame(alquant))
# co<-AddMetaData(co, a3a)
# co<-AddMetaData(co, a3b)
co@meta.data$seurat_clusters<-as.numeric(co@meta.data$seurat_clusters)
plot2<-DimPlot(co)
plot1<-FeaturePlot(co, features = c('APOBEC3BAS1_FPKM'), reduction = "tsne")
plot3<- FeaturePlot(co, features = c("NODAL"),reduction = "tsne")
plot4<-FeaturePlot(co, features = c("NODAL","NANOG", "SOX2", "APOBEC3A", "APOBEC3B", "APOBEC3H", "APOBEC3C", "APOBEC3D", "APOBEC3G"), reduction = "tsne")
png("figs/dataex/ALBEXtsneplot.png", width = 2000, height = 750)
CombinePlots(list(plot1, plot4), nrow =1) + theme(plot.title = element_text(size = 10)))
dev.off()
cluster9.markers <- FindMarkers(counts, ident.1 = 10, min.pct = 0.25)
head(cluster9.markers, n = 20)
FeaturePlot(co, features = c("NANOG", "NODAL","APOBEC3A", "APOBEC3C"))
png
plot1<-FeaturePlot(co, features = c('APOBEC3BAS1_FPKM'), reduction = "tsne")
plot2<-DimPlot(counts, reduction = "tsne")

plot3<- FeaturePlot(co, features = c("NODAL"), reduction = "tsne")
plot4<-FeaturePlot(co, features = c("NANOG", "NODAL","APOBEC3A", "APOBEC3C"), reduction = "tsne")
png("figs/dataex/ALBEXtsneplot.png", width = 2000, height = 750)
plot2 + plot1 + plot4
dev.off()
png("figs/dataex/APOBECumapplot.png", width = 2000, height = 750)
DimPlot(co)+FeaturePlot(co, features = c("APOBEC3A","APOBEC3B","APOBEC3H","APOBEC3C","APOBEC3D","APOBEC3F","APOBEC3G"))
dev.off()
png("figs/dataex/Pluripotenttsneplot.png", width = 2000, height = 750)
FeaturePlot(co, features = c("NANOG","SOX2", "POU5F1", "NANOGP1"), reduction =tsne)
dev.off()
png("figs/dataex/Pluripotentumapplot.png", width = 2000, height = 750)
FeaturePlot(co, features = c("NANOG","SOX2", "POU5F1", "NANOGP1"))
dev.off()
png("figs/dataex/APOBECtsneplot.png", width = 2000, height = 750)
DimPlot(co, reduction = "tsne")+FeaturePlot(co, features = c("APOBEC3A","APOBEC3B","APOBEC3H","APOBEC3C","APOBEC3D","APOBEC3F","APOBEC3G"), reduction = "tsne")
dev.off()
f<-co[rownames(co) =="APOBEC3B",]
write(colnames(f)[data.frame(f@assays$RNA@counts)>1], "A3Bexpressors")
co@meta.data$A3B<-FALSE
co@meta.data$A3B[rownames(co@meta.data) %in% colnames(f)[data.frame(f@assays$RNA@counts)> 1]]<-TRUE
co@meta.data$Run[co@meta.data$A3B & !co@meta.data$ALBEX1]
A3Bonly<-co@meta.data$Run[co@meta.data$A3B & !co@meta.data$ALBEX1]
ALBEXA3Bcoexpress<-co@meta.data$Run[co@meta.data$A3B & co@meta.data$ALBEX1]
write(paste0(ALBEXA3Bcoexpress,".bam"), "ALBEX_A3B_coexpress")
write(paste0(A3Bonly,".bam"), "A3B_only")

intronquant<-read.delim("introns.fpkm_tracking.1", header=FALSE)

co@meta.data$A3Bintron<-na.omit(intronquant[match(colnames(co), intronquant$V14),])$V9
metas<-co@meta.data
metasA3B<-subset(metas, A3B)



# clust9genes<-rownames(cluster9.markers)[cluster9.markers$p_val_adj <0.05]
# write(paste(noquote(clust9genes), sep = "\n", collapse = "\n"), "clust9genes.txt")
# clust9FC<-cluster9.markers[abs(cluster9.markers$avg_log2FC) >1 & cluster9.markers$p_val_adj < 0.05,]
#
# all.markers <- FindAllMarkers(counts, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# allmarktop2<-all.markers %>%
#     group_by(cluster) %>%
#     top_n(n = 2, wt = avg_log2FC) %>% data.frame()


png("figs/dataex/clustex1-4.png", width = 1500)
FeaturePlot(co, features = allmarktop2$gene[1:4], reduction = "tsne")
dev.off()
png("figs/dataex/clustex5-8.png", width = 1500)
FeaturePlot(co, features = allmarktop2$gene[5:8], reduction = "tsne")
dev.off()
png("figs/dataex/clustex9-12.png", width = 1500)
FeaturePlot(co, features = allmarktop2$gene[9:12], reduction = "tsne")
dev.off()
png("figs/dataex/clustex13-16.png", width = 1500)
FeaturePlot(co, features = allmarktop2$gene[13:16], reduction = "tsne")
dev.off()
png("figs/dataex/clustex17-20.png", width = 1500)
FeaturePlot(co, features = allmarktop2$gene[17:20], reduction = "tsne")
dev.off()






# albexmeans<-unlist(lapply(seq(a3a$ALBEX1_A3A), function(x){
#    mean(a3a$ALBEX1_A3A[x], a3b$ALBEX1_A3B[x])}))
emat<-data.frame(co@assays$RNA@counts)
metadata<-subset(co@meta.data, !is.na(APOBEC3BAS1_FPKM))
emat<-emat[,match(rownames(metadata), colnames(emat))]
corsA3A<-lapply(seq(nrow(emat)), function(x){
  g<-cor.test(metadata$ALBEX1_A3A, as.numeric(emat[x,]))
  gene<-rownames(emat)[x]
  pval<-g[3]
  cor<-g[4]
  data.frame(gene,pval,cor)
})
corsA3B<-lapply(seq(nrow(emat)), function(x){
  g<-cor.test(metadata$ALBEX1_A3B, as.numeric(emat[x,]))
  gene<-rownames(emat)[x]
  pval<-g[3]
  cor<-g[4]
  data.frame(gene,pval,cor)
})
corsALBEX1<-lapply(seq(nrow(emat)), function(x){
  g<-cor.test(metadata$APOBEC3BAS1_FPKM, as.numeric(emat[x,]))
  gene<-rownames(emat)[x]
  pval<-g[3]
  cor<-g[4]
  data.frame(gene,pval,cor)
})
corsA3A<-bind_rows(corsA3A[1:54507])
corsA3A<-corsA3A[order(corsA3A$estimate),]
corsA3A$gene<-fct_relevel(corsA3A$gene, corsA3A$gene)
corsA3A$sig<-corsA3A$p.value<0.05
corsA3B<-bind_rows(corsA3B[1:54507])
corsA3B<-corsA3B[order(corsA3B$estimate),]
corsA3B$gene<-fct_relevel(corsA3B$gene, corsA3B$gene)
corsA3B$sig<-corsA3B$p.value<0.05
corsALBEX1<-bind_rows(corsALBEX1[1:54507])
corsALBEX1<-corsALBEX1[order(corsALBEX1$estimate),]
corsALBEX1$gene<-fct_relevel(corsALBEX1$gene, corsALBEX1$gene)
corsALBEX1$sig<-corsALBEX1$p.value<0.05
p<-ggplot(corsA3A, aes(gene, estimate, colour = -log10(p.value))) + geom_point()+
  scale_y_continuous(breaks= seq(-1,1, by =0.05))+scale_colour_gradient(low="blue", high = "red")
ggsave("corplotA3A.png", p)
p<-ggplot(corsA3B, aes(gene, estimate, colour = -log10(p.value))) + geom_point()+
  scale_y_continuous(breaks= seq(-1,1, by =0.05))+scale_colour_gradient(low="blue", high = "red")
ggsave("corplotA3B.png", p)
p<-ggplot(corsALBEX1, aes(gene, estimate, colour = -log10(p.value))) + geom_point()+
  scale_y_continuous(breaks= seq(-1,1, by =0.05))+scale_colour_gradient(low="blue", high = "red")
ggsave("corplotALBEX1.png", p)
corsA3Asigpos<-subset(corsA3A, estimate > 0 & p.value<0.05)
corsA3Bsigpos<-subset(corsA3B, estimate > 0 & p.value<0.05)
intpos<-intersect(corsA3Asigpos$gene, corsA3Bsigpos$gene)
write(noquote(intpos), "global_corgenes_pos.txt")
corsA3Asigneg<-subset(corsA3A, estimate < 0 & p.value<0.05)
corsA3Bsigneg<-subset(corsA3B, estimate < 0 & p.value<0.05)
intneg<-intersect(corsA3Asigneg$gene, corsA3Bsigneg$gene)
write(noquote(intneg), "global_corgenes_neg.txt")

corsALBEX1sigpos<-subset(corsALBEX1, estimate > 0 & p.value<0.05)
corsALBEX1signeg<-subset(corsALBEX1, estimate < 0 & p.value<0.05)
ALBEX1neg_cor<-as.character(corsALBEX1signeg$gene)
ALBEX1pos_cor<-as.character(corsALBEX1sigpos$gene)
write(noquote(ALBEX1neg_cor), "global_corgenes_ALBEX1_neg.txt")
write(noquote(ALBEX1pos_cor), "global_corgenes_ALBEX1_pos.txt")


# corsA3A<-unlist(lapply(seq(nrow(emat)), function(x){
#   cor(metadata$ALBEX1_A3A, as.numeric(emat[x,]))
# }))
# corsA3B<-unlist(lapply(seq(nrow(emat)), function(x){
#   cor(metadata$ALBEX1_A3B, as.numeric(emat[x,]))
# }))
# names(corsA3A)<-rownames(emat)
# names(corsA3B)<-rownames(emat)
# corsA3AnoNA<-na.omit(corsA3A)
# corsA3BnoNA<-na.omit(corsA3B)
# corsA3AnoNA<-corsA3AnoNA[order(corsA3AnoNA)]
# corsA3BnoNA<-corsA3BnoNA[order(corsA3BnoNA)]
# corsA3AnoNA<-data.frame(cor= corsA3AnoNA, gene = names(corsA3AnoNA))
# corsA3BnoNA<-data.frame(cor =corsA3BnoNA, gene = names(corsA3BnoNA))
# corsA3AnoNA$gene<-fct_relevel(corsA3AnoNA$gene, corsA3AnoNA$gene)
# corsA3BnoNA$gene<-fct_relevel(corsA3BnoNA$gene, corsA3BnoNA$gene)
# saveRDS(corsA3A, "corsA3A.rds")
# saveRDS(corsA3B, "corsA3B.rds")
#
# png("corplotA3A", width = 1500)
# ggplot(corsA3AnoNA, aes(gene, cor)) + geom_point()+
#   scale_y_continuous(breaks= seq(0,1, by =0.05))
# dev.off()
#
# png("corplotA3B", width = 1500)
# ggplot(corsA3BnoNA, aes(gene, cor)) + geom_point()+
#   scale_y_continuous(breaks= seq(0,1, by =0.05))
# dev.off()

namesA3A<-as.character(na.omit(names(corsA3A)[corsA3A>0.25]))
namesA3B<-as.character(na.omit(names(corsA3B)[corsA3B>0.25]))
int<-intersect(namesA3B, namesA3A)
write(noquote(int), "allcorgenes.txt",sep = "\n")
FeaturePlot(co, features = c("CD44", "ALDH"), reduction = "tsne")
clust10<-rownames(co@meta.data)[co@meta.data$seurat_clusters == 10]
clust10mat<-rawcounts[5:nrow(rawcounts),colnames(rawcounts) %in% clust10]
saveRDS(clust10mat, "cluster10mat.RDS")
