##note if having trouble with hdf5r must install hdf5 first (from source without sudo access)
#  then  run export LD_PRELOAD=~/hdf5-1.12.0/hdf5/lib/libhdf5.so:~/hdf5-1.12.0/hdf5/lib/libhdf5_hl.so
install.packages("hdf5r", configure.args="--with-hdf5=/Volumes/userdata/student_users/conormcguinness/hdf5-1.12.0/hdf5/bin/h5cc")
library(hdf5r)
library(scde)
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(tidyverse)
lf<-list.files(pattern = "ALBEX1_SRR134434*")
#lfs<-paste0(lf, "/outs/filtered_feature_bc_matrix.h5")
# lapply(lf,function(x){
#   lfs<-paste0(x, "/outs/filtered_feature_bc_matrix.h5")
#   h<-Read10X_h5(lfs)
#   counts<-CreateSeuratObject(h)
#   counts[["percent.mt"]] <- PercentageFeatureSet(counts, pattern = "^MT")
#   figdir<-paste0(x,"/figs")
#   figQC<-paste0(figdir, "QC")
#   figex<-paste0(figdir, "dataex")
#   dir.create(figdir)
#   dir.create(figQC)
#   dir.create(figex)
#   png(paste0(figQC,"/VlnPlot.png"), width = 1500)
#   VlnPlot(counts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
#   dev.off()
#   plot1 <- FeatureScatter(counts, feature1 = "nCount_RNA", feature2 = "percent.mt")
#   plot2 <- FeatureScatter(counts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#   plot3 <- FeatureScatter(counts, feature1 = "percent.mt", feature2 = "nFeature_RNA")
#   png(paste0(figQC,"/Featureplots.png"), width = 1500)
#   plot1 + NoLegend() +plot2 + NoLegend() + plot3 + NoLegend()
#   dev.off()
#   counts <- subset(counts, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mt < 5)
#   counts <- NormalizeData(counts)
#   counts <- FindVariableFeatures(counts, selection.method = "vst", nfeatures = 2000)
#
#   # Identify the 10 most highly variable genes
#   top10 <- head(VariableFeatures(counts), 10)
#
#   # plot variable features with and without labels
#   plot1 <- VariableFeaturePlot(counts)
#   plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#   dir.create("figs/dataex/")
#   png(figex,"/VFeaturePlot.png", width =1500)
#   plot2
#   dev.off()
#   all.genes <- rownames(counts)
#   counts <- ScaleData(counts, features = all.genes)
# })
library(tidyverse)
 h<-Read10X_h5("Aggr_Samples/outs/count/filtered_feature_bc_matrix.h5")
 counts<-CreateSeuratObject(h)
# rawcounts<-read.delim("allcounts.txt")
# gn<-read.delim("genenames.txt", header = FALSE)
# rownames(rawcounts)<-gn$V1
# #counts<-counts[-ncol(counts)]
# counts<-rawcounts[5:nrow(rawcounts),1:1502]
# cc<-clean.counts(counts)
# knn <- knn.error.models(cc, k = ncol(cc)/5, n.cores = 4, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 5, save.model.plots = TRUE, verbose = 1)  # turn on verbosity
# varinfo <- pagoda.varnorm(knn, counts = cc, trim = 3/ncol(cc), max.adj.var = 5, n.cores = 4, plot = TRUE)

set.seed(420)
##SEURAT
counts<-CreateSeuratObject(counts)

counts[["percent.mt"]] <- PercentageFeatureSet(counts, pattern = "^MT")
dir.create("figs")
dir.create("figs/QC")
png("figs/QC/VlnPlot", width = 1500)
VlnPlot(counts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
dev.off()
plot1 <- FeatureScatter(counts, feature1 = "nCount_RNA", feature2 = "percent.mt")&
  geom_vline(xintercept=50000, linetype="dashed",col="black")&
  geom_hline(yintercept=5, linetype="dashed",col="black")&
  theme(axis.text=element_text(size=20),title=element_text(size=30))&
  xlab("Number of reads in cell")&
  ylab("Percentage of reads mapped\n to mitochondrial genes")&
  labs(title="")
plot2 <- FeatureScatter(counts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")&
  geom_vline(xintercept=50000, linetype="dashed",col="black")&
  geom_hline(yintercept=200, linetype="dashed",col="blue")&
  geom_hline(yintercept=20000, linetype="dashed",col="black")&
  theme(axis.text=element_text(size=20),title=element_text(size=30))&
    xlab("Number of reads in cell")&
    ylab("Number of unique genes\n mapped in cell")&
    labs(title="")
plot3 <- FeatureScatter(counts, feature1 = "percent.mt", feature2 = "nFeature_RNA")&
  geom_hline(yintercept=200, linetype="dashed",col="blue")&
  geom_hline(yintercept=20000, linetype="dashed",col="black")&
  geom_vline(xintercept=5, linetype="dashed",col="black")&
  theme(axis.text=element_text(size=20),title=element_text(size=30))&
    xlab("Percentage of reads mapped\n to mitochondrial genes")&
    ylab("Number of unique genes\n mapped in cell")&
    labs(title="")
png("figs/QC/Featureplots.png", width = 1500)
plot1 + NoLegend() +plot2 + NoLegend() + plot3 + NoLegend()
dev.off()
counts <- subset(counts, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mt < 5 &nCount_RNA<50000)
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
  xlab("Principal Component 1")&
  ylab("Principal Component 2")&
  labs(title="")
dev.off()

png("figs/dataex/HEATMAPs.png", width = 1500)
DimHeatmap(counts, dims = 1:10, cells = 500, balanced = TRUE)
dev.off()
counts <- JackStraw(counts, num.replicate = 100)
counts <- ScoreJackStraw(counts, dims = 1:20)
png("figs/dataex/JACKSTRAW.png", width = 1500)
JackStrawPlot(counts, dims = 1:15)
dev.off()
png("figs/dataex/ELBOW.png", width = 1500)
ElbowPlot(counts)
dev.off()
counts <- FindNeighbors(counts, dims = 1:15)
counts <- FindClusters(counts, resolution = 0.5)
head(Idents(counts), 5)
counts <- RunUMAP(counts, dims = 1:15)
png("figs/dataex/UMAPplot.png", width = 1500)
DimPlot(counts, reduction = "umap", label = TRUE,pt.size=5)&
  theme(axis.text=element_text(size=20),title=element_text(size=30))
dev.off()
png("figs/dataex/UMAPplot+CD49.png", width = 1500)
DimPlot(counts, reduction = "umap", label = TRUE,pt.size=5)&
  theme(axis.text=element_text(size=20),title=element_text(size=30))&
  FeaturePlot(counts, features = c("PTPRC"),pt.size=5,reduction = "umap")
dev.off()
counts<-RunTSNE(counts, dims = 1:10)
png("figs/dataex/TSNEplot.png", width = 1500)
DimPlot(counts, reduction = "tsne")
dev.off()
saveRDS(counts, "counts.RDS")
png("figs/dataex/APOBECumapplot.png", width = 2000, height = 750)
DimPlot(counts)+FeaturePlot(counts, features = c("APOBEC3A","APOBEC3B","APOBEC3H","APOBEC3C","APOBEC3D","APOBEC3F","APOBEC3G"))
dev.off()

png("figs/dataex/epithelial_plot.png", width = 2000, height = 1000)
epiplot<-FeaturePlot(counts, features = c("PTPRC","ITGA6",
"NFIB","TP63",
"ELF5","EHF",
"FOXA1","ESR1"),pt.size=3,ncol=2,reduction = "umap")
  epiplot&
    theme(axis.text=element_text(size=20),title=element_text(size=30))
dev.off()
png("figs/dataex/ALBEXVlnPlot.png", width = 2000, height = 750)
  albvln<-VlnPlot(counts, features = c("ALBEX1","APOBEC3B"),pt.size=2)
  albvln&
  theme(axis.text=element_text(size=20),axis.text.x=element_text(angle=0),
  title=element_text(size=30),
  legend.text=element_text(size=20))&
  xlab("Cluster")
dev.off()

png("figs/dataex/ALBEXumapplot.png", width = 2000, height = 750)
FeaturePlot(counts, features = c("ALBEX1", "APOBEC3B", "APOBEC3A", "PIK3CA"), reduction = "umap", min.cutoff = 0)
dev.off()
d<-DimPlot(counts)

non_immune<-colnames(counts)[!counts@meta.data$seurat_clusters %in% c(9,13,14,21,0,17,5,24)]
epicounts<-CreateSeuratObject(h[,colnames(h) %in% non_immune])

#nonvar<-CellSelector(d)
epicounts[["percent.mt"]] <- PercentageFeatureSet(epicounts, pattern = "^MT")

dir.create("epifigs")
dir.create("epifigs/QC")
png("epifigs/QC/VlnPlot", width = 1500)
VlnPlot(epicounts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
dev.off()
plot1 <- FeatureScatter(epicounts, feature1 = "nCount_RNA", feature2 = "percent.mt")&
  geom_vline(xintercept=50000, linetype="dashed",col="black")&
  geom_hline(yintercept=5, linetype="dashed",col="black")&
  theme(axis.text=element_text(size=20),title=element_text(size=30))&
  xlab("Number of reads in cell")&
  ylab("Percentage of reads mapped\n to mitochondrial genes")&
  labs(title="")
plot2 <- FeatureScatter(epicounts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")&
  geom_vline(xintercept=50000, linetype="dashed",col="black")&
  geom_hline(yintercept=200, linetype="dashed",col="blue")&
  geom_hline(yintercept=20000, linetype="dashed",col="black")&
  theme(axis.text=element_text(size=20),title=element_text(size=30))&
    xlab("Number of reads in cell")&
    ylab("Number of unique genes\n mapped in cell")&
    labs(title="")
plot3 <- FeatureScatter(epicounts, feature1 = "percent.mt", feature2 = "nFeature_RNA")&
  geom_hline(yintercept=200, linetype="dashed",col="blue")&
  geom_hline(yintercept=20000, linetype="dashed",col="black")&
  geom_vline(xintercept=5, linetype="dashed",col="black")&
  theme(axis.text=element_text(size=20),title=element_text(size=30))&
    xlab("Percentage of reads mapped\n to mitochondrial genes")&
    ylab("Number of unique genes\n mapped in cell")&
    labs(title="")
png("epifigs/QC/Featureplots.png", width = 1500)
plot1 + NoLegend() +plot2 + NoLegend() + plot3 + NoLegend()
dev.off()
epicounts <- subset(epicounts, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mt < 5 &nCount_RNA<50000)
epicounts <- NormalizeData(epicounts)
epicounts <- FindVariableFeatures(epicounts, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(epicounts), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(epicounts)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
dir.create("epifigs/dataex/")
png("epifigs/QC/VFeaturePlot.png", width =1500)
plot2
dev.off()
all.genes <- rownames(epicounts)
epicounts <- ScaleData(epicounts, features = all.genes)
epicounts <- RunPCA(epicounts, features = VariableFeatures(object = epicounts))
png("epifigs/QC/PCAPlot.png", width = 1500)
d<-DimPlot(epicounts, reduction = "pca")
d&NoLegend()&
  theme(axis.text=element_text(size=20),title=element_text(size=30))&
  xlab("Principal Component 1")&
  ylab("Principal Component 2")&
  labs(title="")
dev.off()

png("epifigs/dataex/HEATMAPs.png", width = 1500)
DimHeatmap(epicounts, dims = 1:10, cells = 500, balanced = TRUE)
dev.off()
epicounts <- JackStraw(epicounts, num.replicate = 100)
epicounts <- ScoreJackStraw(epicounts, dims = 1:20)
png("epifigs/dataex/JACKSTRAW.png", width = 1500)
JackStrawPlot(epicounts, dims = 1:15)
dev.off()
png("epifigs/dataex/ELBOW.png", width = 1500)
ElbowPlot(epicounts)
dev.off()
epicounts <- FindNeighbors(epicounts, dims = 1:20)
epicounts <- FindClusters(epicounts, resolution = 0.5)
head(Idents(epicounts), 5)
epicounts <- RunUMAP(epicounts, dims = 1:20)
png("epifigs/dataex/UMAPplot.png", width = 1500)
DimPlot(epicounts, reduction = "umap", label = TRUE,pt.size=3)&
  theme(axis.text=element_text(size=20),title=element_text(size=30))
dev.off()
# epicounts<-RunTSNE(epicounts, dims = 1:10)
# png("epifigs/dataex/TSNEplot.png", width = 1500)
# DimPlot(epicounts, reduction = "tsne")
# dev.off()
# saveRDS(epicounts, "epicounts.RDS")
# epicounts<-readRDS("epicounts.RDS")
png("epifigs/dataex/APOBECumapplot.png", width = 2000, height = 750)
DimPlot(epicounts)+FeaturePlot(epicounts, features = c("APOBEC3A","APOBEC3B","APOBEC3H","APOBEC3C","APOBEC3D","APOBEC3F","APOBEC3G"))
dev.off()

png("epifigs/dataex/epithelial_plot.png", width = 1500, height = 1000)
epiplot<-FeaturePlot(epicounts, features = c("EPCAM","SFRP1",
 "EHF", "ELF5",
 "FOXA1","ESR1",
 "MYC", "VIM"),pt.size=3,ncol=2,reduction = "umap")
  epiplot&
    theme(axis.text=element_text(size=20),title=element_text(size=30))
dev.off()
png("epifigs/dataex/ALBEXumapplot.png", width = 2000, height = 750)
FeaturePlot(epicounts, features = c("ALBEX1", "APOBEC3B"), pt.size=2,reduction = "umap", min.cutoff = 0,max.cutoff=0.5)
dev.off()
FeaturePlot(epicounts, features = c("TNFSF10"), pt.size=2,reduction = "umap", min.cutoff = 0,max.cutoff=0.5)

FeaturePlot(epicounts, features = c("HLA-A"),pt.size=1,ncol=2,reduction = "umap")
png("epifigs/dataex/ALBEXVlnPlot.png", width = 2000, height = 750)
  album<-FeaturePlot(epicounts, "ALBEX1", min.cutoff=0,max.cutoff=0.02,pt.size=3)&
    theme(axis.text=element_text(size=20),title=element_text(size=30),
    legend.text=element_text(size=20))
  album$data<-album$data[order(album$data$ALBEX1),]
  albvln<-VlnPlot(epicounts, features = c("ALBEX1","APOBEC3B"),pt.size=2)&
  theme(axis.text=element_text(size=20),axis.text.x=element_text(angle=0),
  title=element_text(size=30),
  legend.text=element_text(size=20))&
  xlab("Cluster")
  alb<-album+albvln
  albvln
dev.off()

png("epifigs/dataex/ALBEXUMAP.png", width = 2000, height = 750)
  album
dev.off()
# epicounts@meta.data$ALBEX1<-"NO"
# epicounts@meta.data$ALBEX1[as.numeric(epicounts@assays$RNA["ALBEX1",])>0]<-"YES"
metadata<-epicounts@meta.data
png("epifigs/dataex/ALBEXbarplot.png", width = 2000, height = 750)
  bp<-ggplot(metadata, aes(seurat_clusters,fill=ALBEX1))+geom_bar(position="fill")+theme(axis.text=element_text(size=20),title=element_text(size=30),
  legend.text=element_text(size=20))+xlab("Cluster")+ylab("Proportion of cells")
  bp
dev.off()
png("epifigs/dataex/ALBEX_inf.png", width = 2000, height = 750)
plot_grid(album&labs(title=""),
  albvln&labs(title=""),
  bp, nrow=1)
dev.off()

cluster2.markers <- FindMarkers(epicounts, ident.1 = c(2))
png("epifigs/dataex/VolcanoPlot_global.png",width = 500, height = 350)
ggplot(cluster2.markers, aes(avg_log2FC, -log(p_val_adj)))+
  geom_point()+
  geom_hline(yintercept=-log(0.05), linetype="dashed",col="blue")+
  geom_vline(xintercept=1, linetype="dashed",col="black")+
  theme(axis.text=element_text(size=20),title=element_text(size=30),
  legend.text=element_text(size=20))+xlab("Average log2 Fold Change \ncluster 2 vs rest of dataset")+
  ylab("-log adjusted \np-value")
dev.off()
cl.2.m_global<-subset(cluster2.markers, p_val_adj<0.05 & avg_log2FC >1)
cl.2.m_global<-cl.2.m_global[order(cl.2.m_global$avg_log2FC, decreasing=TRUE),]
write.table(cl.2.m_global, "clust2_global.txt", sep ="\t")
cl.2.m_global_noRP<-cl.2.m_global[grep("^RP",rownames(cl.2.m_global),invert=TRUE),]
write.table(cl.2.m_global_noRP, "clust2_global_noRP.txt", sep ="\t")

cluster2.markers_lum <- FindMarkers(epicounts, ident.1 = c(2), ident.2=c(5,18,14,6,3,15))
png("epifigs/dataex/VolcanoPlot_luminal.png",width = 500, height = 350)
ggplot(cluster2.markers_lum, aes(avg_log2FC, -log(p_val_adj)))+
  geom_point()+
  geom_hline(yintercept=-log(0.05), linetype="dashed",col="blue")+
  geom_vline(xintercept=1, linetype="dashed",col="black")+
  theme(axis.text=element_text(size=20),title=element_text(size=30),
  legend.text=element_text(size=20))+xlab("Average log2 Fold Change \ncluster 2 vs luminal clusters")+
  ylab("-log adjusted \np-value")
dev.off()
cl.2.m_luminal<-subset(cluster2.markers_lum, p_val_adj<0.05 & avg_log2FC>1)
cl.2.m_luminal<-cl.2.m_luminal[order(cl.2.m_luminal$avg_log2FC, decreasing=TRUE),]
write.table(cl.2.m_luminal, "clust2_luminal.txt", sep ="\t")
cl.2.m_luminal_noRP<-cl.2.m_luminal[grep("^RP",rownames(cl.2.m_luminal),invert=TRUE),]
write.table(cl.2.m_luminal_noRP, "clust2_luminal_noRP.txt", sep ="\t")
png("epifigs/dataex/TNFplot.png",width = 1000, height = 500)
  FeaturePlot(epicounts, c("TNFRSF11A", "TNFSF11"),max.cutoff=0.5)&
    theme(axis.text=element_text(size=20),title=element_text(size=30),
    legend.text=element_text(size=20))
dev.off()
Idents(epicounts)==0
lumprogcounts<-subset(epicounts, seurat_clusters==0)
dir.create("lumprogfigs")
png("lumprogfigs/ALBEXVlnPlot.png", width = 2000, height = 750)
  album<-FeaturePlot(lumprogcounts, "ALBEX1", min.cutoff=0,max.cutoff=0.02,pt.size=3)&
    theme(axis.text=element_text(size=20),title=element_text(size=30),
    legend.text=element_text(size=20))
  album$data<-album$data[order(album$data$ALBEX1),]
  albvln<-VlnPlot(lumprogcounts, features = c("ALBEX1","APOBEC3B"),pt.size=2)&
  theme(axis.text=element_text(size=20),axis.text.x=element_text(angle=0),
  title=element_text(size=30),
  legend.text=element_text(size=20))&
  xlab("Cluster")
  alb<-album+albvln
  albvln
dev.off()

png("lumprogfigs/ALBEXUMAP.png", width = 2000, height = 750)
  album
dev.off()
png("lumprogfigs/ALBEXumapplot.png", width = 2000, height = 750)
FeaturePlot(lumprogcounts, features = c("ALBEX1", "APOBEC3B"), pt.size=2,reduction = "umap", min.cutoff = 0,max.cutoff=0.5)
dev.off()
subclust<-CellSelector(plot = album)

sc<-lumprogcounts[,colnames(lumprogcounts) %in% subclust]


album<-FeaturePlot(sc, c("ALBEX1"), min.cutoff=0,max.cutoff=0.02,pt.size=3)&
  theme(axis.text=element_text(size=20),title=element_text(size=30),
  legend.text=element_text(size=20))
album$data<-album$data[order(album$data$APOBEC3B),]
ALBEXers<-CellSelector(plot = album, object=sc,ident="ALBEXERS")
ALBEXmarkers<-FindMarkers(ALBEXers, ident.1="ALBEXERS")
cl.0.m_pos<-subset(ALBEXmarkers, p_val_adj<0.05 & avg_log2FC >0)
cl.0.m_pos<-cl.0.m_pos[order(cl.0.m_pos$avg_log2FC, decreasing=TRUE),]
cl.0.m_neg<-subset(ALBEXmarkers, p_val_adj<0.05 & avg_log2FC <0)
cl.0.m_neg<-cl.0.m_neg[order(cl.0.m_neg$avg_log2FC, decreasing=TRUE),]
write.table(cl.0.m_pos, "ALBEXERSPOS.txt", sep ="\t")
write.table(cl.0.m_neg, "ALBEXERSNEG.txt", sep ="\t")
# write.table(cl.3.m, "cluster3table.txt", sep="\t")
#
# lumprog<-colnames(epicounts)[epicounts@meta.data$seurat_clusters == 3]
# lumprogcounts<-CreateSeuratObject(h[,colnames(h) %in% lumprog])
# lumprogcounts[["percent.mt"]] <- PercentageFeatureSet(lumprogcounts, pattern = "^MT")
#
# dir.create("lumprogfigs")
# dir.create("lumprogfigs/QC")
# png("lumprogfigs/QC/VlnPlot", width = 1500)
# VlnPlot(lumprogcounts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
# dev.off()
# plot1 <- FeatureScatter(lumprogcounts, feature1 = "nCount_RNA", feature2 = "percent.mt")&
#   geom_vline(xintercept=50000, linetype="dashed",col="black")&
#   geom_hline(yintercept=5, linetype="dashed",col="black")&
#   theme(axis.text=element_text(size=20),title=element_text(size=30))&
#   xlab("Number of reads in cell")&
#   ylab("Percentage of reads mapped\n to mitochondrial genes")&
#   labs(title="")
# plot2 <- FeatureScatter(lumprogcounts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")&
#   geom_vline(xintercept=50000, linetype="dashed",col="black")&
#   geom_hline(yintercept=200, linetype="dashed",col="blue")&
#   geom_hline(yintercept=20000, linetype="dashed",col="black")&
#   theme(axis.text=element_text(size=20),title=element_text(size=30))&
#     xlab("Number of reads in cell")&
#     ylab("Number of unique genes\n mapped in cell")&
#     labs(title="")
# plot3 <- FeatureScatter(lumprogcounts, feature1 = "percent.mt", feature2 = "nFeature_RNA")&
#   geom_hline(yintercept=200, linetype="dashed",col="blue")&
#   geom_hline(yintercept=20000, linetype="dashed",col="black")&
#   geom_vline(xintercept=5, linetype="dashed",col="black")&
#   theme(axis.text=element_text(size=20),title=element_text(size=30))&
#     xlab("Percentage of reads mapped\n to mitochondrial genes")&
#     ylab("Number of unique genes\n mapped in cell")&
#     labs(title="")
# png("lumprogfigs/QC/Featureplots.png", width = 1500)
# plot1 + NoLegend() +plot2 + NoLegend() + plot3 + NoLegend()
# dev.off()
# lumprogcounts <- subset(lumprogcounts, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mt < 5 &nCount_RNA<50000)
# lumprogcounts <- NormalizeData(lumprogcounts)
# lumprogcounts <- FindVariableFeatures(lumprogcounts, selection.method = "vst", nfeatures = 2000)
#
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(lumprogcounts), 10)
#
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(lumprogcounts)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# dir.create("lumprogfigs/dataex/")
# png("lumprogfigs/QC/VFeaturePlot.png", width =1500)
# plot2
# dev.off()
# all.genes <- rownames(lumprogcounts)
# lumprogcounts <- ScaleData(lumprogcounts, features = all.genes)
# lumprogcounts <- RunPCA(lumprogcounts, features = VariableFeatures(object = lumprogcounts))
# png("lumprogfigs/QC/PCAPlot.png", width = 1500)
# d<-DimPlot(lumprogcounts, reduction = "pca")
# d&NoLegend()&
#   theme(axis.text=element_text(size=20),title=element_text(size=30))&
#   xlab("Principal Component 1")&
#   ylab("Principal Component 2")&
#   labs(title="")
# dev.off()
#
# png("lumprogfigs/QC/HEATMAPs.png", width = 1500)
# DimHeatmap(lumprogcounts, dims = 1:10, cells = 500, balanced = TRUE)
# dev.off()
# lumprogcounts <- JackStraw(lumprogcounts, num.replicate = 100)
# lumprogcounts <- ScoreJackStraw(lumprogcounts, dims = 1:20)
# png("lumprogfigs/QC/JACKSTRAW.png", width = 1500)
# JackStrawPlot(lumprogcounts, dims = 1:15)
# dev.off()
# png("lumprogfigs/QC/ELBOW.png", width = 1500)
# ElbowPlot(lumprogcounts)
# dev.off()
# lumprogcounts <- FindNeighbors(lumprogcounts, dims = 1:15)
# lumprogcounts <- FindClusters(lumprogcounts, resolution = 0.5)
# head(Idents(lumprogcounts), 5)
# lumprogcounts <- RunUMAP(lumprogcounts, dims = 1:15)
# png("lumprogfigs/dataex/UMAPplot.png", width = 1500)
# DimPlot(lumprogcounts, reduction = "umap", label = TRUE,pt.size=3)&
#   theme(axis.text=element_text(size=20),title=element_text(size=30))
# dev.off()
# lumprogcounts<-RunTSNE(lumprogcounts, dims = 1:10)
# png("lumprogfigs/dataex/TSNEplot.png", width = 1500)
# DimPlot(lumprogcounts, reduction = "tsne")
# dev.off()
# saveRDS(lumprogcounts, "lumprogcounts.RDS")
# lumprogcounts<-readRDS("lumprogcounts.RDS")
# png("lumprogfigs/dataex/APOBECumapplot.png", width = 2000, height = 750)
# DimPlot(lumprogcounts)+FeaturePlot(lumprogcounts, features = c("APOBEC3A","APOBEC3B","APOBEC3H","APOBEC3C","APOBEC3D","APOBEC3F","APOBEC3G"))
# dev.off()
#
# png("lumprogfigs/dataex/epithelial_plot.png", width = 1500, height = 1000)
# epiplot<-FeaturePlot(lumprogcounts, features = c("EPCAM","ITGA6",
# "NFIB","SFRP1",
# "ELF5","EHF",
# "FOXA1","ESR1"),pt.size=3,ncol=2,reduction = "umap")
#   epiplot&
#     theme(axis.text=element_text(size=20),title=element_text(size=30))
# dev.off()
# png("lumprogfigs/dataex/ALBEXumapplot.png", width = 2000, height = 750)
# FeaturePlot(lumprogcounts, features = c("ALBEX1", "APOBEC3B", "APOBEC3A", "PIK3CA"), reduction = "umap", min.cutoff = 0,max.cutoff=0.5)
# dev.off()
# FeaturePlot(lumprogcounts, features = c("TRIB3"),pt.size=1,ncol=2,reduction = "umap")
# png("lumprogfigs/dataex/ALBEXVlnPlot.png", width = 2000, height = 750)
#   album<-FeaturePlot(lumprogcounts, "ALBEX1", min.cutoff=0,max.cutoff=0.02,pt.size=3)&
#     theme(axis.text=element_text(size=20),title=element_text(size=30),
#     legend.text=element_text(size=20))
#   album$data<-album$data[order(album$data$ALBEX1),]
#   albvln<-VlnPlot(lumprogcounts, features = c("ALBEX1","APOBEC3B"),pt.size=2)&
#   theme(axis.text=element_text(size=20),axis.text.x=element_text(angle=0),
#   title=element_text(size=30),
#   legend.text=element_text(size=20))&
#   xlab("Cluster")
#   alb<-album+albvln
#   albvln
# dev.off()
#
# png("lumprogfigs/dataex/ALBEXUMAP.png", width = 2000, height = 750)
#   album
# dev.off()
# lumprogcounts@meta.data$ALBEX1<-"NO"
# lumprogcounts@meta.data$ALBEX1[as.numeric(lumprogcounts@assays$RNA["ALBEX1",])>0]<-"YES"
# metadata<-lumprogcounts@meta.data
# png("lumprogfigs/dataex/ALBEXbarplot.png", width = 2000, height = 750)
#   bp<-ggplot(metadata, aes(seurat_clusters,fill=ALBEX1))+geom_bar(position="fill")+theme(axis.text=element_text(size=20),title=element_text(size=30),
#   legend.text=element_text(size=20))+xlab("Cluster")+ylab("Proportion of cells")
#   bp
# dev.off()
# png("lumprogfigs/dataex/ALBEX_inf.png", width = 2000, height = 750)
# plot_grid(album&labs(title=""),
#   albvln&labs(title=""),
#   bp, nrow=1)
# dev.off()
#
# cluster3.markers <- FindMarkers(lumprogcounts, ident.1 = c(3))
# cl.3.m<-subset(cluster3.markers, p_val_adj<0.05 & avg_log2FC >0)
# cl.3.m<-cl.3.m[order(cl.3.m$avg_log2FC, decreasing=TRUE),]
#
# write.table(cl.3.m, "cluster3tablelumprogs.txt", sep="\t")
#

# ##SEURAT
# epicounts<-CreateSeuratObject(epicounts)
#
# epicounts[["percent.mt"]] <- PercentageFeatureSet(epicounts, pattern = "^MT")
# dir.create("epifigs")
# dir.create("epifigs/QC")
# png("epifigs/QC/VlnPlot", width = 1500)
# VlnPlot(epicounts, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
# dev.off()
# plot1 <- FeatureScatter(epicounts, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(epicounts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot3 <- FeatureScatter(epicounts, feature1 = "percent.mt", feature2 = "nFeature_RNA")
# png("epifigs/QC/Featureplots.png", width = 1500)
# plot1 + NoLegend() +plot2 + NoLegend() + plot3 + NoLegend()
# dev.off()
# epicounts <- subset(epicounts, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mt < 5)
# epicounts <- NormalizeData(epicounts)
# epicounts <- FindVariableFeatures(epicounts, selection.method = "vst", nfeatures = 2000)
#
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(epicounts), 10)
#
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(epicounts)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# dir.create("epifigs/dataex/")
# png("epifigs/dataex/VFeaturePlot.png", width =1500)
# plot2
# dev.off()
# all.genes <- rownames(epicounts)
# epicounts <- ScaleData(epicounts, features = all.genes)
# epicounts <- RunPCA(epicounts, features = VariableFeatures(object = epicounts))
# png("epifigs/dataex/PCAPlot.png", width = 1500)
# DimPlot(epicounts, reduction = "pca")
# dev.off()
# png("epifigs/dataex/HEATMAPs.png", width = 1500)
# DimHeatmap(epicounts, dims = 1:10, cells = 500, balanced = TRUE)
# dev.off()
# epicounts <- JackStraw(epicounts, num.replicate = 100)
# epicounts <- ScoreJackStraw(epicounts, dims = 1:20)
# png("epifigs/dataex/JACKSTRAW.png", width = 1500)
# JackStrawPlot(epicounts, dims = 1:15)
# dev.off()
# png("epifigs/dataex/ELBOW.png", width = 1500)
# ElbowPlot(epicounts)
# dev.off()
# epicounts <- FindNeighbors(epicounts, dims = 1:20)
# epicounts <- FindClusters(epicounts, resolution = 0.5)
# head(Idents(epicounts), 5)
# epicounts <- RunUMAP(epicounts, dims = 1:20)
# png("epifigs/dataex/UMAPplot.png", width = 1500)
# DimPlot(epicounts, reduction = "umap", label = TRUE)
# dev.off()
# epicounts<-RunTSNE(epicounts, dims = 1:10)
# png("epifigs/dataex/TSNEplot.png", width = 1500)
# DimPlot(epicounts, reduction = "tsne")
# dev.off()
# png("epifigs/dataex/APOBECumapplot.png", width = 2000, height = 750)
# DimPlot(epicounts)+FeaturePlot(epicounts, features = c("APOBEC3A","APOBEC3B","APOBEC3H","APOBEC3C","APOBEC3D","APOBEC3F","APOBEC3G"),cols=c("white", "black"))
# dev.off()
#
# png("epifigs/dataex/APOBECtsneplot.png", width = 2000, height = 750)
# DimPlot(epicounts, reduction ="tsne")+FeaturePlot(epicounts, features = c("APOBEC3A","APOBEC3B","APOBEC3H","APOBEC3C","APOBEC3D","APOBEC3F","APOBEC3G"), reduction = "tsne")
# dev.off()
# png("epifigs/dataex/ALBEXumapplot.png", width = 2000, height = 750)
# DimPlot(epicounts)+ FeaturePlot(epicounts, features = c("ALBEX1", "APOBEC3B", "APOBEC3A", "PIK3CA"), reduction = "umap", cols=c("white","black"))
# dev.off()
# png("epifigs/dataex/ALBEXtsneplot.png", width = 2000, height = 750)
# FeaturePlot(epicounts, features = c("ALBEX1"), reduction = "tsne")
# dev.off()
# png("epifigs/dataex/PPtsneplot.png", width = 2000, height = 750)
# DimPlot(epicounts, reduction = "tsne")+FeaturePlot(epicounts, features = c("NANOG","NANOGP1","SOX2","POU5F1"), reduction = "tsne")
# dev.off()
# png("epifigs/dataex/NODALumapplot.png", width = 2000, height = 750)
# FeaturePlot(epicounts, features = c("NODAL","LEFTY1","LEFTY2","TDGF1"), reduction = "umap")
# dev.off()
# png("epifigs/dataex/varumapplot.png", width = 2000, height = 750)
# DimPlot(epicounts)+FeaturePlot(epicounts, features = c("LNCPRESS1", "MEG3", "GDF3", "MEG8", "PTPRC"), reduction = "umap")
# dev.off()
# png("epifigs/dataex/epithelial_plot.png", width = 2000, height = 750)
# DimPlot(epicounts)+FeaturePlot(epicounts, features = c("PTPRC","ITGA6", "EPCAM", "ELF5", "EHF", "FOXA1", "ESR1", "TP63", "NFIB"), reduction = "umap")
# dev.off()
# png("epifigs/dataex/lump_plot.png", width = 2000, height = 750)
# FeaturePlot(epicounts, features = c("GATA3","SLPI","PTGDS", "IGF1","SFRP1","KRT15", "KRT17", "KLK5", "ANXA1","STAT3","LIF"), reduction = "umap")
# dev.off()
# png("epifigs/dataex/lump_plot_2.png", width = 2000, height = 750)
# FeaturePlot(epicounts, features = c("ALBEX1","APOBEC3B", "SLC2A1","STC1","PTN", "SERINC2", "TAGLN", "TFF3", "TFF1"), reduction = "umap")
# dev.off()
# png("epifigs/dataex/MMP_TGF_plot.png", width = 2000, height = 750)
# FeaturePlot(epicounts, features = c("MMP2","MMP9","MMP13","MMP2", "TGFB1"), reduction = "umap")
# dev.off()
# png("epifigs/dataex/ALB_.png", width = 2000, height = 750)
# FeatureScatter(epicounts, feature1="ALBEX1", feature2="APOBEC3A")+
#   FeatureScatter(epicounts, feature1="ALBEX1", feature2="APOBEC3B") +
#   FeatureScatter(epicounts, feature1="APOBEC3B", feature2="APOBEC3A")
# dev.off()
# saveRDS(counts, "counts.RDS")
# saveRDS(epicounts, "epicounts.RDS")
# ###monocle
# epicounts.cds<-as.cell_data_set(epicounts)
# epicounts.cds<-cluster_cells(cds = epicounts.cds, reduction.method = "UMAP")
# epicounts.cds<-learn_graph(epicounts.cds,use_partition=TRUE)
# epicounts.cds <- order_cells(epicounts.cds, reduction_method = "UMAP")
# epilums<-colnames(epicounts)[epicounts@meta.data$seurat_clusters %in% c(1,3,4,8,14,17,20)]
# epilums<-CreateSeuratObject(h[,colnames(h) %in% epilums])
#
#
#
# epilums[["percent.mt"]] <- PercentageFeatureSet(epilums, pattern = "^MT")
# dir.create("epilumfigs")
# dir.create("epilumfigs/QC")
# png("epilumfigs/QC/VlnPlot", width = 1500)
# VlnPlot(epilums, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
# dev.off()
# plot1 <- FeatureScatter(epilums, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(epilums, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot3 <- FeatureScatter(epilums, feature1 = "percent.mt", feature2 = "nFeature_RNA")
# png("epilumfigs/QC/Featureplots.png", width = 1500)
# plot1 + NoLegend() +plot2 + NoLegend() + plot3 + NoLegend()
# dev.off()
# epilums <- subset(epilums, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mt < 5)
# epilums <- NormalizeData(epilums)
# epilums <- FindVariableFeatures(epilums, selection.method = "vst", nfeatures = 2000)
#
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(epilums), 10)
#
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(epilums)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# dir.create("epilumfigs/dataex/")
# png("epilumfigs/dataex/VFeaturePlot.png", width =1500)
# plot2
# dev.off()
# all.genes <- rownames(epilums)
# epilums <- ScaleData(epilums, features = all.genes)
# epilums <- RunPCA(epilums, features = VariableFeatures(object = epilums))
# png("epilumfigs/dataex/PCAPlot.png", width = 1500)
# DimPlot(epilums, reduction = "pca")
# dev.off()
# png("epilumfigs/dataex/HEATMAPs.png", width = 1500)
# DimHeatmap(epilums, dims = 1:10, cells = 500, balanced = TRUE)
# dev.off()
# epilums <- JackStraw(epilums, num.replicate = 100)
# epilums <- ScoreJackStraw(epilums, dims = 1:20)
# png("epilumfigs/dataex/JACKSTRAW.png", width = 1500)
# JackStrawPlot(epilums, dims = 1:15)
# dev.off()
# png("epilumfigs/dataex/ELBOW.png", width = 1500)
# ElbowPlot(epilums)
# dev.off()
# epilums <- FindNeighbors(epilums, dims = 1:20)
# epilums <- FindClusters(epilums, resolution = 0.5)
# head(Idents(epilums), 5)
# epilums <- RunUMAP(epilums, dims = 1:20)
# png("epilumfigs/dataex/UMAPplot.png", width = 1500)
# DimPlot(epilums, reduction = "umap", label = TRUE)
# dev.off()
# epilums<-RunTSNE(epilums, dims = 1:10)
# png("epilumfigs/dataex/TSNEplot.png", width = 1500)
# DimPlot(epilums, reduction = "tsne")
# dev.off()
# png("epilumfigs/dataex/APOBECumapplot.png", width = 2000, height = 750)
# DimPlot(epilums)+FeaturePlot(epilums, features = c("APOBEC3A","APOBEC3B","APOBEC3H","APOBEC3C","APOBEC3D","APOBEC3F","APOBEC3G"),cols=c("white", "black"))
# dev.off()
#
# png("epilumfigs/dataex/APOBECtsneplot.png", width = 2000, height = 750)
# DimPlot(epilums, reduction ="tsne")+FeaturePlot(epilums, features = c("APOBEC3A","APOBEC3B","APOBEC3H","APOBEC3C","APOBEC3D","APOBEC3F","APOBEC3G"), reduction = "tsne")
# dev.off()
# png("epilumfigs/dataex/ALBEXumapplot.png", width = 2000, height = 750)
# DimPlot(epilums)+ FeaturePlot(epilums, features = c("ALBEX1"), reduction = "umap")
# dev.off()
# png("epilumfigs/dataex/ALBEXtsneplot.png", width = 2000, height = 750)
# FeaturePlot(epilums, features = c("ALBEX1"), reduction = "tsne")
# dev.off()
# png("epilumfigs/dataex/PPtsneplot.png", width = 2000, height = 750)
# DimPlot(epilums, reduction = "tsne")+FeaturePlot(epilums, features = c("NANOG","NANOGP1","SOX2","POU5F1"), reduction = "tsne")
# dev.off()
# png("epilumfigs/dataex/NODALumapplot.png", width = 2000, height = 750)
# FeaturePlot(epilums, features = c("NODAL","LEFTY1","LEFTY2","TDGF1"), reduction = "umap")
# dev.off()
# png("epilumfigs/dataex/varumapplot.png", width = 2000, height = 750)
# DimPlot(epilums)+FeaturePlot(epilums, features = c("LNCPRESS1", "MEG3", "GDF3", "MEG8", "PTPRC"), reduction = "umap")
# dev.off()
# png("epilumfigs/dataex/epithelial_plot.png", width = 2000, height = 750)
# FeaturePlot(epilums, features = c("PIK3CA","PTPRC","ITGA6", "EPCAM", "ELF5", "EHF", "FOXA1", "ESR1", "TP53", "NFIB", "KRAS"), reduction = "umap")
# dev.off()
# png("epilumfigs/dataex/lump_plot.png", width = 2000, height = 750)
# FeaturePlot(epilums, features = c("GATA3","SLPI","PTGDS", "IGF1","SFRP1","KRT15", "KRT17", "KLK5", "ANXA1","STAT3","LIF"), reduction = "umap")
# dev.off()
# png("epilumfigs/dataex/lump_plot_2.png", width = 2000, height = 750)
# FeaturePlot(epilums, features = c("ALBEX1","APOBEC3B", "SLC2A1","STC1","PTN", "SERINC2", "TAGLN", "TFF3", "TFF1"), reduction = "umap")
# dev.off()
# png("epilumfigs/dataex/MMP_TGF_plot.png", width = 2000, height = 750)
# FeaturePlot(epilums, features = c("MMP2","LIF","KLF5","WNT5A","WNT5B", "TGFB1", "KRT8","ALDH1A3", "JAG2", "HES4"), reduction = "umap")
# dev.off()
# png("epilumfigs/dataex/lump_plot_2.png", width = 2000, height = 750)
# FeaturePlot(epilums, features = c("ALBEX1","APOBEC3B", "SLC2A1","STC1","PTN", "SERINC2", "TAGLN", "TFF3", "TFF1", "MKI67"), reduction = "umap")
# dev.off()
# png("epilumfigs/dataex/lump_plot_3.png", width = 2000, height = 750)
# FeaturePlot(epilums, features = c("SAMHD1","IFNGR1","IFNAR1", "IFITM1", "IFITM2", "TNFSF10", "SNAI1", "ZEB1", "KMT2C"), reduction = "umap")
# dev.off()
# png("epilumfigs/dataex/ALB_.png", width = 2000, height = 750)
# FeatureScatter(epilums, feature1="ALBEX1", feature2="APOBEC3A")+
#   FeatureScatter(epilums, feature1="ALBEX1", feature2="APOBEC3B") +
#   FeatureScatter(epilums, feature1="APOBEC3B", feature2="APOBEC3A")
# dev.off()
# saveRDS(epilums,"epilums.RDS")
# epicountlums.cds<-as.cell_data_set(epilums)
# epilums<-readRDS("epilums.RDS")
# epicountlums.cds<-readRDS("epicountlums_cds.RDS")
# rowData(epicountlums.cds)$gene_name<-rownames(epicountlums.cds)
# rowData(epicountlums.cds)$gene_short_name <- rowData(epicountlums.cds)$gene_name
# epicountlums.cds<-cluster_cells(cds = epicountlums.cds, reduction.method = "UMAP")
# epicountlums.cds<-learn_graph(epicountlums.cds, use_partition = TRUE)
# epicountlums.cds <- order_cells(epicountlums.cds, reduction_method = "UMAP")
# plot_cells(epicountlums.cds, reduction_method = "UMAP", color_cells_by="pseudotime")
# cds_subset <- choose_cells(epicountlums.cds)
# subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
# pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
# gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=1e-11)
# agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
# module_dendro <- hclust(dist(agg_mat))
# gene_module_df$module <- factor(gene_module_df$module,
#                                 levels = row.names(agg_mat)[module_dendro$order])
#
# plot_cells(cds_subset,
#            genes=data.frame(gene_module_df),
#            label_cell_groups=FALSE,
#            show_trajectory_graph=FALSE)
# plot_cells(cds_subset,
#            genes=c("APOBEC3B","ALBEX1"),
#            label_cell_groups=FALSE,
#            show_trajectory_graph=FALSE)
# write(gene_module_df[gene_module_df$module==1,]$id, "monocle_genes.txt")
# pr_graph_test_res <- graph_test(epicountlums,cds, neighbor_graph="knn", cores=8)
# gene_module_df <- find_gene_modules(epicountlums.cds, resolution=1e-2)
# cell_group_df <- tibble::tibble(cell=row.names(colData(epicountlums.cds)),
#                                 cell_group=partitions(epicountlums.cds))
# agg_mat <- aggregate_gene_expression(epicountlums.cds, gene_module_df, cell_group_df)
# row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
# colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
#
# pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
#                    scale="column", clustering_method="ward.D2",
#                    fontsize=6)
# dplot<-DimPlot(epilums, reduction = "umap", label = TRUE)
# cluster6.markers <- FindMarkers(epilums, ident.1 = "13", min.pct = 0.25)
# ss<-subset(cluster6.markers, avg_log2FC>0 & p_val_adj < 0.05)
# write(rownames(ss), "ss_pos.txt")
# VlnPlot(epilums, features = c("APOBEC3B","ALBEX1"))
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.gene
# epilums<- CellCycleScoring(epilums, s.features = s.genes, g2m.features = g2m.genes, set.idents = FALSE)
#
# VlnPlot(epilums, features="TAGLN")
# taglnhigh<-colnames(subset(epilums, TAGLN > 1))
# epilums@meta.data$tagln<-"LOW"
# epilums@meta.data$tagln[match(taglnhigh, colnames(epilums))]<-"HIGH"
# epilums@meta.data$albclust<-"REST"
# epilums@meta.data$albclust[epilums@meta.data$tagln =="HIGH" & epilums@meta.data$seurat_clusters=="3"]<-"HIGH"
# epilums@meta.data$albclust[epilums@meta.data$tagln =="LOW" & epilums@meta.data$seurat_clusters=="3"]<-"LOW"
# SetIdent(epilums, id = "albclust")
#
# cluster6.markers<-FindMarkers(epilums, ident.1 = "HIGH", ident.2 ="LOW", min.pct =0.1)
# ss<-subset(cluster6.markers, avg_log2FC>0 & p_val_adj < 0.05)
# write(rownames(ss), "ss_pos.txt")
