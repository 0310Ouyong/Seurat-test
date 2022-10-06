rm(list=ls())
library(dplyr)
library(Seurat)
library(patchwork)

#1.Set working path; Load the PBMC dataset; Create a Seurat object
setwd('D:/test')
pbmc.data <- Read10X(data.dir = 'filtered_gene_bc_matrices/hg19/')
pbmc <- CreateSeuratObject(counts = pbmc.data,project = 'pbmc3k',min.cells = 3,min.features = 200)
pbmc
##An object of class Seurat 
##13714 features across 2700 samples within 1 assay 
##Active assay: RNA (13714 features, 0 variable features)

#2.QC and selecting cells for further analysis (visualize QC metrics)
pbmc[['percent.mt']] <- PercentageFeatureSet(pbmc,pattern = '^MT-')
p1 <- VlnPlot(pbmc, features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
png(
  filename = 'p1.png',
  width = 900,
  height = 500
)
plot(p1)
dev.off()

plot1 <- FeatureScatter(pbmc,feature1 = 'nCount_RNA',feature2 = 'percent.mt')
plot2 <- FeatureScatter(pbmc,feature1 = 'nCount_RNA',feature2 = 'nFeature_RNA')
p2 <- plot1+plot2
png(
  filename = 'p2.png',
  width = 700,
  height = 500
)
plot(p2)
dev.off()
pbmc <- subset(pbmc,subset = nFeature_RNA > 200 & nFeature_RNA <2500 & percent.mt < 5)

#3.Normalizing the data
pbmc <- NormalizeData(pbmc)

#4.Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc,selection.method = 'vst',nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc),10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1,points = top10,repel = TRUE)
p3 <- plot1+plot2
png(
  filename = 'p3.png',
  width = 900,
  height = 500
)
plot(p3)
dev.off()

#5.Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.genes)

#6.Perform linear dimensional reduction
pbmc <- RunPCA(pbmc,features = VariableFeatures(object = pbmc))
p4 <- VizDimLoadings(pbmc,dims = 1:2,reduction = 'pca')
png(
  filename = 'p4.png',
  width = 900,
  height = 500
)
plot(p4)
dev.off()

p5 <- DimPlot(pbmc,reduction = 'pca')
png(
  filename = 'p5.png',
  width = 700,
  height = 500
)
plot(p5)
dev.off()


windows()
DimHeatmap(pbmc,dims = 1,cells = 500,balanced = TRUE)

windows()
DimHeatmap(pbmc,dims = 1:15,cells = 500,balanced = TRUE)

#7.Determine the 'dimensionality' of the dataset
pbmc <- JackStraw(pbmc,num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc,dims = 1:20)
p8 <- JackStrawPlot(pbmc,dims = 1:15)
png(
  filename = 'p8.png',
  width = 700,
  height = 500
)
plot(p8)
dev.off()

p9 <- ElbowPlot(pbmc)
png(
  filename = 'p9.png',
  width = 700,
  height = 500
)
plot(p9)
dev.off()

#8.Cluster the cells
pbmc <- FindNeighbors(pbmc,dims = 1:10)
pbmc <- FindClusters(pbmc,resolution = 0.5)

#9.Run non-linear dimensional reduction (UMAP/tSNE)
pbmc <- RunUMAP(pbmc,dims = 1:10)
p10 <- DimPlot(pbmc,reduction = 'umap')
png(
  filename = 'p10.png',
  width = 700,
  height = 500
)
plot(p10)
dev.off()
saveRDS(pbmc,file = 'pbmc_tutorial.rds')

#10.Finding differentially expressed features(cluster biomarkers)
##find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc,ident.1 = 2,min.pct = 0.25)
##find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc,ident.1 = 5,ident.2 = c(0,3),min.pct = 0.25)
head(cluster5.markers,n = 5)
##find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
pbmc.markers%>%
  group_by(cluster)%>%
  slice_max(n = 2,order_by = avg_log2FC)

p11 <- VlnPlot(pbmc,features = c('MS4A1','CD79A'))
png(
  filename = 'p11.png',
  width = 1000,
  height = 500
)
plot(p11)
dev.off()

p12 <- VlnPlot(pbmc,features = c('NKG7','PF4'),slot = 'counts',log = TRUE)
png(
  filename = 'p12.png',
  width = 1000,
  height = 500
)
plot(p12)
dev.off()

p13 <- FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))
png(
  filename = 'p13.png',
  width = 1000,
  height = 800
)
plot(p13)
dev.off()

pbmc.markers%>%
  group_by(cluster)%>%
  top_n(n = 10,wt = avg_log2FC) -> top10
p14 <- DoHeatmap(pbmc,features = top10$gene) + NoLegend()
png(
  filename = 'p14.png',
  width = 1000,
  height = 500
)
plot(p14)
dev.off()

#11.Assigning cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono","NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc,new.cluster.ids)
p15 <- DimPlot(pbmc,reduction = 'umap',label = TRUE,pt.size = 0.5) + NoLegend()
png(
  filename = 'p15.png',
  width = 700,
  height = 400
)
plot(p15)
dev.off()











