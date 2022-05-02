library(Seurat)
library(dplyr)
library(Matrix)
library(gdata)
library(cowplot)

setwd("~/Data/RegevLab/NHLeukemia/")

WTTX.hap.BM.mRNA.files <- "BM_hg19_AR792_allSingletExprs_minGeneCutoff_300.Robj"
load(file = WTTX.hap.BM.mRNA.files,verbose = T)
WTTX.hap.BM.mRNA <<- xxx

WTTX.hap.SP.mRNA.files <- "SP_hg19_AR792_allSingletExprs_minGeneCutoff_300.Robj"
load(file = WTTX.hap.SP.mRNA.files,verbose = T)
WTTX.hap.SP.mRNA <<- xxx

WTVeh.dip.BM.mRNA.files <- "BM_hg19_AR1270_allSingletExprs_minGeneCutoff_300.Robj"
load(file = WTVeh.dip.BM.mRNA.files,verbose = T)
WTVeh.dip.BM.mRNA <<- xxx

WTVeh.dip.SP.mRNA.files <- "SP_hg19_AR1270_allSingletExprs_minGeneCutoff_300.Robj"
load(file = WTVeh.dip.SP.mRNA.files,verbose = T)
WTVeh.dip.SP.mRNA <<- xxx

WTTX.dip.BM.mRNA.files <- "BM_hg19_AR786_allSingletExprs_minGeneCutoff_300.Robj"
load(file = WTTX.dip.BM.mRNA.files,verbose = T)
WTTX.dip.BM.mRNA <<- xxx

WTTX.dip.SP.mRNA.files <- "SP_hg19_AR786_allSingletExprs_minGeneCutoff_300.Robj"
load(file = WTTX.dip.SP.mRNA.files,verbose = T)
WTTX.dip.SP.mRNA <<- xxx

WTVeh.hap.BM.mRNA.files <- "BM_hg19_AR1272_allSingletExprs_minGeneCutoff_200.Robj"
load(file = WTVeh.hap.BM.mRNA.files,verbose = T)
WTVeh.hap.BM.mRNA <<- xxx

WTVeh.hap.SP.mRNA.files <- "SP_hg19_AR1272_allSingletExprs_minGeneCutoff_200.Robj"
load(file = WTVeh.hap.SP.mRNA.files,verbose = T)
WTVeh.hap.SP.mRNA <<- xxx

WT.hap.BM.mRNA.files <- "BM_hg19_AR759_allSingletExprs_minGeneCutoff_200.Robj"
load(file = WT.hap.BM.mRNA.files,verbose = T)
WT.hap.BM.mRNA <<- xxx

WT.hap.SP.mRNA.files <- "SP_hg19_AR759_allSingletExprs_minGeneCutoff_200.Robj"
load(file = WT.hap.SP.mRNA.files,verbose = T)
WT.hap.SP.mRNA <<- xxx

WT.dip.BM.mRNA.files <- "BM_hg19_AR761_allSingletExprs_minGeneCutoff_300.Robj"
load(file = WT.dip.BM.mRNA.files,verbose = T)
WT.dip.BM.mRNA <<- xxx

WT.dip.SP.mRNA.files <- "SP_hg19_AR761_allSingletExprs_minGeneCutoff_300.Robj"
load(file = WT.dip.SP.mRNA.files,verbose = T)
WT.dip.SP.mRNA <<- xxx

PDX.hg19.1906.merged <- merge(WT.hap.BM.mRNA,c(WT.hap.SP.mRNA,
                                              WTVeh.hap.BM.mRNA,WTVeh.hap.SP.mRNA),
                              add.cell.ids = c("WT.hap.BM","WT.hap.SP",
                                               "WTVeh.hap.BM","WTVeh.hap.SP"))
PDX.hg19.1905.merged <- merge(WT.dip.BM.mRNA,c(WT.dip.SP.mRNA,WTVeh.dip.BM.mRNA,WTVeh.dip.SP.mRNA),
                              add.cell.ids = c("WT.dip.BM.","WT.dip.SP.",
                                               "WTVeh.dip.BM.","WTVeh.dip.SP."))
PDX.hg19.1811.merged <- merge(WTTX.hap.BM.mRNA,WTTX.hap.SP.mRNA,
                              add.cell.ids = c("WTTX.hap.BM.","WTTX.hap.SP."))
# PDX.hg19.1810.merged <- merge(WT.BM.mRNA,WT.SP.mRNA,
#                               add.cell.ids = c("WT.BM.","WT.SP."))
PDX.hg19.1809.merged <- merge(WTTX.dip.BM.mRNA,WTTX.dip.SP.mRNA,
                              add.cell.ids = c("WTTX.dip.BM.","WTTX.dip.SP."))

PDX.hg19.1809.merged <- NormalizeData(PDX.hg19.1809.merged)
#PDX.hg19.1809.merged <- ScaleData(PDX.hg19.1809.merged)
PDX.hg19.1809.merged <- FindVariableFeatures(PDX.hg19.1809.merged,selection.method = "vst",nfeatures = 200)
#hvg.hg19.1809 <- rownames(x = head(x = PDX.hg19.1809.merged@hvg.info, n = 200))
PDX.hg19.1809.merged@meta.data[,"batch"] <- "2018_09"

# PDX.hg19.1810.merged <- NormalizeData(PDX.hg19.1810.merged)
# #PDX.hg19.1810.merged <- ScaleData(PDX.hg19.1810.merged)
# PDX.hg19.1810.merged <- FindVariableFeatures(PDX.hg19.1810.merged,selection.method = "vst",nfeatures = 200)
# #hvg.hg19.1810 <- rownames(x = head(x = PDX.hg19.1810.merged@hvg.info, n = 200))
# PDX.hg19.1810.merged@meta.data[,"batch"] <- "2018_10"
# 
PDX.hg19.1811.merged <- NormalizeData(PDX.hg19.1811.merged)
#PDX.hg19.1811.merged <- ScaleData(PDX.hg19.1811.merged)
PDX.hg19.1811.merged <- FindVariableFeatures(PDX.hg19.1811.merged,selection.method = "vst",nfeatures = 200)
#hvg.hg19.1811 <- rownames(x = head(x = PDX.hg19.1811.merged@hvg.info, n = 200))
PDX.hg19.1811.merged@meta.data[,"batch"] <- "2018_11"

PDX.hg19.1905.merged <- NormalizeData(PDX.hg19.1905.merged)
#PDX.hg19.1905.merged <- ScaleData(PDX.hg19.1905.merged)
PDX.hg19.1905.merged <- FindVariableFeatures(PDX.hg19.1905.merged,selection.method = "vst",nfeatures = 200)
#hvg.hg19.1905 <- rownames(x = head(x = PDX.hg19.1905.merged@hvg.info, n = 200))
PDX.hg19.1905.merged@meta.data[,"batch"] <- "2019_05"

PDX.hg19.1906.merged <- NormalizeData(PDX.hg19.1906.merged)
#PDX.hg19.1906.merged <- ScaleData(PDX.hg19.1906.merged)
PDX.hg19.1906.merged <- FindVariableFeatures(PDX.hg19.1906.merged,selection.method = "vst",nfeatures = 200)
#hvg.hg19.1906 <- rownames(x = head(x = PDX.hg19.1906.merged@hvg.info, n = 200))
PDX.hg19.1906.merged@meta.data[,"batch"] <- "2019_06"

# hvg.union <- Reduce("union",list(hvg.hg19.1809,
#                                  hvg.hg19.1810,
#                                  hvg.hg19.1811,
#                                  hvg.hg19.1905,
#                                  hvg.hg19.1906))
# 

PDX.hg19.exp.list <- list(PDX.hg19.1809.merged,
                          #PDX.hg19.1810.merged,
                          PDX.hg19.1811.merged,
                          PDX.hg19.1905.merged,
                          PDX.hg19.1906.merged)

names(PDX.hg19.exp.list) <- c("2018_09","2018_11","2019_05","2019_06")
PDX.hg19.anchors <- FindIntegrationAnchors(PDX.hg19.exp.list)
PDX.hg19.combined <- IntegrateData(PDX.hg19.anchors)

DefaultAssay(PDX.hg19.combined) <- "integrated"
PDX.hg19.combined <- ScaleData(PDX.hg19.combined)
PDX.hg19.combined <- RunPCA(PDX.hg19.combined,npcs=20)
PDX.hg19.combined <- RunUMAP(PDX.hg19.combined,reduction = "pca",dims=1:20)
PDX.hg19.combined <- FindNeighbors(PDX.hg19.combined, reduction = "pca", dims = 1:20)
PDX.hg19.combined <- FindClusters(PDX.hg19.combined, resolution = 0.75)
p1 <- DimPlot(PDX.hg19.combined, reduction = "umap", group.by = "batch")
p2 <- DimPlot(PDX.hg19.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
PDX.hg19.combined <- RunTSNE(PDX.hg19.combined,reduction = "pca",dims=1:20)
p1 <- DimPlot(PDX.hg19.combined, reduction = "tsne", group.by = "batch")
p2 <- DimPlot(PDX.hg19.combined, reduction = "tsne", label = TRUE)
plot_grid(p1, p2)

all.metadata <- PDX.hg19.combined@meta.data
all.metadata[,"sample"] <- sapply(rownames(all.metadata),function(x){
  if (grepl("._",x,fixed=T)) {
    substr(x,1,nchar(x)-20)
  } else {
    substr(x,1,nchar(x)-19)
  }
})
PDX.hg19.combined@meta.data <- all.metadata
# p1 <- DimPlot(PDX.hg19.combined, reduction = "tsne", group.by = "sample")
# plot(p1)
Idents(PDX.hg19.combined) <- "sample"
PDX.hg19.hap <- subset(PDX.hg19.combined,idents=names(table(all.metadata$sample))[grep("hap",
                                                                                names(table(all.metadata$sample)),
                                                                                fixed=T)])
PDX.hg19.dip <- subset(PDX.hg19.combined,idents=names(table(all.metadata$sample))[grep("dip",
                                                                                       names(table(all.metadata$sample)),
                                                                                       fixed=T)])
p1 <- DimPlot(PDX.hg19.hap, reduction = "tsne", group.by = "sample") + xlim(-10,7) + ylim(-6,6)
p2 <- DimPlot(PDX.hg19.dip, reduction = "tsne", group.by = "sample")
plot_grid(p1, p2)
plot_grid(p1)
p1 <- DimPlot(PDX.hg19.combined, reduction = "umap", group.by = "sample") + xlim(-10,7) + ylim(-6,6)
p1 <- DimPlot(PDX.hg19.combined, reduction = "umap", group.by = "seurat_clusters") + xlim(-10,7) + ylim(-6,6)

plot_grid(p1)
PDX.hg19.hap.BM <- subset(PDX.hg19.hap,idents=names(table(PDX.hg19.hap@meta.data$sample))[grep("BM",
                                                                                       names(table(PDX.hg19.hap@meta.data$sample)),
                                                                                       fixed=T)])
PDX.hg19.dip.BM <- subset(PDX.hg19.dip,idents=names(table(PDX.hg19.dip@meta.data$sample))[grep("BM",
                                                                                       names(table(PDX.hg19.dip@meta.data$sample)),
                                                                                       fixed=T)])
PDX.hg19.hap.SP <- subset(PDX.hg19.hap,idents=names(table(PDX.hg19.hap@meta.data$sample))[grep("SP",
                                                                                       names(table(PDX.hg19.hap@meta.data$sample)),
                                                                                       fixed=T)])
PDX.hg19.dip.SP <- subset(PDX.hg19.dip,idents=names(table(PDX.hg19.dip@meta.data$sample))[grep("SP",
                                                                                       names(table(PDX.hg19.dip@meta.data$sample)),
                                                                                       fixed=T)])
p1 <- DimPlot(PDX.hg19.hap.BM, reduction = "umap", group.by = "sample") + xlim(-10,7) + ylim(-6,6)
p2 <- DimPlot(PDX.hg19.dip.BM, reduction = "umap", group.by = "sample") + xlim(-10,7) + ylim(-6,6)
p3 <- DimPlot(PDX.hg19.hap.SP, reduction = "umap", group.by = "sample") + xlim(-10,7) + ylim(-6,6)
p4 <- DimPlot(PDX.hg19.dip.SP, reduction = "umap", group.by = "sample") + xlim(-10,7) + ylim(-6,6)
plot_grid(p1,p2,p3,p4)

allSamples <- names(table(all.metadata$sample))
nCluster <- max(as.numeric(as.character(all.metadata$seurat_clusters)))+1
samples.list <- list()
clusterProportions <- lapply(allSamples,function(s){
  metadata <- all.metadata[grep(s,rownames(all.metadata),fixed=T),]
  sampleSubset <- subset(PDX.hg19.combined,cells = rownames(all.metadata)[grep(s,rownames(all.metadata),fixed=T)])
  samples.list[[s]] <<- sampleSubset
  lookupTable <- rep(0,nCluster)
  names(lookupTable) <- as.character(c(0:(nCluster-1)))
  sampleClusters <- table(metadata$seurat_clusters)
  lookupTable[names(sampleClusters)] <- sampleClusters
  #lookupTable <- lookupTable/sum(lookupTable)
  lookupTable
})
names(clusterProportions) <- allSamples
clusterProportions <- do.call(rbind.data.frame,clusterProportions)
rownames(clusterProportions) <- allSamples
colnames(clusterProportions) <- as.character(0:(nCluster-1))
write.table(clusterProportions/rowSums(clusterProportions),
            "clusterProportions_hg19.txt",sep='\t',quote=F)


comparisonPairs.BM <- list(c("WT.dip.BM","WTVeh.dip.BM"),
                           c("WTVeh.dip.BM","WTTX.dip.BM"),
                           c("WTVeh.hap.BM","WT.hap.BM"),
                           c("WTVeh.hap.BM","WTTX.hap.BM"))
#library(RColorBrewer)
#colvec <- brewer.pal(18,"Spectral")
library(choosecolor)
# Mypalette<-palette.picker(n=4)
# colvec <- Mypalette(9)
# write.table(colvec,"9colorshg19.txt",sep='\t',row.names=F,col.names=F,quote=F)
colvec <- read.delim("9colorshg19.txt",sep='\t',header=F,check.names=F)
colvec <- as.character(colvec[,1])
names(colvec) <- as.character(0:(nCluster-1))
clusterProportions.diff <- lapply(comparisonPairs.BM,function(cp){
  s1 <- cp[1]
  s2 <- cp[2]
  top2clusters <- rev(sort(abs(clusterProportions[[s1]]/sum(clusterProportions[[s1]])-
                                 clusterProportions[[s2]]/sum(clusterProportions[[s2]]))))
  top2clusters[1:2]
})
plot.list <- lapply(comparisonPairs.BM[[1]],function(s){
  xxx <- table(samples.list[[s]]@meta.data$seurat_clusters)
  availClusters <- names(xxx)[xxx>0]
  p1 <- DimPlot(PDX.hg19.combined, reduction = "umap", group.by = "seurat_clusters",
                cells = rownames(all.metadata)[grep(s,rownames(all.metadata),fixed=T)],
                cols = colvec[availClusters]) + xlim(-10,6) + ylim(-6,6)
  p1
})
plot_grid(plotlist = plot.list)

# comparisonPairs.SP <- list(c("WT.dip.SP","WTVeh.dip.SP"),
#                            c("WTVeh.dip.SP","WTTX.dip.SP"),
#                            c("WTVeh.hap.SP","WT.hap.SP"),
#                            c("WTVeh.hap.SP","WTTX.hap.SP"))
Idents(PDX.hg19.combined) <- "seurat_clusters"
hg19.marker.genes <- FindAllMarkers(PDX.hg19.combined)
hg19.marker.genes$gene <- gsub("hg19-","",hg19.marker.genes$gene,fixed=T)
write.table(hg19.marker.genes,"hg19.markerGenes.txt",sep='\t',quote=F)
PDX.hg19.hashtag.markers <- hg19.marker.genes
PDX.hg19.hashtag.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- PDX.hg19.hashtag.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
clusters_to_plot <- c(3,4,6)

DoHeatmap(object = PDX.hg19.combined, features = paste0("hg19-",top10$gene))

save(PDX.hg19.combined,file="PDX.hg19.combined.seurat.Robj")

# Load cluster marker genes

DefaultAssay(PDX.hg19.combined) <- "RNA"
PDX.hg19.combined <- ScaleData(PDX.hg19.combined)
hg19_markers <- read.delim("hg19.markerGenes.txt",header=T,check.names=F)
clusters_to_plot <- c(3,4,6)
sampleGroups <- c("WTVeh.hap.BM","WTTX.hap.BM","WTVeh.dip.BM","WTTX.dip.BM")
genes_to_plot <- rownames(hg19_markers)[hg19_markers$cluster %in% clusters_to_plot]
DotPlot(subset(PDX.hg19.combined, subset = sample %in% sampleGroups), features = genes_to_plot)

Idents(PDX.hg19.combined) <- "sample"
top10_CoI <- PDX.hg19.hashtag.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) %>%
  filter(cluster == 4)
temp <- subset(PDX.hg19.combined, subset = sample %in% sampleGroups)
temp@active.ident <- factor(temp@active.ident,
                            levels = c("WTVeh.hap.BM","WTTX.hap.BM","WTVeh.dip.BM","WTTX.dip.BM"))

DotPlot(temp, features = paste0("hg19-",top10_CoI$gene),
        cols = c("goldenrod1","navyblue")) + RotatedAxis()
DoHeatmap(subset(PDX.hg19.combined, subset = sample %in% sampleGroups[1:2]), features = paste0("hg19-",top10_CoI$gene))






