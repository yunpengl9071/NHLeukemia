rm(list=ls())

setwd("~/Data/RegevLab/NHLeukemia/KBM7_10X_020718/")

#source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")
library(cellrangerRkit)

# Read in cellranger outputs
cellranger_pipestance_path <- "haploid"
haploid.data <- load_cellranger_matrix(cellranger_pipestance_path)
cellranger_pipestance_path <- "diploid"
diploid.data <- load_cellranger_matrix(cellranger_pipestance_path)
#KBM7.list <- list(haploid.data,diploid.data)
#KBM7.data <- concatenate_gene_bc_matrices(KBM7.list)
geneTable <- fData(haploid.data)
# haploid.counts <- data.matrix(exprs(KBM7.data))[,1:ncol(exprs(haploid.data))]
# diploid.counts <- data.matrix(exprs(KBM7.data))[,(ncol(exprs(haploid.data))+1):
#                                                   (ncol(exprs(haploid.data))+ncol(exprs(diploid.data)))]
haploid.counts <- data.matrix(exprs(haploid.data))
diploid.counts <- data.matrix(exprs(diploid.data))
rownames(haploid.counts) <- as.character(geneTable[rownames(haploid.counts),"symbol"])
rownames(diploid.counts) <- as.character(geneTable[rownames(diploid.counts),"symbol"])
write.table(haploid.counts, "haploid_KBM7_counts.csv", sep=',', quote=F)
write.table(diploid.counts, "diploid_KBM7_counts.csv", sep=',', quote=F)

# Read in MAGIC imputation results
haploid.exprs <- read.delim("MAGIC_data_hap.txt",
                            header=T,row.names=1,check.names=F)
haploid.exprs <- data.matrix(haploid.exprs)
diploid.exprs <- read.delim("MAGIC_data_dip.txt",
                            header=T,row.names=1,check.names=F)
diploid.exprs <- data.matrix(diploid.exprs)
haploid.exprs <- t(haploid.exprs)

geneTable <- read.delim('haploid/outs/filtered_gene_bc_matrices/GRCh38/genes.tsv',
                        header=F,check.names=F)
rownames(geneTable) <- geneTable[,1]

rownames(haploid.exprs) <- gsub("x","",rownames(haploid.exprs),fixed=T)

rownames(haploid.exprs) <- as.character(geneTable[,2])

diploid.exprs <- t(diploid.exprs)
rownames(diploid.exprs) <- as.character(geneTable[,2])
gc()

nHaploidCells <- ncol(haploid.exprs)
nDiploidCells <- ncol(diploid.exprs)

hap.nonZeroSD <- haploid.exprs[apply(haploid.exprs,1,sd)!=0,]
dip.nonZeroSD <- diploid.exprs[apply(diploid.exprs,1,sd)!=0,]
gc()
hap.nonZeroSD <- hap.nonZeroSD[log2(rowMeans(hap.nonZeroSD))>-10,]
dip.nonZeroSD <- dip.nonZeroSD[log2(rowMeans(dip.nonZeroSD))>-10,]
gc()

hap.nz.scaled <- t(scale(t(log2(hap.nonZeroSD+1))))
dip.nz.scaled <- t(scale(t(log2(dip.nonZeroSD+1))))

cellCycleStages <- c("MG1","G1S","S","G2","G2M")
cellCycleSig <- lapply(cellCycleStages,function(x){
  fn <- paste0(x,".txt")
  as.character(read.delim(fn,header=F,check.names=F)[,1])
})
names(cellCycleSig) <- cellCycleStages

hap.sigGenes <- lapply(cellCycleSig,function(x){
  availGenes <- intersect(x,rownames(hap.nonZeroSD))
})
names(hap.sigGenes) <- cellCycleStages
dip.sigGenes <- lapply(cellCycleSig,function(x){
  availGenes <- intersect(x,rownames(dip.nonZeroSD))
})
names(dip.sigGenes) <- cellCycleStages
all.sigGenes <- lapply(cellCycleStages,function(x){
  intersect(hap.sigGenes[[x]],
            dip.sigGenes[[x]])
})
names(all.sigGenes) <- cellCycleStages
hap.metagenes <- lapply(cellCycleStages,function(x){
  availGenes <- hap.sigGenes[[x]]
  colMeans(log2(hap.nonZeroSD[availGenes,]+1))
})
names(hap.metagenes) <- cellCycleStages
dip.metagenes <- lapply(cellCycleStages,function(x){
  availGenes <- dip.sigGenes[[x]]
  colMeans(log2(dip.nonZeroSD[availGenes,]+1))
})
names(dip.metagenes) <- cellCycleStages

corrThreshold <- 0.5
hapSig <- "hap"
dipSig <- "dip"

param_string <- paste0("cor",as.character(corrThreshold),
                       "_hapSig_",hapSig,
                       "_dipSig_",dipSig)
hap.coreSigGenes <- lapply(cellCycleStages,function(x){
  sigGenes <- hap.sigGenes[[x]]
  sigGenesCor <- sapply(sigGenes,function(g){
    cor(log2(hap.nonZeroSD[g,]+1),hap.metagenes[[x]])
  })
  coreSigGenes <- sigGenes[sigGenesCor > corrThreshold]
})
names(hap.coreSigGenes) <- cellCycleStages
sapply(hap.coreSigGenes,length)
dip.coreSigGenes <- lapply(cellCycleStages,function(x){
  sigGenes <- dip.sigGenes[[x]]
  sigGenesCor <- sapply(sigGenes,function(g){
    cor(log2(dip.nonZeroSD[g,]+1),dip.metagenes[[x]])
  })
  coreSigGenes <- sigGenes[sigGenesCor > corrThreshold]
})
names(dip.coreSigGenes) <- cellCycleStages
sapply(dip.coreSigGenes,length)

all.coreSigGenes <- lapply(cellCycleStages,function(x){
  sigGenes <- all.sigGenes[[x]]
  metagene <- c(colMeans(hap.nz.scaled[sigGenes,]),colMeans(dip.nz.scaled[sigGenes,]))
  sigGenesCor <- sapply(sigGenes,function(g){
    cor(c(hap.nz.scaled[g,],
          dip.nz.scaled[g,]),
        metagene)
  })
  coreSigGenes <- sigGenes[sigGenesCor > corrThreshold]
})
names(all.coreSigGenes) <- cellCycleStages
sapply(all.coreSigGenes,length)

# Generate metagene scores using core signature genes and plot cell cycle scores

hap.coreMetagenes <- lapply(cellCycleStages,function(x){
  if (hapSig == "hap") {
    colMeans((hap.nz.scaled[hap.coreSigGenes[[x]],,drop=F]))
  } else if (hapSig == "dip") {
    colMeans((hap.nz.scaled[intersect(dip.coreSigGenes[[x]],
                                      rownames(hap.nz.scaled)),,drop=F]))
  } else if (hapSig == "all") {
    colMeans((hap.nz.scaled[all.coreSigGenes[[x]],,drop=F]))
  }
})
names(hap.coreMetagenes) <- cellCycleStages
hap.coreMetagenes.df <- do.call(cbind.data.frame,hap.coreMetagenes)
hap.coreMetagenes.mat <- t(scale(t(data.matrix(hap.coreMetagenes.df))))

dip.coreMetagenes <- lapply(cellCycleStages,function(x){
  if (dipSig == "hap") {
    colMeans((dip.nz.scaled[intersect(hap.coreSigGenes[[x]],
                                      rownames(dip.nz.scaled)),,drop=F]))
  } else if (dipSig == "dip") {
    colMeans((dip.nz.scaled[dip.coreSigGenes[[x]],,drop=F]))
  } else if (dipSig == "all") {
    colMeans((dip.nz.scaled[all.coreSigGenes[[x]],,drop=F]))
  }
})
names(dip.coreMetagenes) <- cellCycleStages
dip.coreMetagenes.df <- do.call(cbind.data.frame,dip.coreMetagenes)
dip.coreMetagenes.mat <- t(scale(t(data.matrix(dip.coreMetagenes.df))))

# Assign cell cycle stages by looking at the highest core for each cell
hap.cellCycleStages <- apply(hap.coreMetagenes.mat,1,function(x){
  cellCycleStages[which(x==max(x))[1]]
})
dip.cellCycleStages <- apply(dip.coreMetagenes.mat,1,function(x){
  cellCycleStages[which(x==max(x))[1]]
})
hap.table <- table(hap.cellCycleStages)/length(hap.cellCycleStages)
dip.table <- table(dip.cellCycleStages)/length(dip.cellCycleStages)

# for (stage in cellCycleStages) {
#   print(stage)
#   hap.fn <- paste0("haploid_",stage,"_MAGIC_exprs_nz_scaled.txt")
#   write.table(hap.nz.scaled[,hap.cellCycleStages==stage],hap.fn,sep='\t',
#               quote=F)
#   dip.fn <- paste0("diploid_",stage,"_MAGIC_exprs_nz_scaled.txt")
#   write.table(dip.nz.scaled[,dip.cellCycleStages==stage],dip.fn,sep='\t',
#               quote=F)
# }

hap.col <- rep("violet",length(hap.cellCycleStages))
names(hap.col) <- names(hap.cellCycleStages)
hap.col[hap.cellCycleStages=="S"] <- "lightblue"
hap.col[hap.cellCycleStages=="G2"] <- "navyblue"
hap.col[hap.cellCycleStages=="G2M"] <- "orange"
hap.col[hap.cellCycleStages=="MG1"] <- "pink"
dip.col <- rep("violet",length(dip.cellCycleStages))
names(dip.col) <- names(dip.cellCycleStages)
dip.col[dip.cellCycleStages=="S"] <- "lightblue"
dip.col[dip.cellCycleStages=="G2"] <- "navyblue"
dip.col[dip.cellCycleStages=="G2M"] <- "orange"
dip.col[dip.cellCycleStages=="MG1"] <- "pink"
hap.noG2MG1 <- which(hap.cellCycleStages!= "G2" & hap.cellCycleStages != "MG1")
dip.noG2MG1 <- which(dip.cellCycleStages!= "G2" & dip.cellCycleStages != "MG1")
hap.noMG1 <- which(hap.cellCycleStages != "MG1")
dip.noMG1 <- which(dip.cellCycleStages != "MG1")
 
hap.fn <- paste0("haploid_MAGIC_t_3_cellCycle_",param_string,".pdf")
dip.fn <- paste0("diploid_MAGIC_t_3_cellCycle_",param_string,".pdf")

pdf(file = hap.fn,width = 7,height = 8)
plot(hap.coreMetagenes.mat[,"G1S"],hap.coreMetagenes.mat[,"G2M"],pch=16,col=hap.col,
     xlab = "Normalized G1S score", ylab = "Normalized G2M score")
dev.off()
pdf(file = dip.fn,width = 7,height = 8)
plot(dip.coreMetagenes.mat[,"G1S"],dip.coreMetagenes.mat[,"G2M"],pch=16,col=dip.col,
     xlab = "Normalized G1S score", ylab = "Normalized G2M score")
dev.off()

hap.fn <- paste0("haploid_MAGIC_t_3_cellCycle_",param_string,"_MG1.pdf")
dip.fn <- paste0("diploid_MAGIC_t_3_cellCycle_",param_string,"_MG1.pdf")

pdf(file = hap.fn,width = 7,height = 8)
plot(hap.coreMetagenes.mat[hap.cellCycleStages=="MG1",
                           "G1S"],
     hap.coreMetagenes.mat[hap.cellCycleStages=="MG1",
                           "G2M"],pch=16,col=hap.col[hap.cellCycleStages=="MG1"],
     xlab = "Normalized G1S score", ylab = "Normalized G2M score")
dev.off()
pdf(file = dip.fn,width = 7,height = 8)
plot(dip.coreMetagenes.mat[dip.cellCycleStages=="MG1",
                           "G1S"],
     dip.coreMetagenes.mat[dip.cellCycleStages=="MG1",
                           "G2M"],pch=16,col=dip.col[dip.cellCycleStages=="MG1"],
     xlab = "Normalized G1S score", ylab = "Normalized G2M score")
dev.off()


stage <- "G2M"
stage_color <- "navyblue"
hap.data <- cbind.data.frame(hap.coreMetagenes.mat[hap.cellCycleStages==stage,
                                                   "G1S"],
                             hap.coreMetagenes.mat[hap.cellCycleStages==stage,
                                                   "G2M"])
colnames(hap.data) <- c("G1S","G2M")
# ggplot(hap.data, aes(x=G1S,y=G2M)) + geom_bin2d(bins=20) + theme_bw() + 
#   scale_fill_gradientn(limits=c(0,7), breaks=seq(0, 40, by=10), colours=c("white","pink"))
# ggplot(hap.data, aes(x=G1S,y=G2M)) + geom_density_2d(bins=10) + theme_bw() 


dip.data <- cbind.data.frame(dip.coreMetagenes.mat[dip.cellCycleStages==stage,
                                                   "G1S"],
                             dip.coreMetagenes.mat[dip.cellCycleStages==stage,
                                                   "G2M"])
colnames(dip.data) <- c("G1S","G2M")

ggplot(hap.data, aes(x=G1S,y=G2M)) + geom_density_2d(bins=10) + theme_bw() 
#+ scale_fill_gradientn(limits=c(0,15), breaks=seq(0, 40, by=10), colours=c("white","pink"))
library(spatstat)
pppo=ppp(x=hap.data$G1S,y=hap.data$G2M,window = owin(c(-2,2),c(-2,2)))
den=density(pppo,kernel="gaussian",edge=T,diggle=T,adjust=0.6)
plot(den,main='haploid',col=colorRampPalette(c("white",stage_color))(10),xlim=c(-2,2),
     ylim=c(-2,2))

ggplot(dip.data, aes(x=G1S,y=G2M)) + geom_density_2d(bins=10) + theme_bw() 
#+ scale_fill_gradientn(limits=c(0,15), breaks=seq(0, 40, by=10), colours=c("white","pink"))
library(spatstat)
pppo=ppp(x=dip.data$G1S,y=dip.data$G2M,window = owin(c(-2,2),c(-2,2)))
den=density(pppo,kernel="gaussian",edge=T,diggle=T,adjust=0.6)
plot(den,main='diploid',col=colorRampPalette(c("white",stage_color))(10),xlim=c(-2,2),
     ylim=c(-2,2))

# Fit an ellipse to the G1S-G2M plot and project data onto the fitted pseudotrajectory
# Use diploid data as reference standard.
hap.core2D <- hap.coreMetagenes.mat[,c("G1S","G2M")]
dip.core2D <- dip.coreMetagenes.mat[,c("G1S","G2M")]
source("ellipseFit.R")
dip.eFit <- fit.ellipse(dip.core2D[,"G1S"],dip.core2D[,"G2M"])
dip.e <- get.ellipse(dip.eFit)
plot(dip.core2D,pch=16,col="grey")
lines(dip.e,lwd=5)

# Use an epsilon-neighborhood circle method to average points and smooth
# out expression data along the fitted pseudotrajectory.
epsilon <- 0.5
goi <- ""
goi.set <- c("CCND1","CCNE1","CCNA1","CCNB1","RAD51","RAD51B")
goi.set <- c("RAD51")
for (goi in goi.set) {
  fn <- paste0(goi,"_elliptical_pseudotime.pdf")
  pdf(file = fn,width = 6,height = 5)
  dip.e.intervals <- dip.e[seq(1,360,3),]
  dip.goi.smthExprs <- t(apply(dip.e.intervals,1,function(x){
    epsilon.points <- which(apply(dip.core2D,1,function(y){
      ((x[1]-y[1])^2+(x[2]-y[2])^2)^0.5
    })<epsilon)
    avg <- mean(dip.nonZeroSD[goi,epsilon.points])
    sd <- sd(dip.nonZeroSD[goi,epsilon.points])
    n <- length(epsilon.points)
    c(avg,sd,n)
  }))
  plot(dip.goi.smthExprs[rev(c(41:120,1:40)),1],pch=".",ylim=c(0.1,0.35),main=goi,
       xlab="pseudotime",ylab="average expression")
  lines(dip.goi.smthExprs[rev(c(41:120,1:40)),1],lwd=5,col="orange")
  CI.y <- c((dip.goi.smthExprs[,1] - 1.96*dip.goi.smthExprs[,2]/(dip.goi.smthExprs[,3])^0.5)[rev(c(41:120,1:40))],
            (dip.goi.smthExprs[,1] + 1.96*dip.goi.smthExprs[,2]/(dip.goi.smthExprs[,3])^0.5)[c(41:120,1:40)])
  CI.x <- c(1:120,120:1)
  polygon(CI.x,CI.y,col=adjustcolor("orange",alpha.f = 0.2),
          lty='blank')
  
  hap.goi.smthExprs <- t(apply(dip.e.intervals,1,function(x){
    epsilon.points <- which(apply(hap.core2D,1,function(y){
      ((x[1]-y[1])^2+(x[2]-y[2])^2)^0.5
    })<epsilon)
    avg <- mean(2*hap.nonZeroSD[goi,epsilon.points])
    sd <- sd(2*hap.nonZeroSD[goi,epsilon.points])
    n <- length(epsilon.points)
    c(avg,sd,n)
  }))
  #plot(hap.goi.smthExprs,pch=".")
  lines(hap.goi.smthExprs[rev(c(41:120,1:40)),1],lwd=5,col="navyblue")
  CI.y <- c((hap.goi.smthExprs[,1] - 1.96*hap.goi.smthExprs[,2]/(hap.goi.smthExprs[,3])^0.5)[rev(c(41:120,1:40))],
            (hap.goi.smthExprs[,1] + 1.96*hap.goi.smthExprs[,2]/(hap.goi.smthExprs[,3])^0.5)[c(41:120,1:40)])
  CI.x <- c(1:120,120:1)
  polygon(CI.x,CI.y,col=adjustcolor("navyblue",alpha.f = 0.2),
          lty='blank')
  dev.off()
}


hap.sigExprsByStage <- lapply(cellCycleStages,function(s){
  sigGenes <- hap.coreSigGenes[[s]]
  mat <- hap.nz.scaled[sigGenes,hap.cellCycleStages==s]
  hclust.row <- hclust(dist(mat),method = "complete")
  #hclust.col <- hclust(dist(t(mat)),method = "complete")
  #colIdx.list[[s]] <- hclust.col$order
  # Determine column order using a measurement of how 'pure' the cell is in
  # terms of metagene expression
  # Now using just expression of the label metagene
  mat.all <- hap.coreMetagenes.mat[hap.cellCycleStages==s,]
  idx <- which(cellCycleStages==s)
  purity <- apply(mat.all,1,function(x){
    x[idx]
  })
  colOrder <- rev(order(purity))
  list(hap.nz.scaled[sigGenes,][hclust.row$order,],colOrder)
})
colIdx.list <- lapply(hap.sigExprsByStage,function(x){x[[2]]})
names(colIdx.list) <- cellCycleStages
hap.sigExprsByStage <- lapply(hap.sigExprsByStage,function(x){x[[1]]})
hap.sigLen <- sapply(hap.sigExprsByStage,nrow)
hap.sigExprsByStage.mat <- do.call(rbind.data.frame,hap.sigExprsByStage)
hap.sigExprsByStage.mat <- data.matrix(hap.sigExprsByStage.mat)
hap.colID <- unlist(lapply(cellCycleStages,function(s){
  which(hap.cellCycleStages==s)[colIdx.list[[s]]]
}))

dip.sigExprsByStage <- lapply(cellCycleStages,function(s){
  sigGenes <- dip.coreSigGenes[[s]]
  mat <- dip.nz.scaled[sigGenes,dip.cellCycleStages==s]
  hclust.row <- hclust(dist(mat),method = "complete")
  #hclust.col <- hclust(dist(t(mat)),method = "complete")
  #colIdx.list[[s]] <- hclust.col$order
  # Determine column order using a measurement of how 'pure' the cell is in
  # terms of metagene expression
  # Now using just expression of the label metagene
  mat.all <- dip.coreMetagenes.mat[dip.cellCycleStages==s,]
  idx <- which(cellCycleStages==s)
  purity <- apply(mat.all,1,function(x){
    x[idx]
  })
  colOrder <- rev(order(purity))
  list(dip.nz.scaled[sigGenes,][hclust.row$order,],colOrder)
})
colIdx.list <- lapply(dip.sigExprsByStage,function(x){x[[2]]})
names(colIdx.list) <- cellCycleStages
dip.sigExprsByStage <- lapply(dip.sigExprsByStage,function(x){x[[1]]})
dip.sigLen <- sapply(dip.sigExprsByStage,nrow)
dip.sigExprsByStage.mat <- do.call(rbind.data.frame,dip.sigExprsByStage)
dip.sigExprsByStage.mat <- data.matrix(dip.sigExprsByStage.mat)
dip.colID <- unlist(lapply(cellCycleStages,function(s){
  which(dip.cellCycleStages==s)[colIdx.list[[s]]]
}))


library(gplots)

# Now show expression of cyclins etc.

cyclins <- c("CCND1","CCNE1","CCNA1","CCNA2","CCNB1","CCNB2")

colIdx.list <- lapply(cellCycleStages,function(s){
  mat.all <- hap.coreMetagenes.mat[hap.cellCycleStages==s,]
  idx <- which(cellCycleStages==s)
  purity <- apply(mat.all,1,function(x){
    x[idx]
  })
  colOrder <- rev(order(purity))
})
names(colIdx.list) <- cellCycleStages
hap.colID <- unlist(lapply(cellCycleStages,function(s){
  which(hap.cellCycleStages==s)
}))
hap.cyclinExprs <- hap.nz.scaled[cyclins,hap.colID]

pdf(file=paste0("hap_cyclins_",param_string,".pdf"),width = 10,height = 6)
hmRes <- heatmap3(hap.cyclinExprs[rev(1:nrow(hap.cyclinExprs)),],
                  Colv=NA,
                  Rowv=NA,
                  scale="none",
                  density.info="none", trace="none",
                  col=myheatmapfun(75),
                  reorderfun = function(d,w) { d },
                  cexCol = 0.5,
                  srtCol = 45,
                  cexRow = 0.7,
                  offsetCol = -0.5,
                  #RowSideColors = rscVec,
                  dendrogram="none",
                  key.title=NA,
                  key.xlab=NA,
                  breaks=seq(-2,2,length.out=76),
                  #margins=c(5,5),
                  lhei=c(1,8),
                  labCol="")
dev.off()

colIdx.list <- lapply(cellCycleStages,function(s){
  mat.all <- dip.coreMetagenes.mat[dip.cellCycleStages==s,]
  idx <- which(cellCycleStages==s)
  purity <- apply(mat.all,1,function(x){
    x[idx]
  })
  colOrder <- rev(order(purity))
})
names(colIdx.list) <- cellCycleStages
dip.colID <- unlist(lapply(cellCycleStages,function(s){
  which(dip.cellCycleStages==s)[colIdx.list[[s]]]
}))
dip.cyclinExprs <- dip.nz.scaled[cyclins,dip.colID]

pdf(file=paste0("dip_cyclins_",param_string,".pdf"),width = 10,height = 6)
hmRes <- heatmap3(dip.cyclinExprs[rev(1:nrow(dip.cyclinExprs)),],
                  Colv=NA,
                  Rowv=NA,
                  scale="none",
                  density.info="none", trace="none",
                  col=myheatmapfun(75),
                  reorderfun = function(d,w) { d },
                  cexCol = 0.5,
                  srtCol = 45,
                  cexRow = 0.7,
                  offsetCol = -0.5,
                  #RowSideColors = rscVec,
                  dendrogram="none",
                  key.title=NA,
                  key.xlab=NA,
                  breaks=seq(-2,2,length.out=76),
                  #margins=c(5,5),
                  lhei=c(1,8),
                  labCol="")
dev.off()

# Also generate scatterplots with color gradients

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

for (goi in cyclins) {
  goi.exprs <- hap.nz.scaled[goi,]
  crp <- colorRampPalette((c(
    rgb(255,255,255,maxColorValue=255),
    rgb(234,145,118,maxColorValue=255),
    rgb(102,3,32,maxColorValue=255))))
  ii <- cut(goi.exprs, breaks = length(goi.exprs), 
            include.lowest = TRUE)
  colors <- crp(length(goi.exprs))[ii]
  colors <- map2color(goi.exprs,crp(length(goi.exprs)),limits=c(-1.5,1.5))
  
  pdf(file=paste0("hap_",goi,"_",param_string,".pdf"),width=7,height=8)
  plot(hap.coreMetagenes.mat[,"G1S"],
       hap.coreMetagenes.mat[,"G2M"],
       pch=16,col=colors,main=goi,
       xlab = "Normalized G1S score", ylab = "Normalized G2M score")
  dev.off()
  
  goi.exprs <- dip.nz.scaled[goi,]
  crp <- colorRampPalette((c(
    rgb(255,255,255,maxColorValue=255),
    rgb(234,145,118,maxColorValue=255),
    rgb(102,3,32,maxColorValue=255))))
  ii <- cut(goi.exprs, breaks = length(goi.exprs), 
            include.lowest = TRUE)
  colors <- crp(length(goi.exprs))[ii]
  colors <- map2color(goi.exprs,crp(length(goi.exprs)),limits=c(-1.5,1.5))
  
  pdf(file=paste0("dip_",goi,"_",param_string,".pdf"),width=7,height=8)
  plot(dip.coreMetagenes.mat[,"G1S"],
       dip.coreMetagenes.mat[,"G2M"],
       pch=16,col=colors,main=goi,
       xlab = "Normalized G1S score", ylab = "Normalized G2M score")
  dev.off()
}

for (s in cellCycleStages) {
  goi.exprs <- hap.coreMetagenes.mat[,s]
  crp <- colorRampPalette((c(
    rgb(255,255,255,maxColorValue=255),
    rgb(234,145,118,maxColorValue=255),
    rgb(102,3,32,maxColorValue=255))))
  ii <- cut(goi.exprs, breaks = length(goi.exprs), 
            include.lowest = TRUE)
  colors <- crp(length(goi.exprs))[ii]
  colors <- map2color(goi.exprs,crp(length(goi.exprs)),limits=c(-1.5,1.5))
  
  pdf(file=paste0("hap_",s,"_metagene_",param_string,".pdf"),width=7,height=8)
  plot(hap.coreMetagenes.mat[,"G1S"],
       hap.coreMetagenes.mat[,"G2M"],
       pch=16,col=colors,main=s,
       xlab = "Normalized G1S score", ylab = "Normalized G2M score")
  dev.off()
  
  goi.exprs <- dip.coreMetagenes.mat[,s]
  crp <- colorRampPalette((c(
    rgb(255,255,255,maxColorValue=255),
    rgb(234,145,118,maxColorValue=255),
    rgb(102,3,32,maxColorValue=255))))
  ii <- cut(goi.exprs, breaks = length(goi.exprs), 
            include.lowest = TRUE)
  colors <- crp(length(goi.exprs))[ii]
  colors <- map2color(goi.exprs,crp(length(goi.exprs)),limits=c(-1.5,1.5))
  
  pdf(file=paste0("dip_",s,"_metagene_",param_string,".pdf"),width=7,height=8)
  plot(dip.coreMetagenes.mat[,"G1S"],
       dip.coreMetagenes.mat[,"G2M"],
       pch=16,col=colors,main=s,
       xlab = "Normalized G1S score", ylab = "Normalized G2M score")
  dev.off()
}

# Now look at disproportionate genes

hap.exprs.byStage <- lapply(cellCycleStages,function(x){
  haploid.counts[,hap.cellCycleStages==x]
})
names(hap.exprs.byStage) <- cellCycleStages
dip.exprs.byStage <- lapply(cellCycleStages,function(x){
  diploid.counts[,dip.cellCycleStages==x]
})

hap.norm.exprs.byStage <- lapply(cellCycleStages,function(x){
  temp <- haploid.exprs[,hap.cellCycleStages==x]
  temp <- apply(temp,2,function(x){
    x/sum(x) * 1E6
  })
  temp
})
names(hap.norm.exprs.byStage) <- cellCycleStages

dip.exprs.byStage <- lapply(cellCycleStages,function(x){
  diploid.counts[,dip.cellCycleStages==x]
})
names(dip.exprs.byStage) <- cellCycleStages

dip.norm.exprs.byStage <- lapply(cellCycleStages,function(x){
  temp <- diploid.exprs[,dip.cellCycleStages==x]
  temp <- apply(temp,2,function(x){
    x/sum(x) * 1E6
  })
  temp
})
names(dip.norm.exprs.byStage) <- cellCycleStages 

dip.hap.ratios.byStage <- lapply(cellCycleStages,function(x){
  rowMeans(dip.exprs.byStage[[x]])/rowMeans(hap.exprs.byStage[[x]])
})
names(dip.hap.ratios.byStage) <- cellCycleStages
for (s in cellCycleStages) {
  pdf(file = paste0(s,"_log2_dip_hap_ratios_hist.pdf"),width=6,height=5)
  hist(log2(dip.hap.ratios.byStage[[s]]),75,xlim=c(-4,4),col='gray',
       xlab = "log2 diploid:haploid expression ratio",main = s)
  dev.off()
}

# Also plot overall ratios
dip.hap.ratios.all <- rowMeans(diploid.counts)/rowMeans(haploid.counts)
pdf(file = paste0("log2_dip_hap_ratios_hist.pdf"),width=6,height=5)
hist(log2(dip.hap.ratios.all)[log2(dip.hap.ratios.all) > -3 & log2(dip.hap.ratios.all) < 4],75,xlim=c(-3,4),col='gray',
     xlab = "log2 diploid:haploid expression ratio",main = "all cells")
dev.off()
# Compare ECDF curves of the ratios in each stage
pdf(file="ECDF.pdf",width=6,height=5)
xxx <- log2(dip.hap.ratios.all)
xxx <- xxx[!is.na(xxx)]
xxx <- xxx[is.finite(xxx)]
plot(ecdf(xxx),xlim=c(-3,4),col="black",lwd=3,
     xlab="log2 diploid:haploid expression ratio (x)",
     ylab="empirical CDF(x)",
     main="")
cellCycleStages.col <- c("pink","violet","lightblue","navyblue","orange")
all.ecdf <- ecdf(xxx)
names(cellCycleStages.col) <- cellCycleStages
for (s in cellCycleStages) {
  xxx <- log2(dip.hap.ratios.byStage[[s]])
  xxx <- xxx[!is.na(xxx)]
  xxx <- xxx[is.finite(xxx)]
  plot(ecdf(xxx),xlim=c(-3,4),col=cellCycleStages.col[s],lwd=3,add=T)
}
dev.off()
ecdf.diff <- sapply(cellCycleStages,function(s){
  xxx <- log2(dip.hap.ratios.byStage[[s]])
  xxx <- xxx[!is.na(xxx)]
  xxx <- xxx[is.finite(xxx)]
  ecdf(xxx)(log2(1.5)) - all.ecdf(log2(1.5))
})

dip.hap.ratios.byStage.df <- do.call(cbind.data.frame,dip.hap.ratios.byStage)
dip.hap.ratios.byStage.mat <- data.matrix(dip.hap.ratios.byStage.df)
rownames(dip.hap.ratios.byStage.mat) <- rownames(dip.exprs.byStage[[1]])
dip.hap.ratios.byStage.mat <- dip.hap.ratios.byStage.mat[complete.cases(dip.hap.ratios.byStage.mat),]
abnormGenes.list <- lapply(1:length(cellCycleStages),function(i){
  abnormGenes <- apply(dip.hap.ratios.byStage.mat,1,function(x){
    if (mean(x[-i])>1.5 & x[i] < 1) {
      T
    } else {
      F
    }
  })
  rownames(dip.hap.ratios.byStage.mat)[abnormGenes]
})
names(abnormGenes.list) <- cellCycleStages
temp <- lapply(cellCycleStages,function(x){
  xx <- abnormGenes.list[[x]]
  AS.LNC.idx <- union(grep("-",xx,fixed=T),grep(".",xx,fixed=T))
  AS.LNC.idx <- union(AS.LNC.idx,grep("LINC",xx,fixed=T))
  write.table(xx[-AS.LNC.idx],paste0("abnormGenes_",param_string,x,"noLINCnoAS_1.5_1.txt"),
              sep="\t",row.names=F,col.names=F,quote=F)
})
temp <- lapply(cellCycleStages,function(x){
  xx <- abnormGenes.list[[x]]
  write.table(xx,paste0("abnormGenes_",param_string,"_",x,"_1.5_1.txt"),
              sep="\t",row.names=F,col.names=F,quote=F)
})
write.table(dip.hap.ratios.byStage.mat,paste0("dip.hap.ratios.byStage_",param_string,".txt"),
            sep="\t",quote=F)

# Look at the genes commonly overexpressed in G2-G2M-MG1
commonSet <- list()
for (stage in cellCycleStages[c(4,5,1)]) {
  diffGenes <- rownames(KBM7.diffExprs.byStage[[stage]])[KBM7.diffExprs.byStage[[stage]]$p_val_adj < 0.001 & 
                                                           KBM7.diffExprs.byStage[[stage]]$avg_logFC > log2(1.5)]
  commonSet[[stage]] <- diffGenes
}
commonSet <- Reduce('intersect',commonSet)

# Try calculating correlations among genes along elliptical manifold

epsilon <- 0.5
dip.e.intervals <- dip.e[seq(1,360,3),]

goi <- "RAD51"

hap.manifoldExprs <- lapply(rownames(hap.nonZeroSD),function(goi){
  hap.goi.smthExprs <- t(apply(dip.e.intervals,1,function(x){
    epsilon.points <- which(apply(hap.core2D,1,function(y){
      ((x[1]-y[1])^2+(x[2]-y[2])^2)^0.5
    })<epsilon)
    avg <- mean(2*hap.nonZeroSD[goi,epsilon.points])
    avg
  }))
  cat(".")
  hap.goi.smthExprs
})

hap.manifoldExprs <- data.matrix(do.call(rbind.data.frame,hap.manifoldExprs))
rownames(hap.manifoldExprs) <- rownames(hap.nonZeroSD)
write.table(hap.manifoldExprs,"hap.manifoldExprs.txt",
            sep='\t',quote=F)

dip.manifoldExprs <- lapply(rownames(dip.nonZeroSD),function(goi){
  dip.goi.smthExprs <- t(apply(dip.e.intervals,1,function(x){
    epsilon.points <- which(apply(dip.core2D,1,function(y){
      ((x[1]-y[1])^2+(x[2]-y[2])^2)^0.5
    })<epsilon)
    avg <- mean(dip.nonZeroSD[goi,epsilon.points])
    avg
  }))
  cat(".")
  dip.goi.smthExprs
})

dip.manifoldExprs <- data.matrix(do.call(rbind.data.frame,dip.manifoldExprs))
rownames(dip.manifoldExprs) <- rownames(dip.nonZeroSD)
write.table(dip.manifoldExprs,"dip.manifoldExprs.txt",
            sep='\t',quote=F)

hap.manifoldCorWithRAD51B <- sapply(rownames(hap.manifoldExprs),function(g){
  cor(hap.manifoldExprs[g,],hap.manifoldExprs["RAD51B",])
})

head(rev(sort(hap.manifoldCorWithRAD51B)),100)

write.table(names(hap.manifoldCorWithRAD51B)[hap.manifoldCorWithRAD51B > 0.5],
            "~/Data/RegevLab/NHLeukemia/KBM7_hap_10X_RAD51B_manifold_signature.txt",
            sep='\t',row.names=F,col.names=F,quote=F)

dip.manifoldCorWithRAD51B <- sapply(rownames(dip.manifoldExprs),function(g){
  cor(dip.manifoldExprs[g,],dip.manifoldExprs["RAD51B",])
})

head(rev(sort(dip.manifoldCorWithRAD51B)),100)

write.table(names(dip.manifoldCorWithRAD51B)[dip.manifoldCorWithRAD51B > 0.5],
            "~/Data/RegevLab/NHLeukemia/KBM7_dip_10X_RAD51B_manifold_signature.txt",
            sep='\t',row.names=F,col.names=F,quote=F)



