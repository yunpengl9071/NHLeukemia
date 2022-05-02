# NHLeukemia
Code used for processing and analyzing 10X scRNA-seq data collected from near-haploid and diploidized KBM7 cell lines and B-ALL cells from patient-derived xenograft (PDX) models.

KBM7_10X_cellRanger_MAGIC_pipeline.R describes the processing of single-cell expression data (10X Genomics cellranger software output), inference of cell cycle stages for each single cell, and estimation of pseudo-time trajectories for cell cycle progression.

DemuxEM_Seurat_analysis_hg19_comparisons.R contains the code for processing single-cell expression data from B-ALL cells in a patient-derived xenograft (PDX) model of near-haploid leukemia. The code combines multiple batches of samples collected across different experiment dates, removes batch effects using the Seurat single cell data integration method (canonical correlation analysis, CCA) and produces a single, integrated Seurat object of normalized gene expression data in the entire study.
