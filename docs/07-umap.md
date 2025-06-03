# Section 7: Clustering Results Visualizations with UMAP/TSNE

This section provides code to visualize integrated and clustered single-cell data using UMAP and t-SNE, 
along with violin plots for metadata-based quality control.

## Dimensionality Reduction Plots

```r
seuratObj.int <- readRDS(file = '/gpfs/data/biocore-workshop/scRNA-seq_2025_workshop3/testData/data2_seurat/part2_demo.rds')
# UMAP and t-SNE plots to visualize clusters
pdf(file = file.path(getwd(), 'integratedSamp_umap.pdf'), width = 4, height = 3)
DimPlot(seuratObj.int, reduction = "umap")
dev.off()

pdf(file = file.path(getwd(), 'integratedSamp_tsne.pdf'), width = 4, height = 3)
DimPlot(seuratObj.int, reduction = "tsne")
dev.off()

# Split view by experimental condition
pdf(file = file.path(getwd(), 'integratedSamp_umap_sampSep.pdf'), width = 4, height = 3)
DimPlot(seuratObj.int, reduction = "umap", split.by = 'expCond2')
dev.off()

# Grouped and split view of UMAP
pdf(file = file.path(getwd(), 'integratedSamp_umap_sampSep2.pdf'), width = 4, height = 3)
DimPlot(seuratObj.int, reduction = "umap", group.by = 'expCond3', split.by = 'expCond1')
dev.off()
```

<!--
## Violin Plots on Metadata Information

```R
# Violin plots for general QC metrics
VlnPlot(seuratObj.int, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "rRNA.content"), 
        pt.size = 0, ncol = 4, group.by = 'orig.ident')

# Violin plots grouped by experimental condition
VlnPlot(seuratObj.int, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "rRNA.content"), 
        pt.size = 0, ncol = 4, group.by = 'expCond3')

```
-->