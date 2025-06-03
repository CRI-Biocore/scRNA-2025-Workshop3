# 08. Visualization with UMAP/TSNE

This section provides code to visualize integrated and clustered single-cell data using UMAP and t-SNE, 
along with violin plots for metadata-based quality control.

## Dimensionality Reduction Plots

```r
# UMAP and t-SNE plots to visualize clusters
DimPlot(seuratObj.int, reduction = "umap")
DimPlot(seuratObj.int, reduction = "tsne")

# Split view by experimental condition
DimPlot(seuratObj.int, reduction = "umap", split.by = 'expCond2')

# Grouped and split view of UMAP
DimPlot(seuratObj.int, reduction = "umap", group.by = 'expCond3', split.by = 'expCond1')

```

## Violin Plots on Metadata Information

```R
# Violin plots for general QC metrics
VlnPlot(seuratObj.int, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "rRNA.content"), 
        pt.size = 0, ncol = 4, group.by = 'orig.ident')

# Violin plots grouped by experimental condition
VlnPlot(seuratObj.int, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "rRNA.content"), 
        pt.size = 0, ncol = 4, group.by = 'expCond3')

```
