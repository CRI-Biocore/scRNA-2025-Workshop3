## Clustering Workflow in Seurat

```r
# Step 1: Rescale the integrated Seurat object
seuratObj.int <- ScaleData(object = seuratObjIntegrated)  # Re-scaling if Assay has changed

# Step 2: Principal Component Analysis (PCA)
seuratObj.int <- RunPCA(object = seuratObj.int)

# Step 3: Construct a Shared Nearest Neighbor (SNN) Graph
seuratObj.int <- FindNeighbors(seuratObj.int, reduction = "pca", dims = 1:20)

# Step 4: Identify Clusters using Louvain algorithm
seuratObj.int <- FindClusters(seuratObj.int, resolution = 0.5)

# Step 5: Run UMAP for visualization
seuratObj.int <- RunUMAP(seuratObj.int, reduction = "pca", dims = 1:20)

# Step 6: Run t-SNE for alternate visualization
seuratObj.int <- RunTSNE(object = seuratObj.int, dims = 1:20)
```