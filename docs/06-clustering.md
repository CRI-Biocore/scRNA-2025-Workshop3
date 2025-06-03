# Section 6: Scaling, and Clustering for Integrated Multi-sample Data 

Similarly to the analysis performed on individual samples previously as shown in the below workflow

![](./images/seurat_norm_clustering_workflow.png)

**The integrated object is ready for downstream analysis steps, as outlined below:**

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

> **Note:** The analysis steps should complete within several minutes.

