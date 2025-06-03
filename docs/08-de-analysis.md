# Section 8: Differential Expression Analysis

Differential expression (DE) analysis identifies genes that are differentially expressed across different conditions, clusters, or populations. This section walks through the process of performing DE analysis using Seurat’s built-in functions.

## Methods Used:
1. **`FindMarkers()`**: Identifies differentially expressed genes between two groups (e.g., clusters or conditions).
2. **`FindConservedMarkers()`**: Identifies markers that are conserved across conditions or groups.
3. **`FindAllMarkers()`**: Identifies markers for all clusters.

This allows the identification of genes that are significant in distinguishing between groups, and visualizations like feature plots, dot plots, and heatmaps can be used to interpret the results.

### DE Analysis Code

```r
## 3.1 FindMarkers(): DE genes identification
print(DefaultAssay(seuratObj.int))  # Check the current assay being used
DefaultAssay(seuratObj.int) = 'RNA'  # Set the assay to RNA
table(seuratObj.int$expCond2)  # Check distribution of condition labels

deRes1 <- FindMarkers(subset(seuratObj.int, idents = c('4')), 
                      ident.1 = 'A', ident.2 = 'F', group.by = 'expCond2', 
                      test.use = 'MAST', 
                      min.cells.group = 3,  # Define test cells
                      logfc.threshold = 0.25, min.pct = 0.1)  # Define test genes

# Display the DE results for condition 'A' vs 'F'
print(deRes1)

## 4.3 Gene expression check: Feature plot
FeaturePlot(seuratObj.int, features = c("CFD", "IGKC", "SLCO2A1"), 
            split.by = "expCond2", max.cutoff = 3, cols = c("grey","red"), reduction = "umap")

## 3.2 FindConservedMarkers(): Finds markers that are conserved between the groups
deRes2 <- FindConservedMarkers(seuratObj.int, 
                               ident.1 = '4', grouping.var = 'expCond2', 
                               min.cells.group = 3)  # Set minimum cells per group
# Print the top conserved markers
print(head(deRes2))

## 3.3 FindAllMarkers(): Identify highly expressed cluster marker genes
deRes3 <- FindAllMarkers(seuratObj.int, 
                         assay = 'integrated', slot = "data", 
                         only.pos = TRUE,  # Only positive markers
                         test.use = 'MAST', 
                         min.cells.group = 3, 
                         logfc.threshold = 0.25, min.pct = 0.1)  # Define test genes
# Print the top markers for each cluster
print(head(deRes3))

## Dotplot for the top markers
markers.to.plot <- deRes3 %>% 
  group_by(cluster) %>% 
  arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
  slice(1:3) %>% 
  pull(gene)

# Create a DotPlot using the top markers
DefaultAssay(seuratObj.int) = 'integrated'
DotPlot(seuratObj.int, features = unique(markers.to.plot)) + RotatedAxis()

## Heatmap of the top markers
DoHeatmap(seuratObj.int, assay = 'integrated', slot = "scale.data", features = markers.to.plot)
```

## Explanation of the Functions:

+  FindMarkers():
Compares two groups (e.g., conditions or clusters) and identifies differentially expressed genes using statistical tests. Here, we use MAST as the testing method.

+ FindConservedMarkers():
Finds genes that are consistently expressed across multiple conditions or groups. This helps identify conserved markers that are robust across different experimental settings.

+ FindAllMarkers():
Identifies the most significant markers for all clusters in the dataset. This is useful for annotating and understanding the identity of each cluster.

## Visualizing DE Results:

+ FeaturePlot(): Visualizes the expression of specific genes, 
such as CFD, IGKC, and SLCO2A1, on UMAP plots, allowing you to observe their distribution across the dataset.

+ DotPlot(): Displays the expression of the top markers across clusters, 
helping to visualize the relative expression of genes within different cell populations.

+ Heatmap(): Shows the expression patterns of the top markers across cells in the dataset, 
providing a more detailed view of gene expression variations across conditions or clusters.

By identifying DE genes and visualizing their expression, 
we can gain insights into cell-type-specific functions, biological pathways,
and experimental conditions affecting gene expression.