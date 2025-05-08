# Step 5: Doublet Detection Using DoubletDecon

In this step, we perform doublet detection on single-cell RNA-seq data using the **DoubletDecon** package. Two methods—**centroids** and **medoids**—are used to identify potential doublets. The resulting annotations are then integrated into a Seurat object for downstream analysis.

---

## 5.1 Doublet Detection using Centroids

```r
doubletDeconResCentroids <- DoubletDecon::Main_Doublet_Decon(
  rawDataFile = paste0(resProcessDir, '/', filename, "_expression"),
  groupsFile = paste0(resProcessDir, '/', filename , "_groups"),
  filename = filename,
  location = paste(resProcessDir, 'Centroids_', sep = '/'),
  removeCC = FALSE,  # default is FALSE
  species = as.character(doubletDeconSpecies(genomeSpecies)),  # default is 'mmu'
  rhop = doubletDeconRhop,  # x in mean + x*SD to determine upper cutoff for correlation in the blacklist
  PMF = doubletDeconPMF,  # use step 2 (unique gene expression) in doublet determination
  useFull = FALSE,  # use full gene list for PMF analysis
  heatmap = FALSE,
  centroids = TRUE,
  nCores = doubletDeconNoCore
)

print('Doublet centroids detection results table:')
print(table(doubletDeconResCentroids$DRS_doublet_table$isADoublet))
print('END Step 3.1: DoubletDecon centroid detection')
```

---

## 5.2 Doublet Detection using Medoids

```r
print('Start Step 3.2: DoubletDecon medoids detection')
print('=========')

doubletDeconResMedoids <- DoubletDecon::Main_Doublet_Decon(
  rawDataFile = paste0(resProcessDir, '/', filename, "_expression"),
  groupsFile = paste0(resProcessDir, '/', filename , "_groups"),
  filename = filename,
  location = paste(resProcessDir, 'Medoids_', sep = '/'),
  removeCC = FALSE,
  species = as.character(doubletDeconSpecies(genomeSpecies)),
  rhop = doubletDeconRhop,
  PMF = doubletDeconPMF,
  useFull = FALSE,
  heatmap = FALSE,
  centroids = FALSE,
  nCores = doubletDeconNoCore
)
```

---

## 5.3 Integrate Doublet Annotation into Seurat Object

```r
doubletDf <- data.frame(cells = rownames(seuratObjOrg@meta.data))
doubletDf <- doubletDf %>% dplyr::mutate(doublet = ifelse(
  rownames(seuratObjOrg@meta.data) %in% doubletCellsUpdate,
  'TRUE',
  'FASLE'  # Typo kept from original code; consider correcting to 'FALSE'
))

seuratObjOrg <- Seurat::AddMetaData(
  object = seuratObjOrg,
  col.name = 'doublet',
  metadata = as.factor(doubletDf$doublet)
)

seuratObj <- subset(seuratObjOrg, doublet == 'FASLE')
```

---

## Notes

- Ensure all required objects and paths (`resProcessDir`, `filename`, `doubletCellsUpdate`, etc.) are defined in your environment before executing this script.
- Consider fixing the typo `'FASLE'` to `'FALSE'` for consistency and clarity.
- Doublet removal is crucial for reducing noise and improving downstream clustering and differential expression accuracy.
