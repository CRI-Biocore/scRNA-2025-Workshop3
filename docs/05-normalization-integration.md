# Seurat Object Normalization and Integration Execution

## Step 1: Normalize and preprocess each individual sample

```r
seuratObj = seurat1.filter
seuratObj <- NormalizeData(seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
seuratObj <- FindVariableFeatures(seuratObj, selection.method = 'vst', nfeatures = 1000)
seuratObj <- ScaleData(object = seuratObj)
seuratObj <- RunPCA(object = seuratObj)
seuratObjList[[i]] = seuratObj
names(seuratObjList) = names(inputPath)
```

> Repeat this for each sample and store in `seuratObjList`.

---

## Step 2: Select integration features

```r
features <- SelectIntegrationFeatures(object.list = seuratObjList, nfeatures = 500)
anchors  <- FindIntegrationAnchors(object.list = seuratObjList,
                                   anchor.features = features,
                                   reduction = 'rpca')
```

### Alternative:

```r
anchors <- FindIntegrationAnchors(object.list = seuratObjList,
                                  anchor.features = 500,
                                  reduction = 'rpca')
```

---

## Step 3: Integrate datasets using the precomputed anchor set

```r
seuratObjIntegrated <- IntegrateData(anchorset = anchors, k.weight = 30)
seuratObjIntegrated$expCond3 = paste(seuratObjIntegrated$expCond1,
                                     seuratObjIntegrated$expCond2,
                                     sep = '_')
print(head(seuratObjIntegrated))
print(table(seuratObjIntegrated$expCond3))
print(table(seuratObjIntegrated$expCond2))
print(table(seuratObjIntegrated$expCond1))
```

> Use the integrated object for clustering, visualization, and further analysis.
