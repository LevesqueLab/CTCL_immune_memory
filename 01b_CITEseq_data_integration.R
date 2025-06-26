
# This R script has been used to process CITEseq data

library(xlsx)
library(Seurat)
library(Azimuth)
library(pheatmap)
library(cowplot)
library(ggplot2)
library(SeuratData)
library(SeuratDisk)
library(SCpubr)
library(Nebulosa)
library(RColorBrewer)

# Combine all Seurat objects for CITEseq blood samples
scList <- list(pbmc1, pbmc2, pbmc3, pbmc4, pbmc5, pbmc6, pbmc7, pbmc8, pbmc9, pbmc10, pbmc11)

# Select variable features that are common across data sets for integration
features <- SelectIntegrationFeatures(object.list = scList)

# Identify RNA based integration anchors to run canonical correlation analysis
intAnchors <- FindIntegrationAnchors(object.list = scList, anchor.features = features, dims = 1:30,
                                      scale = FALSE, reduction = "cca")

# Integrate all CITEseq data sets based on RNA assay
pbmc.RNA.combined <- IntegrateData(anchorset = intAnchors)

# Run the standard workflow for visualization and clustering
pbmc.RNA.combined <- ScaleData(pbmc.RNA.combined, verbose = FALSE)
pbmc.RNA.combined <- RunPCA(pbmc.RNA.combined, npcs = 20, verbose = FALSE)

ElbowPlot(pbmc.RNA.combined, ndims = 30)

pbmc.RNA.combined <- RunUMAP(pbmc.RNA.combined,
                             dims = 1:20, reduction = "pca",
                             assay = "RNA", reduction.name = "rna.umap",
                             reduction.key = "rnaUMAP_")

pbmc.RNA.combined <- RunTSNE(pbmc.RNA.combined,
                             dims = 1:20,reduction = "pca",
                             assay = "RNA", reduction.name = "rna.tsne",
                             reduction.key = "rnaTSNE_")

pbmc.RNA.combined <- FindNeighbors(pbmc.RNA.combined, reduction = "pca", dims = 1:20)
pbmc.RNA.combined <- FindClusters(pbmc.RNA.combined, resolution = 0.5)

DimPlot(pbmc.RNA.combined, reduction = "rna.umap", split.by = "patient", ncol=4)
DimPlot(pbmc.RNA.combined, reduction = "rna.tsne", split.by = "patient", ncol=4)

DimPlot(pbmc.RNA.combined, reduction = "rna.umap", split.by = "condition")
DimPlot(pbmc.RNA.combined, reduction = "rna.tsne", split.by = "condition")

#######################################################################################################################################

# Change assay for all Seurat obejcts
scList <- lapply(X = scList_ADT, FUN = function(x) {
   DefaultAssay(x) <- "ADT"
})

# Select variable features that are common across data sets for integration
features <- SelectIntegrationFeatures(object.list = scList, assay=c("ADT", "ADT", "ADT", "ADT", "ADT", "ADT","ADT", "ADT", "ADT", "ADT", "ADT"))

# Identify ADT based integration anchors to run canonical correlation analysis
tebeAnchorsADT <- FindIntegrationAnchors(object.list = scList_ADT, anchor.features = features, dims = 1:10, 
                                         assay=c("ADT", "ADT", "ADT", "ADT", "ADT", "ADT", "ADT", "ADT", "ADT", "ADT", "ADT"),
                                         scale = FALSE, reduction = "cca")

# Integrate all CITEseq data sets based on ADT assay
pbmc.ADT.combined <- IntegrateData(anchorset = tebeAnchorsADT)

# Run the standard workflow for visualization and clustering
pbmc.ADT.combined <- ScaleData(pbmc.ADT.combined, verbose = FALSE)
pbmc.ADT.combined <- RunPCA(pbmc.ADT.combined, npcs = 10, verbose = FALSE)
ElbowPlot(pbmc.ADT.combined)

pbmc.ADT.combined <- RunUMAP(pbmc.ADT.combined,
                             dims = 1:10, reduction = "pca",
                             assay = "ADT", reduction.name = "adt.umap",
                             reduction.key = "adtUMAP_")

pbmc.ADT.combined <- RunTSNE(pbmc.ADT.combined,
                             dims = 1:10, reduction = "pca",
                             assay = "ADT", reduction.name = "adt.tsne",
                             reduction.key = "adtTSNE_")

pbmc.ADT.combined <- FindNeighbors(pbmc.ADT.combined, reduction = "pca", dims = 1:10)
pbmc.ADT.combined <- FindClusters(pbmc.ADT.combined, resolution = 0.5)

DimPlot(pbmc.ADT.combined, reduction = "adt.umap", split.by = "patient", ncol=4)
DimPlot(pbmc.ADT.combined, reduction = "adt.tsne", split.by = "patient", ncol=4)

DimPlot(pbmc.ADT.combined, reduction = "adt.umap", split.by = "condition")
DimPlot(pbmc.ADT.combined, reduction = "adt.tsne", split.by = "condition")

#######################################################################################################################################

# Combine the two integrated assays
pbmc.combined <- pbmc.RNA.combined

pbmc.combined[["integrated.adt"]] <- pbmc.ADT.combined[["integrated"]]
pbmc.combined[["apca"]] <- pbmc.ADT.combined[["pca"]]
pbmc.combined[["adt.umap"]] <- pbmc.ADT.combined[["adt.umap"]]
pbmc.combined[["adt.tsne"]] <- pbmc.ADT.combined[["adt.tsne"]]

# Run WNN analysis to integrate different assays
pbmc.combined <- FindMultiModalNeighbors(pbmc.combined, reduction.list = list("pca", "apca"), 
                                         dims.list = list(1:20, 1:10), modality.weight.name = "RNA.weight")

pbmc.combined <- RunUMAP(pbmc.combined, nn.name = "weighted.nn", reduction.name = "umap", reduction.key = "UMAP_")
pbmc.combined <- RunTSNE(pbmc.combined, nn.name = "weighted.nn", reduction.name = "tsne", reduction.key = "TSNE_")
pbmc.combined <- FindClusters(pbmc.combined, graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = FALSE) # resolution = 0.8 for CITEseq skin samples

DimPlot(pbmc.combined, reduction = 'umap', label = TRUE) + NoLegend()
DimPlot(pbmc.combined, reduction = 'tsne', label = TRUE) + NoLegend()

DimPlot(pbmc.combined, reduction = 'umap', label = FALSE, split.by = "patient", ncol=4) + NoLegend()
DimPlot(pbmc.combined, reduction = 'tsne', label = TRUE, split.by = "patient", ncol=4) + NoLegend()

DimPlot(pbmc.combined, reduction = 'umap', group.by = "patient") + NoLegend()
DimPlot(pbmc.combined, reduction = 'tsne', group.by = "patient") + NoLegend()

DimPlot(pbmc.combined, reduction = 'umap', label = TRUE, split.by = "condition", ncol=3) + NoLegend()
DimPlot(pbmc.combined, reduction = 'tsne', label = TRUE, split.by = "condition", ncol=3) + NoLegend()

DimPlot(pbmc.combined, reduction = 'umap', group.by = "condition") + NoLegend()
DimPlot(pbmc.combined, reduction = 'tsne', group.by = "condition") + NoLegend()

# Reference mapping against the Azimuth reference data set
DefaultAssay(scData_sub) <- "RNA"

scData_sub <- RunAzimuth(scData_sub, reference = "pbmcref")

# Predicted cell types level 1
p1 <- DimPlot(scData_sub, reduction = "umap",group.by = "predicted.celltype.l1", label = TRUE,
              label.size = 3, repel = TRUE)  + ggtitle("Query transfered labels") + NoLegend()

p2 <- DimPlot(scData_sub, reduction = "tsne",group.by = "predicted.celltype.l1", label = TRUE,
              label.size = 3, repel = TRUE)  + ggtitle("Query transfered labels") + NoLegend()

p1 + p2

# Predicted cell types level 2
DimPlot(scData_sub, reduction = "umap",group.by = "predicted.celltype.l2", label = TRUE,
        label.size = 3, repel = TRUE)  + ggtitle("Query transfered labels") + NoLegend()

DimPlot(scData_sub, reduction = "umap",group.by = "predicted.celltype.l3", label = TRUE,
        label.size = 3, repel = TRUE)  + ggtitle("Query transfered labels") + NoLegend()

p1 <- FeaturePlot(scData_sub,features = "mapping.score", reduction = "umap")
p2 <- FeaturePlot(scData_sub,features = "predicted.celltype.l2.score", reduction = "umap")

p1 | p2

# Identify cluster markers
markers <- FindAllMarkers(scData_sub, 
                          only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, return.thresh = 0.01)

write_tsv(as_tibble(markers), file = "/srv/GT/analysis/adhideb/p28409-Andrea/SS_HD_PBMCs/allCells/pos_markers_perClust.tsv")

# Visualize top 10 cluster markers
top10 <- markers %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC)

DotPlot(scData_sub, features = unique(top10$gene), cols = c("lightgrey","red"), col.min = 0, col.max = 1) +
  theme(axis.text.x = element_text(angle = 90))

