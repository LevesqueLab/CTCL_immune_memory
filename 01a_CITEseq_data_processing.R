
# This R script has been used to process CITEseq data

library(ezRun)
library(cowplot)
library(kableExtra)
library(fastICA)
library(Matrix)
library(RColorBrewer)
library(scDblFinder)
library(SingleCellExperiment)
library(scater)
library(Seurat)
library(SingleR)
library(plotly)
library(tidyverse)
library(dsb)
library(magrittr)
library(Azimuth)

# Following code chunk has been illustrated for single CITEseq sample
# Same code chunk has been used for remaining all CITEseq samples

# Read the antibody derived tags
fileADT <- "../feature_ref.csv"
ADTs <- read_csv(fileADT)$id

# Load raw and filtered count matrices 
raw <- Seurat::Read10X("../raw_feature_bc_matrix")
cells <- Seurat::Read10X("../filtered_feature_bc_matrix")

# Define cell-containing barcodes and separate cells and empty drops
stained_cells <- colnames(cells$`Gene Expression`)
background <- setdiff(colnames(raw$`Gene Expression`), stained_cells)

# Split the raw data into RNA and ADT matrices
prot <- raw$`Antibody Capture`
rna <- raw$`Gene Expression`

# Update ADT names
rownames(prot) <- ADTs

# Remove self conjugating antibodies
prot <- prot[-which(rownames(prot) %in% c("Aquaporin1", "GLUT1", "IGF1R")),]

# Create metadata of droplet QC stats used in standard scRNAseq processing
mtgene <- grep(pattern = "^MT-", rownames(rna), value = TRUE) # used below
ribogene <- grep(pattern = "^RP[SL]", rownames(rna), value = TRUE) # used below

md <- data.frame(
  rna.size = log10(Matrix::colSums(rna)), 
  prot.size = log10(Matrix::colSums(prot)), 
  n.gene = Matrix::colSums(rna > 0), 
  mt.prop = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna),
  ribo.prop = Matrix::colSums(rna[ribogene, ]) / Matrix::colSums(rna)
)

# Add indicator for barcodes that Cell Ranger called as cells
md$drop.class <- ifelse(rownames(md) %in% stained_cells, 'cell', 'background')

# Remove barcodes with no evidence of capture in the experiment
md <- md[md$rna.size > 0 & md$prot.size > 0, ] 
qplot(prot.size, rna.size, data = md, facets = . ~ drop.class, colour = n.gene) 
qplot(prot.size, rna.size, data = md, facets = . ~ drop.class, colour = mt.prop) 
qplot(prot.size, rna.size, data = md, facets = . ~ drop.class, colour = ribo.prop) 

background_drops = rownames(
  md[ md$prot.size > 1 & 
        md$prot.size < 3 & 
        md$rna.size < 3, ]
) 
background.adt.mtx <- as.matrix(prot[ , background_drops]) 

# Calculate statistical thresholds for droplet filtering.
cellmd <- md[md$drop.class == 'cell', ]

# Filter drops with + / - 3 median absolute deviations from the median library size
rna.mult <- (3*mad(cellmd$rna.size))
prot.mult <- (3*mad(cellmd$prot.size))
rna.lower <- median(cellmd$rna.size) - rna.mult
rna.upper <- median(cellmd$rna.size) + rna.mult
prot.lower <- median(cellmd$prot.size) - prot.mult
prot.upper <- median(cellmd$prot.size) + prot.mult

# Filter rows based on droplet qualty control metrics
qc_cells = rownames(
  cellmd[cellmd$prot.size > prot.lower & 
           cellmd$prot.size < prot.upper & 
           cellmd$rna.size > rna.lower & 
           cellmd$rna.size < rna.upper, ]
)

length(qc_cells) 

cell.adt.raw <- as.matrix(prot[ , qc_cells])
cell.rna.raw <- rna[ ,qc_cells]
cellmd <- cellmd[qc_cells, ]

# Define isotype controls (specific to antibodies)
isotype.controls = c("ahIgG", "mIgG1", "mIgG2a", "mIgG2b", "rIgG1", "rIgG2a", "rIgG2b")

# Initialize Seurat object with raw filtered data
scData <- CreateSeuratObject(counts = cell.rna.raw, project = "PBMC")

# Add ADT assay
scData[["ADT"]] <- CreateAssayObject(cell.adt.raw)

pvalue_allMarkers <- 0.01
pvalue_all2allMarkers <- 0.01

# Calculate Mito and Ribo gene %
scData[["percent.mt"]] <- PercentageFeatureSet(scData,  pattern = "^MT-")
scData[["percent.ribo"]] <- PercentageFeatureSet(scData,  pattern = "^RP[SL]")

# The following cut-offs are for CITEseq blood samples
cellsToKeep <- scData$nFeature_RNA < 4000 &
  scData$nCount_RNA < 15000 & 
  scData$percent.mt < 7.5 

# # The following cut-offs are for CITEseq skin samples
# cellsToKeep <- scData$nFeature_RNA < 6000 &
#   scData$nCount_RNA < 40000 & 
#   scData$percent.mt < 20

scData_sub <- scData[, cellsToKeep]

#### RNA assay: Pre-processing & dimensional reduction
  
# Normalize data 
scData_sub <- NormalizeData(scData_sub, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify variable features
scData_sub <- FindVariableFeatures(scData_sub)

# Scale data
scData_sub <- ScaleData(scData_sub)

# Perform linear dimensional reduction
scData_sub <- RunPCA(scData_sub)

# Determine the 'dimensionality' 
  
p <- ElbowPlot(scData_sub, ndims = 50)
p

# Cluster cells based on RNA
scData_sub <- FindNeighbors(scData_sub, dims = 1:30)
scData_sub <- FindClusters(scData_sub, resolution = 0.5)

scData_sub <- RunUMAP(scData_sub,
                      dims = 1:30, reduction = "pca",
                      assay = "RNA", reduction.name = "rna.umap",
                      reduction.key = "rnaUMAP_",
                      n.neighbors = getPerplexity(ncol(scData_sub))
)
scData_sub <- RunTSNE(scData_sub,
                      dims = 1:30,reduction = "pca",
                      assay = "RNA", reduction.name = "rna.tsne",
                      reduction.key = "rnaTSNE_",
                      perplexity = getPerplexity(ncol(scData_sub))
)


#### ADT assay: Pre-processing & dimensional reduction
  
# Change assay
DefaultAssay(scData_sub) <- "ADT"

# Update the raw ADT matrix
qc_cells2 <- names(cellsToKeep)[which(cellsToKeep==TRUE)]
cell.adt.raw <- as.matrix(prot[ , qc_cells2])

# Normalize and denoise with dsb  
cells.dsb.norm <- DSBNormalizeProtein(
  cell_protein_matrix = cell.adt.raw, 
  empty_drop_matrix = background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls,
  scale.factor = 'mean.subtract')

# Update Seurat object with norm ADT counts
scData_sub[["ADT"]] = Seurat::CreateAssayObject(data = cells.dsb.norm)

# Identify variable features
VariableFeatures(scData_sub) <- rownames(scData_sub[["ADT"]])[-which(rownames(scData_sub[["ADT"]]) %in% isotype.controls)]

# Scale the data
scData_sub <- ScaleData(scData_sub)

# Perform linear dimensional reduction
scData_sub <- RunPCA(scData_sub, reduction.name = 'apca')

# Determine the 'dimensionality' 
  
p <- ElbowPlot(scData_sub, ndims = 50, reduction = "apca")
p

# Cluster cells based on ADT
scData_sub <- FindNeighbors(scData_sub, dims = 1:18, reduction = "apca")
scData_sub <- FindClusters(scData_sub, resolution = 0.5)

scData_sub <- RunUMAP(scData_sub,
                      dims = 1:18, reduction = "apca",
                      assay = "ADT", reduction.name = "adt.umap",
                      reduction.key = "adtUMAP_",
                      n.neighbors = getPerplexity(ncol(scData_sub))
)

scData_sub <- RunTSNE(scData_sub,
                      dims = 1:18, reduction = "apca",
                      assay = "ADT", reduction.name = "adt.tsne",
                      reduction.key = "adtTSNE_",
                      perplexity = getPerplexity(ncol(scData_sub))
)


#### WNN: Cluster the cells based on multimodal neighbors
  
# WNN analysis based on RNA + ADT assays
scData_sub <- FindMultiModalNeighbors(scData_sub, reduction.list = list("pca", "apca"), 
                                      dims.list = list(1:30, 1:18))

scData_sub <- FindClusters(scData_sub, graph.name = "wsnn", algorithm = 3, resolution = 0.6, verbose = FALSE) 

scData_sub <- RunUMAP(scData_sub, nn.name = "weighted.nn", 
                      reduction.name = "umap", reduction.key = "UMAP_")

scData_sub <- RunTSNE(scData_sub, nn.name = "weighted.nn", 
                      reduction.name = "tsne", reduction.key = "TSNE_")

# RNA clustering colored by WNN clusters 
p1 <- DimPlot(scData_sub,  reduction = "rna.umap", label = TRUE)
p2 <- DimPlot(scData_sub,  reduction = "rna.tsne", label = TRUE)
p <- plot_grid(p1, p2, labels = "AUTO", ncol = 2)
p

# ADT clustering colored by WNN clusters 
p1 <- DimPlot(scData_sub,  reduction = "adt.umap", label = TRUE)
p2 <- DimPlot(scData_sub,  reduction = "adt.tsne", label = TRUE)
p <- plot_grid(p1, p2, labels = "AUTO", ncol = 2)
p

# Modality weight distribution
v1 <- VlnPlot(scData_sub, features = "RNA.weight", pt.size = 0.1)
v2 <- VlnPlot(scData_sub, features = "ADT.weight", pt.size = 0.1)

v <- plot_grid(v1, v2, labels = "AUTO", ncol = 2)
v

# Reference mapping against the Azimuth reference data set
scData_sub <- RunAzimuth(scData_sub, reference = "pbmcref")

# Query mapping on the reference UMAP
DimPlot(scData_sub, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE,
        label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")

# Predicted cell types
p1 <- DimPlot(scData_sub, reduction = "umap",group.by = "predicted.celltype.l2", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend()

p2 <- DimPlot(scData_sub, reduction = "tsne",group.by = "predicted.celltype.l2", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend()

p1 + p2


# Query mapping confidence
FeaturePlot(scData_sub,features = c("mapping.score", "predicted.celltype.l2.score"), reduction = "umap")

# Marker genes
DefaultAssay(scData_sub) <- "RNA"

markers <- FindAllMarkers(scData_sub, 
                          only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, return.thresh = pvalue_allMarkers
)

posMarkersFn <- "pos_markers.tsv"

write_tsv(as_tibble(markers), file = posMarkersFn)

# Heatmap of top 10 cluster marker genes
top10 <- markers %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC)
DoHeatmap(scData_sub,  features = top10$gene) + NoLegend()

# Feature plots to visualize multiple modalities side-by-side
featureADT <- read.csv("../CITEseq_ADTs.csv")

for(i in 1:nrow(featureADT)) {
  DefaultAssay(scData_sub) <- "ADT"
  
  p1 <- FeaturePlot(scData_sub, featureADT$id[i], cols = c("lightgrey", "darkgreen"), reduction = "umap") + 
    ggtitle(paste0(featureADT$id[i], " protein"))
  
  DefaultAssay(scData_sub) <- "RNA"
  
  p2 <- DimPlot(scData_sub, cols = rep("gray", length(levels(Idents(scData_sub)))), reduction = "umap") + theme(legend.position = "none") + 
    ggtitle(paste0(featureADT$id[i], " RNA"))
  tryCatch({
    p2 = FeaturePlot(scData_sub, featureADT$Gene[i], reduction = "umap") + ggtitle(paste0(featureADT$Gene[i], " RNA"))},
    error=function(e){
      return()
    }
  )
  print(p1 | p2)
}

# Ridge plots to visualize ADTs
DefaultAssay(scData_sub) <- "ADT"

for (p in featureADT$id) 
{ 
  rp <- RidgePlot(scData_sub, features= p, assay = "ADT") + theme(legend.position = "none")
  print(rp)
}
