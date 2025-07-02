#!/usr/bin/env Rscript

#' CTCL Malignant T-cell PCA and Cluster Similarity Analysis
#'
#' This script analyzes heterogeneity and similarity in malignant T-cells
#' from CTCL blood and skin samples. It includes PCA on pseudobulked samples
#' and pairwise MSE distance computations using all genes.
#'
#' Dependencies: Seurat, edgeR, ggfortify, pheatmap, qs

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggfortify)
  library(edgeR)
  library(pheatmap)
  library(qs)
})

#' Load malignant T-cell Seurat object and preprocess
load_malignant_tcells <- function(path) {
  seu <- qread(path)
   seu <- SetIdent(seu, value = seu$TCM)
  DefaultAssay(seu) <- "RNA"
  seu$tpc <- paste0(seu$TCM, " ", seu$patient)
  seu$tissue_cells <- paste0(seu$tissue, " ", seu$TCM)
  seu <- SetIdent(seu, value = seu$tissue_cells)
  return(seu)
}

#' Perform PCA on pseudobulked gene expression
run_pseudobulk_pca <- function(seu) {
  pb_expr <- AggregateExpression(seu, group.by = "tpc", assays = "RNA")$RNA
  pb_expr <- pb_expr[rowSums(pb_expr) > 0, ]
  samples <- colnames(pb_expr)

  metadata <- data.frame(
    Sample = samples,
    TCM = factor(sub(" .*", "", samples), levels = c("mTCM", "mTRM", "mCyto")),
    Patient = sub(".* ", "", samples),
    row.names = samples
  )

  dge <- edgeR::DGEList(counts = pb_expr, group = metadata$TCM)
  dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  normCounts <- cpm(dge, log = TRUE)

  pca <- prcomp(t(normCounts), scale. = TRUE)

  plot <- autoplot(pca, data = metadata, colour = "TCM", shape = "Patient", frame = TRUE, size = 4) +
    scale_colour_manual(values = c("mTCM" = "darkgreen", "mTRM" = "darkred", "mCyto" = "darkblue")) +
    theme_classic()

  return(plot)
}

#' Compute MSE distance between cluster centroids
compute_mse_distance <- function(expr, groups) {
  group_ids <- unique(groups)
  mse <- function(x, y) mean((x - y)^2)

  centroids <- sapply(group_ids, function(id) {
    rowMeans(expr[, groups == id, drop = FALSE])
  })

  mat <- outer(group_ids, group_ids, Vectorize(function(a, b) {
    mse(centroids[, a], centroids[, b])
  }))

  rownames(mat) <- colnames(mat) <- group_ids
  return(mat)
}

#' Plot distance heatmap
plot_distance_heatmap <- function(dist_matrix, title, color_palette = viridis::viridis) {
  pheatmap(dist_matrix,
           color = color_palette(100),
           display_numbers = TRUE,
           number_color = "black",
           main = title)
}

# Main analysis workflow
main <- function() {
  # Load Seurat object
  
  # PCA on pseudobulk
  cat("Running pseudobulk PCA...\n")
  p <- run_pseudobulk_pca(seu)
  print(p)

  # MSE distance using all genes
  cat("Computing MSE distances (all genes)...\n")
  expr <- GetAssayData(seu, slot = "data")
  dist_all <- compute_mse_distance(expr, seu$tissue_cells)
  plot_distance_heatmap(dist_all, "MSE Distance Between Clusters - All Genes", viridis::plasma)
}

# Execute
main()

