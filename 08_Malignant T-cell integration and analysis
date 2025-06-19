#!/usr/bin/env Rscript

#' Malignant T-cell Integration (Blood and Skin)
#'
#' This script integrates malignant T-cell populations from CTCL blood and skin samples
#' using multiple strategies (unintegrated reclustering, CCA, and BBKNN). It supports
#' visualization of tissue origin, patient identity, and malignant subtype annotations.
#'
#' Dependencies: Seurat, SCpubr, bbknnR, dplyr, future, RColorBrewer, qs

# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(bbknnR)
  library(SCpubr)
  library(future)
  library(RColorBrewer)
  library(qs)
})

# Set parallelization and random seed
plan("multicore", workers = 4)
set.seed(123)

#' Load and preprocess Seurat objects
load_malignant_samples <- function(blood_path, skin_path) {
  cat("Loading malignant T-cell Seurat objects...\n")
  
  blood <- readRDS(blood_path)
  skin <- readRDS(skin_path)
  
  # Set identity to malignant subtype column
  blood <- SetIdent(blood, value = blood$TCM)
  skin <- SetIdent(skin, value = skin$TCM)
  
  DefaultAssay(blood) <- "RNA"
  DefaultAssay(skin) <- "RNA"
  
  return(list(blood = blood, skin = skin))
}

#' Harmonize metadata across samples
harmonize_metadata <- function(blood, skin) {
  cat("Aligning metadata and tagging tissue origin...\n")
  
  meta_b <- blood@meta.data %>% mutate(condition = ifelse(status == "CTCL", "SS", as.character(status)))
  meta_s <- skin@meta.data
  
  common_cols <- intersect(names(meta_b), names(meta_s))
  meta_b_common <- meta_b[, common_cols]
  meta_s_common <- meta_s[, common_cols]
  
  meta_b_common$tissue <- "blood"
  meta_s_common$tissue <- "skin"
  
  blood <- AddMetaData(blood, meta_b_common)
  skin <- AddMetaData(skin, meta_s_common)
  
  return(list(blood = blood, skin = skin))
}

#' Merge and normalize data
normalize_and_merge <- function(blood, skin) {
  cat("Merging and normalizing data with SCTransform...\n")
  
  blood <- SCTransform(blood, verbose = FALSE)
  skin <- SCTransform(skin, verbose = FALSE)
  
  merged <- merge(skin, blood)
  return(merged)
}

#' Run unintegrated reclustering by tissue
recluster_unintegrated <- function(seu_obj) {
  cat("Performing unintegrated clustering split by tissue...\n")
  
  DefaultAssay(seu_obj) <- "RNA"
  seu_obj[["RNA"]] <- JoinLayers(seu_obj[["RNA"]])
  seu_obj[["RNA"]] <- split(seu_obj[["RNA"]], f = seu_obj$tissue)
  
  seu_obj <- SCTransform(seu_obj)
  seu_obj <- RunPCA(seu_obj)
  seu_obj <- FindNeighbors(seu_obj, dims = 1:30)
  seu_obj <- FindClusters(seu_obj, resolution = 1, cluster.name = "unintegrated_clusters")
  seu_obj <- RunUMAP(seu_obj, dims = 1:30, reduction.name = "umap.unintegrated")
  
  return(seu_obj)
}

#' Run CCA-based integration
run_cca_integration <- function(seu_obj) {
  cat("Running CCA-based integration...\n")
  
  seu_obj <- IntegrateLayers(seu_obj, method = CCAIntegration, normalization.method = "SCT")
  seu_obj[["RNA"]] <- JoinLayers(seu_obj[["RNA"]])
  
  seu_obj <- FindNeighbors(seu_obj, reduction = "integrated.dr", dims = 1:30)
  seu_obj <- FindClusters(seu_obj, resolution = seq(0.1, 1, 0.1))
  seu_obj <- RunUMAP(seu_obj, reduction = "integrated.dr", dims = 1:30)
  
  return(seu_obj)
}

#' Run BBKNN integration
run_bbknn_integration <- function(seu_obj) {
  cat("Running BBKNN integration...\n")
  
  seu_obj$patient <- paste(seu_obj$tissue, seu_obj$patient, sep = "_")
  seu_obj <- subset(seu_obj, subset = patient != "CTCL_10")  # remove problematic sample
  
  seu_obj <- RunBBKNN(seu_obj, batch_key = "patient")
  seu_obj <- FindClusters(seu_obj, graph.name = "bbknn", resolution = seq(0.1, 1, 0.1))
  
  return(seu_obj)
}

#' Visualize UMAPs
plot_umaps <- function(seu_obj, label_group = "tissue", split_group = "Malignant_T") {
  p1 <- DimPlot(seu_obj, reduction = "umap", group.by = label_group, label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
  p2 <- DimPlot(seu_obj, reduction = "umap", group.by = split_group)
  print(p1 + p2)
}

# Main pipeline
main <- function() {
  # Load and merge data
  paths <- list(
    blood = "/BLOOD_Tcells.rds",
    skin = "/SKIN_Malignant_T.rds"
  )
  
  data <- load_malignant_samples(paths$blood, paths$skin)
  data <- harmonize_metadata(data$blood, data$skin)
  merged <- normalize_and_merge(data$blood, data$skin)
  
  # --- Unintegrated Reclustering ---
  reclustered <- recluster_unintegrated(merged)
  DimPlot(reclustered, reduction = "umap.unintegrated", group.by = c("tissue", "patient"))
  
  # --- CCA Integration ---
  integrated <- run_cca_integration(merged)
  DimPlot(integrated, reduction = "integrated.dr", group.by = c("tissue", "patient"))
  
  # --- BBKNN Integration ---
  bbknn_result <- run_bbknn_integration(integrated)
  plot_umaps(bbknn_result)
  
  # Optional save
  # saveRDS(bbknn_result, "/CTCL_Malignant_T.rds")
}

# Run pipeline
main()

