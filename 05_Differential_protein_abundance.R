#!/usr/bin/env Rscript

#' CTCL Blood: Differential Protein Abundance Analysis
#'
#' This script performs differential abundance testing of surface proteins (ADT)
#' between malignant T-cells and non-malignant CD4 T-cells in CTCL blood samples.
#' It also explores marker proteins within malignant subpopulations using LR tests
#' with patient-level covariate correction.
#'
#' Dependencies: Seurat, dplyr, ggplot2, qs (for qread)

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(future)
  library(qs)
})

# Set up parallel processing
plan("multicore", workers = 4)

# ------------------------------------------------------------------------------
# Section 1: Between-group ADT comparison (Malignant TCM vs. CD4 T-cells)
# ------------------------------------------------------------------------------

#' Load and preprocess Seurat object for blood T-cell protein comparison
load_blood_tcells <- function(path) {
  cat("Loading and subsetting blood T-cell data...\n")
  seu <- readRDS(path)
  seu_sub <- subset(seu, subset = Tcells %in% c("Malignant TCM", "CD4 naive T", "CD4 TCM", "CD4 TEM"))
  
  # Merge CD4 subtypes
  meta <- seu_sub@meta.data
  meta$Tcells[meta$Tcells %in% c("CD4 naive T", "CD4 TCM", "CD4 TEM")] <- "CD4 T"
  seu_sub <- AddMetaData(seu_sub, meta)
  return(seu_sub)
}

#' Perform differential ADT analysis (Malignant TCM vs CD4 T)
run_bulk_marker_test <- function(seu_obj) {
  cat("Running differential protein abundance (Malignant TCM vs CD4 T)...\n")
  DefaultAssay(seu_obj) <- "ADT"
  seu_obj <- SetIdent(seu_obj, value = seu_obj$Tcells)
  markers <- FindAllMarkers(
    object = seu_obj,
    test.use = "LR",
    latent.vars = "patient",
    return.thresh = 0.05
  )
  markers$Protein <- rownames(markers)
  return(markers)
}

#' Plot top protein markers between Malignant TCM vs CD4 T
plot_top_bulk_markers <- function(markers_df, seu_obj, top_n = 8) {
  top_proteins <- markers_df %>%
    filter(p_val_adj < 0.01) %>%
    top_n(top_n, avg_log2FC) %>%
    pull(gene)

  p <- DotPlot(seu_obj, features = top_proteins, group.by = "Tcells") +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_colour_gradient2(low = "lightgreen", mid = "lightgrey", high = "coral") +
    labs(title = "Top Differential Proteins (Malignant vs CD4 T)")
  
  print(p)
}

# ------------------------------------------------------------------------------
# Section 2: Within-malignant subgroup comparison
# ------------------------------------------------------------------------------

#' Load malignant T-cell Seurat object
load_malignant_cells <- function(path) {
  cat("Loading malignant T-cell object...\n")
  seu <- qread(path)
  DefaultAssay(seu) <- "ADT"
  seu <- SetIdent(seu, value = seu$Malignant_T)
  return(seu)
}

#' Run differential protein analysis within malignant subtypes
run_malignant_marker_test <- function(seu_obj) {
  cat("Running intra-malignant subtype marker discovery...\n")
  markers <- FindAllMarkers(
    object = seu_obj,
    test.use = "LR",
    latent.vars = "patient",
    return.thresh = 0.05
  )
  markers$Protein <- rownames(markers)
  return(markers)
}

#' Plot top proteins from a selected malignant subtype
plot_malignant_markers <- function(markers_df, seu_obj, subtype = "mTCM", top_n = 5) {
  top_proteins <- markers_df %>%
    filter(cluster == subtype) %>%
    group_by(cluster) %>%
    top_n(top_n, avg_log2FC) %>%
    pull(Protein) %>%
    gsub("\\.1$", "", .) %>%
    unique()

  p <- DotPlot(seu_obj, features = top_proteins, group.by = "Malignant_T") +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_colour_gradient2(low = "lightgreen", mid = "lightgrey", high = "coral") +
    labs(title = paste("Top ADT Markers in", subtype))
  
  print(p)
}

# ------------------------------------------------------------------------------
# Run Analysis
# ------------------------------------------------------------------------------

main <- function() {
  # Part 1: Malignant vs CD4 T
  seu_blood <- load_blood_tcells("/BLOOD_Tcells.rds")
  markers_bulk <- run_bulk_marker_test(seu_blood)
  plot_top_bulk_markers(markers_bulk, seu_blood)

  # Part 2: Within malignant subtypes
  seu_malignant <- load_malignant_cells("/CTCLblood_Malignant_T.qs")
  markers_malignant <- run_malignant_marker_test(seu_malignant)
  plot_malignant_markers(markers_malignant, seu_malignant)
}

# Execute
main()


