#!/usr/bin/env Rscript

#' CTCL Malignancy Gene Signature Scoring and ROC Evaluation
#'
#' This script scores T-cells from CTCL blood and skin samples based on known 
#' malignant gene signatures from published datasets. It computes UCell scores 
#' for consensus gene sets and evaluates classification performance using ROC curves.
#'
#' Based on the CTCL single-cell analysis pipeline.

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(pROC)
  library(ggplot2)
  library(dplyr)
  library(readxl)
  library(ggVennDiagram)
  library(viridis)
  library(patchwork)
})

# Set analysis parameters
params <- list(
  skin_path = "/SKIN_Tcells.rds",
  blood_path = "/BLOOD_Tcells.rds",
  cores = 20,
  target_label = "Malignant TCM"
)

#' Load and preprocess skin Seurat object
load_skin_data <- function(path) {
  cat("Loading skin T-cell data...\n")
  seu <- readRDS(path)
  meta <- seu@meta.data
  meta$Tcells[meta$Tcells %in% c("Malignant TCM", "Malignant IS")] <- "Malignant TCM"
  seu@meta.data <- meta
  seu <- SetIdent(seu, value = meta$Tcells)
  return(seu)
}

#' Load blood Seurat object
load_blood_data <- function(path) {
  cat("Loading blood T-cell data...\n")
  seu <- readRDS(path)
  return(seu)
}

#' Load and filter malignant gene signatures from published datasets
load_signatures <- function() {
  cat("Loading malignant gene signatures...\n")
  
  Rindler <- read.csv("/rindler_DEGs.csv") %>%
    filter(p_val_adj < 0.001 & avg_log2FC > 0) %>%
    pull(X)
  
  Jonak <- read_excel("/Jonak_Signature.xlsx") %>%
    filter(P_val_adj < 0.001 & avg_logFC > 0) %>%
    pull(Gene)
  
  Li <- read_excel("/Li_signature.xlsx", sheet = "Supplementary Table 2")
  colnames(Li) <- Li[1,]
  Li <- Li[-1,] %>%
    mutate(across(c(padj, log2FoldChange), as.numeric)) %>%
    filter(padj < 0.001 & log2FoldChange > 1) %>%
    pull(Gene)
  
  Borcherding <- read_excel("/Borcherding_signature.xlsx") %>%
    filter(p_val_adj < 0.001 & avg_logFC > 1) %>%
    pull(genes)
  
  mal.list <- list(
    'Jonak (s)' = Jonak,
    'Rindler (s)' = Rindler,
    'Borcherding (b)' = Borcherding,
    'Li (s)' = Li
  )
  
  return(mal.list)
}

#' Identify consensus gene sets shared across multiple studies
get_shared_signatures <- function(gene_lists) {
  cat("Finding overlapping genes across signatures...\n")
  item_counts <- table(unlist(gene_lists))
  list(
    '2 genesets' = names(item_counts[item_counts >= 2]),
    '3 genesets' = names(item_counts[item_counts >= 3]),
    '4 genesets' = names(item_counts[item_counts >= 4])
  )
}

#' Score UCell signatures and compute ROC curves
score_and_plot_signatures <- function(seu_obj, sample_label, gene_sets, positive_class, ncores = 4) {
  cat("Scoring UCell and computing ROC for:", sample_label, "\n")
  
  seu_score <- AddModuleScore_UCell(seu_obj, features = gene_sets, assay = "RNA", ncores = ncores)
  metadata <- seu_score@meta.data
  metadata$Target <- ifelse(metadata$Tcells == positive_class, 1, 0)
  
  sig_cols <- grep("_UCell$", colnames(metadata), value = TRUE)
  roc_list <- lapply(sig_cols, function(sig) roc(metadata$Target, metadata[[sig]], quiet = TRUE))
  
  roc_df <- do.call(rbind, lapply(seq_along(roc_list), function(i) {
    data.frame(
      Sensitivity = roc_list[[i]]$sensitivities,
      Specificity = 1 - roc_list[[i]]$specificities,
      Signature = sig_cols[i]
    )
  }))
  
  roc_df$Signature <- recode(roc_df$Signature,
                              "2 genesets_UCell" = "shared 2 genesets",
                              "3 genesets_UCell" = "shared 3 genesets",
                              "4 genesets_UCell" = "shared 4 genesets")
  
  # Define colors
  signature_colors <- c(
    "shared 2 genesets" = "#FF033E",
    "shared 3 genesets" = "coral",
    "shared 4 genesets" = "#DCAE96"
  )
  
  roc_plot <- ggplot(roc_df, aes(x = Specificity, y = Sensitivity, color = Signature)) +
    geom_line(size = 1) +
    scale_color_manual(values = signature_colors) +
    theme_minimal() +
    labs(title = sample_label, x = "1 - Specificity", y = "Sensitivity", color = "Signature") +
    theme(legend.position = "right")
  
  feat_plot <- FeaturePlot(seu_score, features = "2 genesets_UCell", ncol = 1, pt.size = 0.8, 
                           order = TRUE, reduction = "umap", max.cutoff = 'q99') + 
    scale_color_viridis(option = "plasma")
  
  print(roc_plot | feat_plot)
}

#' Plot Venn diagram of gene signature overlap
plot_signature_venn <- function(gene_lists) {
  cat("Plotting Venn diagram...\n")
  ggVennDiagram(gene_lists, label = "count", set_size = 6, label_size = 7) +
    scale_x_continuous(expand = expansion(mult = .2)) +
    scale_fill_gradientn(
      colors = c("darkred", "red", "pink", "white"),
      values = c(0, 0.3, 0.6, 1),
      limits = c(0, 450)
    )
}

#' Main pipeline function
main <- function() {
  # Load data
  skin <- load_skin_data(params$skin_path)
  blood <- load_blood_data(params$blood_path)
  
  # Load signatures
  signatures <- load_signatures()
  comb_signatures <- get_shared_signatures(signatures)
  
  # Score and plot
  score_and_plot_signatures(skin, "Skin", comb_signatures, params$target_label, params$cores)
  score_and_plot_signatures(blood, "Blood", comb_signatures, params$target_label, params$cores)
  
  # Venn diagram
  plot_signature_venn(signatures)
}

# Run the pipeline
main()



