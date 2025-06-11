#!/usr/bin/env Rscript

#' T-cell Integration and Clustering Analysis (without TCR genes)
#' 
#' This script performs T-cell clustering and annotation without TCR genes
#' to avoid bias from clonal expansion. It includes both blood and skin
#' T-cell analysis from the CTCL single-cell memory manuscript.
#'
#' Based on the CTCL manuscript analysis pipeline

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SCpubr)
  library(ggplot2)
  library(viridis)
  library(paletteer)
  library(RColorBrewer)
  library(dplyr)
  library(UCell)
  library(SignatuR)
  library(patchwork)
  library(future)
})

# Set up parallel processing
param <- list(cores = 4, ram = 60, npcs = 20, seed = 38)
plan("multicore", workers = param$cores)
set.seed(param$seed)
options(future.globals.maxSize = param$ram * 1024^3)

#' Remove TCR genes from variable features and re-cluster
#'
#' @param seurat_obj Seurat object to process
#' @param assay Assay to use (default: "SCT")
#' @param reduction_name Name for new PCA reduction
#' @param umap_name Name for new UMAP reduction
#' @return Seurat object with TCR genes removed and re-clustered
remove_tcr_genes_and_recluster <- function(seurat_obj, assay = "SCT", 
                                          reduction_name = "pca.notcr",
                                          umap_name = "umap.notcr") {
  
  cat("Removing TCR genes and re-clustering...\n")
  
  # Get TCR genes from SignatuR
  tcr_genes <- GetSignature(SignatuR$Hs$Compartments$TCR)
  cat("Found", length(tcr_genes), "TCR genes to remove\n")
  
  # Set default assay
  DefaultAssay(seurat_obj) <- assay
  
  # Remove TCR genes from variable features
  v_features <- VariableFeatures(seurat_obj)
  v_features_no_tcr <- v_features[!v_features %in% tcr_genes]
  VariableFeatures(seurat_obj) <- v_features_no_tcr
  
  cat("Removed", length(v_features) - length(v_features_no_tcr), "TCR genes from variable features\n")
  
  # Re-run PCA, clustering, and UMAP
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE, npcs = param$npcs, 
                      reduction.name = reduction_name)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction_name, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = reduction_name, 
                       reduction.name = umap_name, dims = 1:param$npcs, verbose = FALSE)
  
  return(seurat_obj)
}

#' Annotate T-cell subtypes in blood samples
#'
#' @param seurat_obj Seurat object containing blood T-cells
#' @return Seurat object with T-cell annotations
annotate_blood_tcells <- function(seurat_obj) {
  
  cat("Annotating blood T-cell subtypes...\n")
  
  # Prepare metadata
  meta <- seurat_obj@meta.data
  
  # Clean up patient names
  meta$patient <- gsub("^P(\\d+)$", "blood \\1", meta$patient)
  meta$patient <- gsub("^HD(\\d+)$", "blood HD \\1", meta$patient)
  
  # Clean up cell type names
  meta$cellType <- as.character(meta$cellType)
  meta$cellType[meta$cellType == "Prolif. T"] <- "Prolif T"
  meta$cellType[meta$cellType == "Other T"] <- "NK TOX+"
  meta$cellType[meta$cellType == "CD4 T p2"] <- "dn T"
  
  # Create main T-cell annotation
  meta$Tcells <- meta$cellType
  
  # Assign malignant clusters based on clustering results
  malignant_clusters <- c("0", "5", "7", "8", "9", "13", "14", "18", "20")
  meta$Tcells[meta$seurat_clusters %in% malignant_clusters] <- "Malignant TCM"
  
  # Consolidate cell types
  meta$Tcells[meta$Tcells %in% c("CD16- NK", "CD16+ NK")] <- "NK"
  meta$Tcells[meta$Tcells %in% c("CD4 TCM1", "CD4 TCM2")] <- "CD4 TCM"
  meta$Tcells[meta$Tcells %in% c("CD4 CTL")] <- "CTL"
  meta$Tcells[meta$Tcells %in% c("CD8 TEM")] <- "CTL"
  meta$Tcells[meta$Tcells %in% c("CD8 TOX+")] <- "NK TOX+"
  meta$Tcells[meta$Tcells %in% c("Tregs")] <- "Treg"
  
  # Add metadata back to Seurat object
  seurat_obj <- AddMetaData(seurat_obj, meta)
  
  # Remove clusters with too few cells
  seurat_obj <- subset(seurat_obj, subset = Tcells != "CD4 T p1")
  
  return(seurat_obj)
}

#' Annotate T-cell subtypes in skin samples
#'
#' @param seurat_obj Seurat object containing skin T-cells
#' @return Seurat object with T-cell annotations
annotate_skin_tcells <- function(seurat_obj) {
  
  cat("Annotating skin T-cell subtypes...\n")
  
  # Prepare metadata
  meta <- seurat_obj@meta.data
  
  # Assign T-cell types based on clustering resolution
  meta$Tcells[meta$RNA_snn_res.1 %in% c("2", "3", "6", "10", "11")] <- "Malignant TCM"
  meta$Tcells[meta$RNA_snn_res.1 %in% c("5", "9")] <- "Malignant IS"
  meta$Tcells[meta$RNA_snn_res.1 %in% c("4")] <- "Treg"
  meta$Tcells[meta$RNA_snn_res.1 %in% c("0")] <- "CD4 TRM"
  meta$Tcells[meta$RNA_snn_res.1 %in% c("1")] <- "CD4 TRMr"
  meta$Tcells[meta$RNA_snn_res.1 %in% c("7")] <- "CD8 T"
  meta$Tcells[meta$RNA_snn_res.1 %in% c("8")] <- "NK"
  
  # Clean up patient names
  meta$patient <- gsub("CTCL_(\\d+)", "skin \\1", meta$patient)
  
  # Add metadata back to Seurat object
  seurat_obj@meta.data <- meta
  
  # Remove clusters with too few cells
  clusters_to_remove <- c("12", "13", "14")
  seurat_obj <- subset(seurat_obj, subset = RNA_snn_res.1 %in% clusters_to_remove, invert = TRUE)
  
  return(seurat_obj)
}

#' Create color palette for T-cell types
#'
#' @param tissue_type Either "blood" or "skin"
#' @return Named vector of colors for each cell type
create_tcell_color_palette <- function(tissue_type = "blood") {
  
  if (tissue_type == "blood") {
    col_palette <- c(
      "Malignant TCM" = "#000080",
      "CD4 naive T" = "#1CBE4FFF",
      "CD4 TCM" = "#90AD1CFF", 
      "CD4 TEM" = "#1C8356FF",
      "Treg" = "#006600FF",
      "CTL" = "#664466FF",
      "CD8 naive T" = "#CC99CCFF",
      "CD8 TEM" = "#FF9966FF",
      "gdT" = "#BB6622FF",
      "Prolif T" = "#FFCC99FF",
      "NK TOX+" = "#CC6666FF",
      "NK" = "#882211FF"
    )
  } else if (tissue_type == "skin") {
    col_palette <- c(
      "Malignant TCM" = "#000080",
      "Malignant IS" = "#0011EEFF",
      "CD4 TRM" = "#66CC33FF",
      "CD4 TRMr" = "#BBAA55FF",
      "Treg" = "#006600FF",
      "CD8 T" = "#DD6644FF",
      "NK" = "#882211FF"
    )
  }
  
  return(col_palette)
}

#' Generate UMAP plot for T-cell types
#'
#' @param seurat_obj Seurat object with T-cell annotations
#' @param tissue_type Either "blood" or "skin"
#' @param reduction UMAP reduction to use
#' @return ggplot object
plot_tcell_umap <- function(seurat_obj, tissue_type = "blood", reduction = "umap") {
  
  # Get color palette
  col_palette <- create_tcell_color_palette(tissue_type)
  
  # Set cell type order for plotting
  if (tissue_type == "blood") {
    cell_order <- c("NK", "NK TOX+", "Prolif T", "gdT", "CD8 TEM", "CD8 naive T", 
                   "CTL", "Treg", "CD4 TEM", "CD4 TCM", "CD4 naive T", "Malignant TCM")
  } else {
    cell_order <- c("NK", "CD8 T", "Treg", "CD4 TRM", "CD4 TRMr", "Malignant IS", "Malignant TCM")
  }
  
  # Set identities
  Idents(seurat_obj) <- seurat_obj$Tcells
  
  # Create UMAP plot
  umap_plot <- do_DimPlot(
    seurat_obj,
    pt.size = 0.3,
    reduction = reduction,
    plot.axes = TRUE,
    label = FALSE,
    label.box = FALSE,
    label.size = 4,
    repel = TRUE,
    legend.position = "right",
    plot.title = paste(stringr::str_to_title(tissue_type), "T-cells"),
    order = cell_order,
    colors.use = col_palette
  )
  
  return(umap_plot)
}

#' Generate patient-wise UMAP plot
#'
#' @param seurat_obj Seurat object with patient annotations
#' @param reduction UMAP reduction to use
#' @return ggplot object
plot_patient_umap <- function(seurat_obj, reduction = "umap") {
  
  # Create color palette for patients
  n_patients <- length(unique(seurat_obj$patient))
  colors <- paletteer_d("trekcolors::lcars_series")
  colors <- colors[1:n_patients]
  col_palette <- setNames(colors, unique(seurat_obj$patient))
  
  # Create split UMAP plot by patient
  patient_plot <- do_DimPlot(
    seurat_obj,
    pt.size = 0.3,
    reduction = reduction,
    group.by = "patient",
    split.by = "patient",
    plot.axes = TRUE,
    label = FALSE,
    label.box = FALSE,
    label.size = 4,
    repel = TRUE,
    legend.position = "right",
    plot.title = "T-cells by Patient",
    colors.use = col_palette
  )
  
  return(patient_plot)
}

#' Create stacked barplot of T-cell proportions
#'
#' @param seurat_obj Seurat object with T-cell annotations
#' @param tissue_type Either "blood" or "skin"
#' @return ggplot object
plot_tcell_proportions <- function(seurat_obj, tissue_type = "blood") {
  
  # Extract metadata
  metadata <- seurat_obj@meta.data
  
  # Get color palette
  col_palette <- create_tcell_color_palette(tissue_type)
  
  # Set factor levels for consistent ordering
  if (tissue_type == "blood") {
    cell_levels <- c("Malignant TCM", "CD4 naive T", "CD4 TCM", "CD4 TEM", "Treg", 
                    "CTL", "CD8 naive T", "CD8 TEM", "gdT", "Prolif T", "NK TOX+", "NK")
  } else {
    cell_levels <- c("Malignant TCM", "Malignant IS", "CD4 TRM", "CD4 TRMr", 
                    "Treg", "CD8 T", "NK")
  }
  
  metadata$Tcells <- factor(metadata$Tcells, levels = cell_levels)
  
  # Count cells per patient and T-cell type
  cell_counts <- metadata %>%
    group_by(patient, Tcells) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(patient) %>%
    mutate(proportion = count / sum(count))
  
  # Create stacked barplot
  prop_plot <- ggplot(cell_counts, aes(x = patient, y = proportion, fill = Tcells)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = col_palette) +
    labs(x = "Patient", y = "Proportion", fill = "T Cell Type",
         title = paste(stringr::str_to_title(tissue_type), "T-cell Proportions")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  return(prop_plot)
}

#' Create dotplot of marker genes
#'
#' @param seurat_obj Seurat object with T-cell annotations
#' @param tissue_type Either "blood" or "skin"
#' @return ggplot object
plot_tcell_markers <- function(seurat_obj, tissue_type = "blood") {
  
  # Define marker genes based on tissue type
  if (tissue_type == "blood") {
    markers <- c(
      "CD3E", "CD3D", "CD4",  # T-cell markers
      "TOX", "TSHZ2", "SESN3", "DNM3",  # Malignant markers
      "TCF7", "LEF1", "CCR7",  # Naive/central memory
      "ITGB1", "AQP3",  # TCM
      "IL32", "CCL5", "GZMK",  # TEM
      "FOXP3", "IL2RA",  # Treg
      "CD8A", "CD8B", "GZMB", "NKG7",  # CTL/CD8
      "TRGC1", "TRGC2",  # gdT
      "KLRF1",  # NK TOX+
      "NCAM1", "FCGR3A",  # NK
      "PCNA", "MKI67"  # Proliferation
    )
    
    # Set factor levels
    cell_levels <- rev(c("Malignant TCM", "CD4 naive T", "CD4 TCM", "CD4 TEM", 
                        "Treg", "CTL", "CD8 naive T", "CD8 TEM", "gdT", 
                        "NK TOX+", "NK", "Prolif T"))
    
  } else {
    markers <- c(
      "CD3E", "CD3D", "CD4",  # T-cell markers
      "TOX", "TSHZ2", "SESN3",  # Malignancy
      "DNM3", "CD28",  # Malignant subtypes
      "IL2RA", "IKZF2", "ICOS", "CTLA4",  # Immunosuppressive
      "FOXP3", "TIGIT", "CXCR4",  # Treg
      "CD69",  # TRM
      "LGALS3", "IL32",  # Reactive
      "GZMK", "CD8A",  # CD8 T
      "NKG7", "GZMB", "GNLY", "NCAM1"  # NK
    )
    
    # Set factor levels
    cell_levels <- c("NK", "CD8 T", "CD4 TRMr", "CD4 TRM", "Treg", 
                    "Malignant IS", "Malignant TCM")
  }
  
  # Set cell type factor levels
  seurat_obj$Tcells <- factor(seurat_obj$Tcells, levels = cell_levels)
  Idents(seurat_obj) <- seurat_obj$Tcells
  
  # Set default assay
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Create dotplot
  dot_plot <- DotPlot(seurat_obj, group.by = "Tcells", features = markers) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1)) +
    scale_colour_gradient2(low = "lightgreen", mid = "lightgrey", high = "coral") +
    labs(title = paste(stringr::str_to_title(tissue_type), "T-cell Markers"))
  
  return(dot_plot)
}

#' Main T-cell integration workflow
#'
#' @param blood_tcells_path Path to blood T-cells Seurat object
#' @param skin_tcells_path Path to skin T-cells Seurat object
#' @param output_dir Output directory for results
#' @param remove_tcr_genes Whether to remove TCR genes before clustering
#' @return List containing processed Seurat objects and plots
process_tcell_integration <- function(blood_tcells_path, skin_tcells_path,
                                     output_dir = "Tcell_integration_output",
                                     remove_tcr_genes = TRUE) {
  
  cat("=== Starting T-cell Integration Workflow ===\n")
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  results <- list()
  
  # Process blood T-cells
  if (!is.null(blood_tcells_path) && file.exists(blood_tcells_path)) {
    cat("\n--- Processing Blood T-cells ---\n")
    
    blood_tcells <- readRDS(blood_tcells_path)
    DefaultAssay(blood_tcells) <- "RNA"
    
    # Remove TCR genes and re-cluster if requested
    if (remove_tcr_genes) {
      blood_tcells <- remove_tcr_genes_and_recluster(blood_tcells)
    }
    
    # Annotate T-cell subtypes
    blood_tcells <- annotate_blood_tcells(blood_tcells)
    
    # Generate plots
    blood_plots <- list()
    blood_plots$umap_celltypes <- plot_tcell_umap(blood_tcells, "blood")
    blood_plots$umap_patients <- plot_patient_umap(blood_tcells)
    blood_plots$proportions <- plot_tcell_proportions(blood_tcells, "blood")
    blood_plots$markers <- plot_tcell_markers(blood_tcells, "blood")
    
    # Save blood results
    saveRDS(blood_tcells, file.path(output_dir, "blood_tcells_annotated.rds"))
    results$blood_tcells <- blood_tcells
    results$blood_plots <- blood_plots
  }
  
  # Process skin T-cells
  if (!is.null(skin_tcells_path) && file.exists(skin_tcells_path)) {
    cat("\n--- Processing Skin T-cells ---\n")
    
    skin_tcells <- readRDS(skin_tcells_path)
    DefaultAssay(skin_tcells) <- "RNA"
    
    # Remove TCR genes and re-cluster if requested
    if (remove_tcr_genes) {
      skin_tcells <- remove_tcr_genes_and_recluster(skin_tcells)
    }
    
    # Annotate T-cell subtypes
    skin_tcells <- annotate_skin_tcells(skin_tcells)
    
    # Generate plots
    skin_plots <- list()
    skin_plots$umap_celltypes <- plot_tcell_umap(skin_tcells, "skin")
    skin_plots$umap_patients <- plot_patient_umap(skin_tcells)
    skin_plots$proportions <- plot_tcell_proportions(skin_tcells, "skin")
    skin_plots$markers <- plot_tcell_markers(skin_tcells, "skin")
    
    # Save skin results
    saveRDS(skin_tcells, file.path(output_dir, "skin_tcells_annotated.rds"))
    results$skin_tcells <- skin_tcells
    results$skin_plots <- skin_plots
  }
  
  # Save all plots
  all_plots <- c(results$blood_plots, results$skin_plots)
  for (plot_name in names(all_plots)) {
    ggsave(
      filename = file.path(output_dir, paste0(plot_name, ".pdf")),
      plot = all_plots[[plot_name]],
      width = 12, height = 8
    )
  }
  
  cat("\n=== T-cell Integration Workflow Completed ===\n")
  cat("Results saved to:", output_dir, "\n")
  
  return(results)
}

# Example usage (uncomment and modify paths as needed):
# 
# # Process T-cell integration
# results <- process_tcell_integration(
#   blood_tcells_path = "path/to/BLOOD_Tcells.rds",
#   skin_tcells_path = "path/to/SKIN_Tcells.rds",
#   output_dir = "Tcell_integration_results",
#   remove_tcr_genes = TRUE
# )

cat("T-cell integration script loaded successfully!\n")
cat("Use process_tcell_integration() to run the complete analysis pipeline.\n")
