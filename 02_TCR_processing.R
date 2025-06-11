#!/usr/bin/env Rscript

#' TCR Processing and Analysis Pipeline
#' 
#' This script processes T-cell receptor (TCR) sequencing data from CellRanger output
#' and performs comprehensive TCR repertoire analysis including clonotype identification,
#' diversity analysis, and integration with single-cell RNA-seq data.
#'
#' Based on the CTCL single-cell memory manuscript analysis pipeline

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(scRepertoire)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(kableExtra)
  library(SignatuR)
  library(Trex)
  library(pals)
  library(future)
})

# Set up parallel processing
param <- list(cores = 4, ram = 60, seed = 38)
plan("multicore", workers = param$cores)
set.seed(param$seed)
options(future.globals.maxSize = param$ram * 1024^3)

#' Load TCR data from CellRanger output directories
#'
#' @param sample_paths Vector of paths to CellRanger output directories
#' @param sample_names Vector of sample names corresponding to paths
#' @return List of TCR data frames
load_tcr_data <- function(sample_paths, sample_names) {
  cat("Loading TCR data for", length(sample_paths), "samples...\n")
  
  tcr_list <- list()
  
  for (i in seq_along(sample_paths)) {
    tcr_file <- file.path(sample_paths[i], "filtered_contig_annotations.csv")
    
    if (file.exists(tcr_file)) {
      tcr_data <- read.csv(tcr_file, stringsAsFactors = FALSE)
      tcr_list[[sample_names[i]]] <- tcr_data
      cat("✓ Loaded TCR data for sample:", sample_names[i], 
          "- contigs:", nrow(tcr_data), "\n")
    } else {
      warning("✗ TCR file not found for sample:", sample_names[i], "\n")
    }
  }
  
  return(tcr_list)
}

#' Combine and process TCR data using scRepertoire
#'
#' @param tcr_list List of TCR data frames from load_tcr_data
#' @param remove_na Remove cells without TCR information
#' @param remove_multi Remove cells with multiple productive chains
#' @return Combined TCR object from scRepertoire
combine_tcr_data <- function(tcr_list, remove_na = FALSE, remove_multi = FALSE) {
  cat("Combining TCR data across", length(tcr_list), "samples...\n")
  
  # Combine using scRepertoire
  combined_tcr <- combineTCR(
    tcr_list, 
    samples = names(tcr_list),
    removeNA = remove_na, 
    removeMulti = remove_multi, 
    filterMulti = FALSE
  )
  
  # Add clone size categories
  combined_tcr <- categorize_clone_sizes(combined_tcr)
  
  return(combined_tcr)
}

#' Categorize clones by expansion level
#'
#' @param combined_tcr Combined TCR object from scRepertoire
#' @return Combined TCR object with clone size categories
categorize_clone_sizes <- function(combined_tcr) {
  cat("Categorizing clones by expansion level...\n")
  
  for (i in seq_along(combined_tcr)) {
    combined_tcr[[i]]$cloneSize <- case_when(
      combined_tcr[[i]]$Frequency == 1 ~ "Single",
      combined_tcr[[i]]$Frequency >= 2 & combined_tcr[[i]]$Frequency <= 5 ~ "Small",
      combined_tcr[[i]]$Frequency >= 6 & combined_tcr[[i]]$Frequency <= 20 ~ "Medium", 
      combined_tcr[[i]]$Frequency >= 21 & combined_tcr[[i]]$Frequency <= 50 ~ "Large",
      combined_tcr[[i]]$Frequency > 50 ~ "Hyperexpanded"
    )
    
    # Set factor levels for consistent plotting
    combined_tcr[[i]]$cloneSize <- factor(
      combined_tcr[[i]]$cloneSize,
      levels = c("Single", "Small", "Medium", "Large", "Hyperexpanded")
    )
  }
  
  return(combined_tcr)
}

#' Add TCR data to Seurat object
#'
#' @param seurat_obj Seurat object containing single-cell RNA-seq data
#' @param combined_tcr Combined TCR object from scRepertoire
#' @param group_by Grouping variable for combining (default: "sample")
#' @return Seurat object with TCR information added
integrate_tcr_with_seurat <- function(seurat_obj, combined_tcr, group_by = "sample") {
  cat("Integrating TCR data with Seurat object...\n")
  
  # Update cell barcodes to match TCR data format if needed
  if (!is.null(seurat_obj$sample)) {
    cell_barcodes <- Cells(seurat_obj)
    cell_barcodes <- gsub("_.", "", cell_barcodes)
    cell_barcodes <- paste0(seurat_obj$sample, "_", cell_barcodes)
    seurat_obj <- RenameCells(seurat_obj, new.names = cell_barcodes)
  }
  
  # Combine with expression data
  seurat_obj <- combineExpression(
    combined_tcr, 
    seurat_obj, 
    cloneCall = "strict", 
    group.by = group_by, 
    proportion = FALSE, 
    cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 50, Hyperexpanded = 100)
  )
  
  return(seurat_obj)
}

#' Analyze clonal diversity using multiple metrics
#'
#' @param seurat_obj Seurat object with TCR data
#' @param group_by Grouping variable for analysis
#' @return ggplot object showing diversity metrics
analyze_clonal_diversity <- function(seurat_obj, group_by = "sample") {
  cat("Analyzing clonal diversity...\n")
  
  diversity_plot <- clonalDiversity(
    seurat_obj,
    cloneCall = "gene",
    group.by = group_by
  ) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Clonal Diversity by", group_by))
  
  return(diversity_plot)
}

#' Analyze clonal homeostasis
#'
#' @param seurat_obj Seurat object with TCR data
#' @param group_by Grouping variable for analysis
#' @return ggplot object showing clonal homeostasis
analyze_clonal_homeostasis <- function(seurat_obj, group_by = "sample") {
  cat("Analyzing clonal homeostasis...\n")
  
  homeostasis_plot <- clonalHomeostasis(
    seurat_obj, 
    cloneCall = "gene",
    group.by = group_by
  ) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Clonal Space Homeostasis by", group_by))
  
  return(homeostasis_plot)
}

#' Analyze clonal overlap between groups
#'
#' @param seurat_obj Seurat object with TCR data
#' @param group_by Grouping variable for analysis
#' @param method Overlap calculation method
#' @return ggplot object showing clonal overlap heatmap
analyze_clonal_overlap <- function(seurat_obj, group_by = "sample", method = "morisita") {
  cat("Analyzing clonal overlap using", method, "index...\n")
  
  overlap_plot <- clonalOverlap(
    seurat_obj,
    cloneCall = "strict",
    group.by = group_by,
    method = method
  ) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Clonal Overlap (", method, ") by", group_by))
  
  return(overlap_plot)
}

#' Visualize TCR gene usage patterns
#'
#' @param seurat_obj Seurat object with TCR data
#' @param chains Vector of TCR chains to analyze (e.g., c("TRAV", "TRBV"))
#' @return List of ggplot objects showing gene usage
visualize_tcr_gene_usage <- function(seurat_obj, chains = c("TRAV", "TRBV")) {
  cat("Visualizing TCR gene usage for chains:", paste(chains, collapse = ", "), "\n")
  
  plots <- list()
  
  for (chain in chains) {
    plots[[chain]] <- vizGenes(
      seurat_obj, 
      x.axis = chain,
      y.axis = NULL,
      plot = "barplot",  
      scale = TRUE
    ) + 
      labs(title = paste(chain, "Gene Usage")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  return(plots)
}

#' Perform antigen annotation using public databases
#'
#' @param seurat_obj Seurat object with TCR data
#' @param chains TCR chains to annotate
#' @param edit_distance Maximum edit distance for matching
#' @param min_cells Minimum number of cells for keeping antigen annotation
#' @return Seurat object with antigen annotations
annotate_antigens <- function(seurat_obj, chains = c("TRA", "TRB"), 
                             edit_distance = 0, min_cells = 5) {
  cat("Annotating antigens for chains:", paste(chains, collapse = ", "), "\n")
  
  for (chain in chains) {
    seurat_obj <- annotateDB(
      seurat_obj, 
      chains = chain, 
      edit.distance = edit_distance
    )
    
    # Filter rare epitopes
    epitope_col <- paste0(chain, "_Epitope.species")
    if (epitope_col %in% colnames(seurat_obj@meta.data)) {
      epitope_counts <- table(seurat_obj@meta.data[[epitope_col]])
      rare_epitopes <- names(epitope_counts[epitope_counts < min_cells])
      seurat_obj@meta.data[[epitope_col]][
        seurat_obj@meta.data[[epitope_col]] %in% rare_epitopes
      ] <- "NA"
    }
  }
  
  return(seurat_obj)
}

#' Generate comprehensive TCR analysis plots
#'
#' @param seurat_obj Seurat object with TCR data
#' @param group_by Primary grouping variable
#' @param output_dir Directory to save plots
#' @return List of generated plots
generate_tcr_plots <- function(seurat_obj, group_by = "sample", output_dir = NULL) {
  cat("Generating comprehensive TCR analysis plots...\n")
  
  plots <- list()
  
  # Clone size distribution
  plots$clone_occupy <- clonalOccupy(seurat_obj, x.axis = group_by) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Clone Size Distribution by", group_by))
  
  # Clone quantification
  plots$clone_quant <- clonalQuant(
    seurat_obj, 
    cloneCall = "strict", 
    chain = "both", 
    scale = TRUE
  ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Clonal Quantification (Scaled)")
  
  # Clone abundance
  plots$clone_abundance <- clonalAbundance(
    seurat_obj, 
    cloneCall = "strict", 
    scale = FALSE
  ) +
    labs(title = "Clonal Abundance")
  
  # Diversity analysis
  plots$diversity <- analyze_clonal_diversity(seurat_obj, group_by)
  
  # Homeostasis analysis
  plots$homeostasis <- analyze_clonal_homeostasis(seurat_obj, group_by)
  
  # Overlap analysis
  plots$overlap <- analyze_clonal_overlap(seurat_obj, group_by)
  
  # TCR gene usage
  gene_plots <- visualize_tcr_gene_usage(seurat_obj)
  plots <- c(plots, gene_plots)
  
  # Save plots if output directory specified
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    for (plot_name in names(plots)) {
      ggsave(
        filename = file.path(output_dir, paste0("TCR_", plot_name, ".pdf")),
        plot = plots[[plot_name]],
        width = 10, height = 8
      )
    }
    cat("Plots saved to:", output_dir, "\n")
  }
  
  return(plots)
}

#' Main TCR processing workflow
#'
#' @param sample_paths Vector of paths to CellRanger output directories
#' @param sample_names Vector of sample names
#' @param seurat_obj Seurat object (optional, for integration)
#' @param output_dir Output directory for results
#' @param perform_antigen_annotation Whether to perform antigen annotation
#' @return List containing processed TCR data and analysis results
process_tcr_workflow <- function(sample_paths, sample_names, seurat_obj = NULL, 
                                output_dir = "TCR_analysis_output",
                                perform_antigen_annotation = TRUE) {
  
  cat("=== Starting TCR Processing Workflow ===\n")
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Step 1: Load TCR data
  tcr_list <- load_tcr_data(sample_paths, sample_names)
  
  # Step 2: Combine TCR data
  combined_tcr <- combine_tcr_data(tcr_list)
  
  # Step 3: Integrate with Seurat object if provided
  if (!is.null(seurat_obj)) {
    seurat_obj <- integrate_tcr_with_seurat(seurat_obj, combined_tcr)
    
    # Step 4: Perform antigen annotation if requested
    if (perform_antigen_annotation) {
      seurat_obj <- annotate_antigens(seurat_obj)
    }
    
    # Step 5: Generate analysis plots
    plots <- generate_tcr_plots(seurat_obj, output_dir = output_dir)
    
    # Save processed Seurat object
    saveRDS(seurat_obj, file.path(output_dir, "seurat_with_tcr.rds"))
  } else {
    plots <- NULL
  }
  
  # Save TCR data
  saveRDS(combined_tcr, file.path(output_dir, "combined_tcr.rds"))
  
  cat("=== TCR Processing Workflow Completed ===\n")
  
  return(list(
    combined_tcr = combined_tcr,
    seurat_obj = seurat_obj,
    plots = plots
  ))
}

# Example usage (uncomment and modify paths as needed):
# 
# # Define sample paths and names
# sample_paths <- c(
#   "/path/to/sample1/outs",
#   "/path/to/sample2/outs"
# )
# sample_names <- c("Sample1", "Sample2")
# 
# # Load existing Seurat object (optional)
# seurat_obj <- readRDS("path/to/seurat_object.rds")
# 
# # Run TCR processing workflow
# results <- process_tcr_workflow(
#   sample_paths = sample_paths,
#   sample_names = sample_names,
#   seurat_obj = seurat_obj,
#   output_dir = "TCR_analysis_results"
# )

cat("TCR processing script loaded successfully!\n")
cat("Use process_tcr_workflow() to run the complete analysis pipeline.\n")
