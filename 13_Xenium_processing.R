#!/usr/bin/env Rscript

#' Xenium Spatial Transcriptomics Processing Pipeline
#' 
#' This script processes Xenium spatial transcriptomics data including
#' preprocessing, segmentation, RCTD deconvolution, and analysis.
#' It reproduces the spatial analysis from the CTCL single-cell memory manuscript.

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(future)
  library(fs)
  library(BiocParallel)
  library(qs)
  library(ComplexHeatmap)
  library(scales)
})

# Set up parallel processing
param <- list(cores = 20, ram = 60, seed = 38)
BPPARAM <- MulticoreParam(workers = param$cores)
register(BPPARAM)
plan("multicore", workers = param$cores)
set.seed(param$seed)
options(future.globals.maxSize = param$ram * 1024^3)

#' Load Xenium data from multiple regions
#'
#' @param base_dir Base directory containing Xenium output folders
#' @param pattern Pattern to match region directories
#' @return List of Seurat objects for each region
load_xenium_regions <- function(base_dir, pattern = "output-XETG00404") {
  
  cat("Loading Xenium regions from:", base_dir, "\n")
  
  # List all directories matching the pattern
  dirs <- dir_ls(base_dir, recurse = FALSE, type = "directory")
  dirs <- dirs[grep(pattern, basename(dirs))]
  
  cat("Found", length(dirs), "potential Xenium directories\n")
  
  regions <- list()
  
  for (dir in dirs) {
    region_name <- str_extract(basename(dir), "reg\\d+")
    
    if (!is.na(region_name)) {
      cat("Loading region:", region_name, "from", dir, "\n")
      
      tryCatch({
        region_obj <- LoadXenium(
          data.dir = dir,
          fov = region_name,
          assay = "Xenium"
        )
        
        regions[[region_name]] <- region_obj
        cat("✓ Successfully loaded", region_name, "with", ncol(region_obj), "cells\n")
        
      }, error = function(e) {
        cat("✗ Error loading", region_name, ":", e$message, "\n")
      })
    }
  }
  
  cat("Successfully loaded", length(regions), "regions\n")
  return(regions)
}

#' Preprocess Xenium data with quality control
#'
#' @param xenium_obj Seurat object with Xenium data
#' @param min_features Minimum number of features per cell
#' @param max_features Maximum number of features per cell
#' @param min_counts Minimum number of counts per cell
#' @param max_counts Maximum number of counts per cell
#' @return Preprocessed Seurat object
preprocess_xenium_data <- function(xenium_obj, min_features = 5, max_features = Inf,
                                  min_counts = 10, max_counts = Inf) {
  
  cat("Preprocessing Xenium data...\n")
  cat("Initial cells:", ncol(xenium_obj), "\n")
  
  # Calculate QC metrics
  xenium_obj[["percent.mt"]] <- PercentageFeatureSet(xenium_obj, pattern = "^MT-")
  xenium_obj[["nCount_Xenium"]] <- colSums(xenium_obj[["Xenium"]]$counts)
  xenium_obj[["nFeature_Xenium"]] <- colSums(xenium_obj[["Xenium"]]$counts > 0)
  
  # Apply quality control filters
  xenium_obj <- subset(
    xenium_obj,
    subset = nFeature_Xenium >= min_features & 
             nFeature_Xenium <= max_features &
             nCount_Xenium >= min_counts &
             nCount_Xenium <= max_counts
  )
  
  cat("Cells after QC:", ncol(xenium_obj), "\n")
  
  return(xenium_obj)
}

#' Create RCTD reference from single-cell RNA-seq data
#'
#' @param reference_seurat Seurat object to use as reference
#' @param cell_type_column Column name containing cell type annotations
#' @param min_cells Minimum number of cells per cell type
#' @return RCTD reference object
create_rctd_reference <- function(reference_seurat, cell_type_column = "celltype", 
                                 min_cells = 25) {
  
  cat("Creating RCTD reference from single-cell data...\n")
  
  # Get expression data and metadata
  counts <- reference_seurat[["RNA"]]$counts
  metadata <- reference_seurat@meta.data
  
  # Filter cell types with sufficient cells
  cell_type_counts <- table(metadata[[cell_type_column]])
  valid_cell_types <- names(cell_type_counts[cell_type_counts >= min_cells])
  
  cat("Keeping", length(valid_cell_types), "cell types with >=", min_cells, "cells\n")
  
  # Filter cells and metadata
  valid_cells <- metadata[[cell_type_column]] %in% valid_cell_types
  counts_filtered <- counts[, valid_cells]
  metadata_filtered <- metadata[valid_cells, ]
  
  # Create cell type factor
  cell_types <- factor(metadata_filtered[[cell_type_column]])
  names(cell_types) <- rownames(metadata_filtered)
  
  # Filter genes (keep genes expressed in at least 3 cells)
  gene_filter <- rowSums(counts_filtered > 0) >= 3
  counts_filtered <- counts_filtered[gene_filter, ]
  
  cat("Reference contains", nrow(counts_filtered), "genes and", ncol(counts_filtered), "cells\n")
  
  # Create RCTD reference
  reference <- Reference(
    counts = counts_filtered,
    cell_types = cell_types,
    n_max_cells = 10000  # Subsample if needed for memory efficiency
  )
  
  return(reference)
}

#' Run RCTD deconvolution on Xenium data
#'
#' @param xenium_obj Seurat object with Xenium data
#' @param reference RCTD reference object
#' @param region_name Name of the region being processed
#' @param umi_min Minimum UMI threshold for RCTD
#' @return List containing RCTD results and processed data
run_rctd_deconvolution <- function(xenium_obj, reference, region_name, umi_min = 20) {
  
  cat("Running RCTD deconvolution for", region_name, "\n")
  
  # Get spatial coordinates and expression data
  coords <- GetTissueCoordinates(xenium_obj)
  rownames(coords) <- coords$cell
  coords <- coords[, c("x", "y")]
  
  counts_query <- xenium_obj[["Xenium"]]$counts
  
  # Ensure common cells between counts and coordinates
  common_cells <- intersect(colnames(counts_query), rownames(coords))
  counts_query <- counts_query[, common_cells]
  coords <- coords[common_cells, ]
  
  cat("Processing", length(common_cells), "cells\n")
  
  # Calculate nUMI
  nUMI <- colSums(counts_query)
  
  # Create SpatialRNA object
  query <- SpatialRNA(coords, counts_query, nUMI)
  
  # Create and run RCTD
  cat("Creating RCTD object...\n")
  RCTD <- create.RCTD(query, reference, max_cores = param$cores, UMI_min = umi_min)
  
  cat("Running RCTD in doublet mode...\n")
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  
  # Process results
  cat("Processing RCTD results...\n")
  weights <- RCTD@results$weights
  norm_weights <- normalize_weights(weights)
  
  # Get cell type assignments (highest weight)
  cell_groups <- data.frame(
    cell_id = rownames(norm_weights),
    cell_type = colnames(norm_weights)[max.col(norm_weights)],
    max_weight = apply(norm_weights, 1, max),
    stringsAsFactors = FALSE
  )
  
  # Add second highest for doublet detection
  cell_groups$second_max_weight <- apply(norm_weights, 1, function(x) sort(x, decreasing = TRUE)[2])
  cell_groups$is_doublet <- cell_groups$second_max_weight > 0.3
  
  return(list(
    RCTD = RCTD,
    weights = weights,
    norm_weights = norm_weights,
    cell_groups = cell_groups,
    query = query
  ))
}

#' Add RCTD results to Seurat object
#'
#' @param xenium_obj Seurat object
#' @param rctd_results RCTD results list
#' @return Seurat object with RCTD annotations
add_rctd_to_seurat <- function(xenium_obj, rctd_results) {
  
  cat("Adding RCTD results to Seurat object...\n")
  
  # Create metadata data frame
  cell_metadata <- rctd_results$cell_groups
  rownames(cell_metadata) <- cell_metadata$cell_id
  
  # Match with Seurat cell names
  seurat_cells <- Cells(xenium_obj)
  matched_metadata <- cell_metadata[seurat_cells, ]
  
  # Add to Seurat object
  xenium_obj$rctd_cell_type <- matched_metadata$cell_type
  xenium_obj$rctd_max_weight <- matched_metadata$max_weight
  xenium_obj$rctd_second_weight <- matched_metadata$second_max_weight
  xenium_obj$rctd_is_doublet <- matched_metadata$is_doublet
  
  # Add individual cell type weights
  norm_weights_df <- as.data.frame(rctd_results$norm_weights)
  norm_weights_df$cell_id <- rownames(norm_weights_df)
  norm_weights_matched <- norm_weights_df[seurat_cells, ]
  
  for (ct in colnames(rctd_results$norm_weights)) {
    xenium_obj@meta.data[[paste0("rctd_weight_", ct)]] <- norm_weights_matched[[ct]]
  }
  
  return(xenium_obj)
}

#' Generate spatial plots for RCTD results
#'
#' @param xenium_obj Seurat object with RCTD results
#' @param region_name Name of the region
#' @return List of ggplot objects
plot_spatial_results <- function(xenium_obj, region_name) {
  
  cat("Generating spatial plots for", region_name, "\n")
  
  plots <- list()
  
  # Get coordinates
  coords <- GetTissueCoordinates(xenium_obj)
  
  # Plot 1: Cell type assignments
  plot_data <- cbind(coords, xenium_obj@meta.data)
  
  plots$cell_types <- ggplot(plot_data, aes(x = x, y = y, color = rctd_cell_type)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_color_viridis_d(name = "Cell Type") +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold")
    ) +
    labs(title = paste("Cell Type Assignments -", region_name)) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  # Plot 2: Assignment confidence (max weight)
  plots$confidence <- ggplot(plot_data, aes(x = x, y = y, color = rctd_max_weight)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_color_viridis_c(name = "Max Weight") +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold")
    ) +
    labs(title = paste("Assignment Confidence -", region_name))
  
  # Plot 3: Doublets
  plots$doublets <- ggplot(plot_data, aes(x = x, y = y, color = rctd_is_doublet)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_color_manual(values = c("FALSE" = "gray80", "TRUE" = "red"), 
                      name = "Doublet") +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold")
    ) +
    labs(title = paste("Predicted Doublets -", region_name))
  
  # Plot 4: Malignant T-cell types (if present)
  malignant_types <- c("Malignant TCM", "Malignant TRM", "malignant_T")
  malignant_present <- any(malignant_types %in% plot_data$rctd_cell_type)
  
  if (malignant_present) {
    plot_data$is_malignant <- plot_data$rctd_cell_type %in% malignant_types
    
    plots$malignant <- ggplot(plot_data, aes(x = x, y = y, color = is_malignant)) +
      geom_point(size = 0.5, alpha = 0.7) +
      scale_color_manual(values = c("FALSE" = "gray80", "TRUE" = "darkred"), 
                        name = "Malignant") +
      theme_void() +
      theme(
        legend.position = "right",
        plot.title = element_text(size = 14, face = "bold")
      ) +
      labs(title = paste("Malignant T-cells -", region_name))
  }
  
  return(plots)
}

#' Generate summary statistics for RCTD results
#'
#' @param xenium_obj Seurat object with RCTD results
#' @param region_name Name of the region
#' @return Data frame with summary statistics
generate_spatial_summary <- function(xenium_obj, region_name) {
  
  # Cell type proportions
  cell_type_props <- table(xenium_obj$rctd_cell_type)
  cell_type_props <- cell_type_props / sum(cell_type_props)
  
  # Doublet rate
  doublet_rate <- mean(xenium_obj$rctd_is_doublet, na.rm = TRUE)
  
  # Mean confidence
  mean_confidence <- mean(xenium_obj$rctd_max_weight, na.rm = TRUE)
  
  # Create summary data frame
  summary_df <- data.frame(
    Region = region_name,
    Total_Cells = ncol(xenium_obj),
    Unique_Cell_Types = length(unique(xenium_obj$rctd_cell_type)),
    Doublet_Rate = round(doublet_rate, 3),
    Mean_Confidence = round(mean_confidence, 3),
    stringsAsFactors = FALSE
  )
  
  # Add cell type proportions as separate columns
  for (ct in names(cell_type_props)) {
    summary_df[[paste0("Prop_", gsub("[^A-Za-z0-9]", "_", ct))]] <- round(cell_type_props[ct], 3)
  }
  
  return(summary_df)
}

#' Analyze malignant T-cell subpopulations in spatial data
#'
#' @param xenium_obj Seurat object with RCTD results
#' @param region_name Name of the region
#' @return List containing malignant analysis results
analyze_malignant_subpopulations <- function(xenium_obj, region_name) {
  
  cat("Analyzing malignant subpopulations in", region_name, "\n")
  
  # Define malignant cell types
  malignant_types <- c("Malignant TCM", "Malignant TRM", "malignant_T")
  
  # Check if malignant cells are present
  malignant_present <- any(malignant_types %in% xenium_obj$rctd_cell_type)
  
  if (!malignant_present) {
    cat("No malignant cell types found in", region_name, "\n")
    return(NULL)
  }
  
  # Filter for malignant cells
  malignant_cells <- xenium_obj$rctd_cell_type %in% malignant_types
  malignant_obj <- subset(xenium_obj, cells = Cells(xenium_obj)[malignant_cells])
  
  # Calculate proportions
  malignant_props <- table(malignant_obj$rctd_cell_type)
  malignant_props <- malignant_props / sum(malignant_props)
  
  # Test for significant enrichment of mTRM vs mTCM
  if (length(malignant_props) >= 2) {
    # Binomial test for mTRM enrichment
    mtcm_count <- sum(grepl("TCM|malignant_T", malignant_obj$rctd_cell_type))
    mtrm_count <- sum(grepl("TRM", malignant_obj$rctd_cell_type))
    total_malignant <- mtcm_count + mtrm_count
    
    if (total_malignant > 0) {
      binom_test <- binom.test(mtrm_count, total_malignant, p = 0.5)
      mtrm_enriched <- binom_test$p.value < 0.05 && mtrm_count > mtcm_count
    } else {
      mtrm_enriched <- FALSE
      binom_test <- NULL
    }
  } else {
    mtrm_enriched <- FALSE
    binom_test <- NULL
  }
  
  return(list(
    region = region_name,
    malignant_proportions = malignant_props,
    mtrm_enriched = mtrm_enriched,
    binom_test = binom_test,
    total_malignant_cells = ncol(malignant_obj)
  ))
}

#' Main Xenium processing workflow
#'
#' @param base_dir Base directory containing Xenium output
#' @param reference_path Path to reference Seurat object or RCTD reference
#' @param output_dir Output directory for results
#' @param custom_genes_file Path to custom gene panel file (optional)
#' @return List containing processed data and results
process_xenium_workflow <- function(base_dir, reference_path, 
                                   output_dir = "Xenium_processing_output",
                                   custom_genes_file = NULL) {
  
  cat("=== Starting Xenium Processing Workflow ===\n")
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load reference
  cat("Loading reference data...\n")
  if (str_ends(reference_path, ".qs")) {
    reference <- qread(reference_path, nthreads = param$cores)
  } else if (str_ends(reference_path, ".rds")) {
    reference_seurat <- readRDS(reference_path)
    reference <- create_rctd_reference(reference_seurat)
  } else {
    stop("Reference file must be .qs (RCTD reference) or .rds (Seurat object)")
  }
  
  # Load Xenium regions
  regions <- load_xenium_regions(base_dir)
  
  if (length(regions) == 0) {
    stop("No Xenium regions found in", base_dir)
  }
  
  # Process each region
  results <- list()
  all_plots <- list()
  summary_stats <- data.frame()
  malignant_analyses <- list()
  
  for (region_name in names(regions)) {
    cat("\n--- Processing", region_name, "---\n")
    
    # Preprocess data
    regions[[region_name]] <- preprocess_xenium_data(regions[[region_name]])
    
    # Run RCTD
    rctd_results <- run_rctd_deconvolution(
      regions[[region_name]], 
      reference, 
      region_name
    )
    
    # Add results to Seurat object
    regions[[region_name]] <- add_rctd_to_seurat(regions[[region_name]], rctd_results)
    
    # Generate plots
    plots <- plot_spatial_results(regions[[region_name]], region_name)
    all_plots[[region_name]] <- plots
    
    # Generate summary statistics
    region_summary <- generate_spatial_summary(regions[[region_name]], region_name)
    summary_stats <- rbind(summary_stats, region_summary)
    
    # Analyze malignant subpopulations
    malignant_analysis <- analyze_malignant_subpopulations(regions[[region_name]], region_name)
    if (!is.null(malignant_analysis)) {
      malignant_analyses[[region_name]] <- malignant_analysis
    }
    
    # Save individual results
    qsave(rctd_results, 
          file = file.path(output_dir, paste0("RCTD_results_", region_name, ".qs")),
          nthreads = param$cores)
    
    write.csv(rctd_results$cell_groups,
              file = file.path(output_dir, paste0(region_name, "_cell_groups.csv")),
              row.names = FALSE)
    
    # Save plots
    for (plot_name in names(plots)) {
      ggsave(
        filename = file.path(output_dir, paste0(region_name, "_", plot_name, ".pdf")),
        plot = plots[[plot_name]],
        width = 10, height = 8
      )
    }
    
    results[[region_name]] <- list(
      xenium_obj = regions[[region_name]],
      rctd_results = rctd_results,
      plots = plots,
      summary = region_summary,
      malignant_analysis = malignant_analysis
    )
  }
  
  # Save comprehensive results
  write.csv(summary_stats, 
            file = file.path(output_dir, "region_summary_stats.csv"), 
            row.names = FALSE)
  
  qsave(summary_stats,
        file = file.path(output_dir, "region_summary_stats.qs"),
        nthreads = param$cores)
  
  # Save malignant analyses
  if (length(malignant_analyses) > 0) {
    saveRDS(malignant_analyses, file.path(output_dir, "malignant_analyses.rds"))
  }
  
  # Create combined summary plot
  if (nrow(summary_stats) > 1) {
    summary_plot <- create_summary_plots(summary_stats, malignant_analyses)
    ggsave(
      filename = file.path(output_dir, "combined_summary.pdf"),
      plot = summary_plot,
      width = 12, height = 8
    )
  }
  
  cat("\n=== Xenium Processing Workflow Completed ===\n")
  cat("Results saved to:", output_dir, "\n")
  
  return(list(
    regions = results,
    summary_stats = summary_stats,
    malignant_analyses = malignant_analyses,
    plots = all_plots
  ))
}

#' Create summary plots across all regions
#'
#' @param summary_stats Summary statistics data frame
#' @param malignant_analyses List of malignant analyses
#' @return Combined ggplot object
create_summary_plots <- function(summary_stats, malignant_analyses) {
  
  # Plot 1: Cell type diversity across regions
  p1 <- ggplot(summary_stats, aes(x = Region, y = Unique_Cell_Types)) +
    geom_col(fill = "steelblue") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Cell Type Diversity by Region", 
         y = "Number of Cell Types")
  
  # Plot 2: Doublet rates
  p2 <- ggplot(summary_stats, aes(x = Region, y = Doublet_Rate)) +
    geom_col(fill = "orange") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Doublet Rate by Region", 
         y = "Doublet Rate")
  
  # Plot 3: Malignant cell proportions (if available)
  if (length(malignant_analyses) > 0) {
    malignant_df <- map_dfr(malignant_analyses, function(x) {
      data.frame(
        Region = x$region,
        mTRM_enriched = x$mtrm_enriched,
        Total_Malignant = x$total_malignant_cells
      )
    })
    
    p3 <- ggplot(malignant_df, aes(x = Region, y = Total_Malignant, fill = mTRM_enriched)) +
      geom_col() +
      scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "darkred"),
                       name = "mTRM\nEnriched") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Malignant Cells by Region", 
           y = "Number of Malignant Cells")
    
    combined_plot <- p1 / p2 / p3
  } else {
    combined_plot <- p1 / p2
  }
  
  combined_plot <- combined_plot +
    plot_annotation(
      title = "Xenium Spatial Analysis Summary",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  return(combined_plot)
}

# Example usage (uncomment and modify paths as needed):
# 
# # Process Xenium data
# results <- process_xenium_workflow(
#   base_dir = "/path/to/xenium/output/",
#   reference_path = "/path/to/reference.qs",
#   output_dir = "Xenium_analysis_results"
# )

cat("Xenium processing script loaded successfully!\n")
cat("Use process_xenium_workflow() to run the complete analysis pipeline.\n")
