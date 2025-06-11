#!/usr/bin/env Rscript

#' SComatic Somatic Mutation Analysis Pipeline
#' 
#' This script processes SComatic pipeline output to analyze somatic mutations
#' in single-cell RNA-seq data, with focus on tumor mutational burden (TMB)
#' analysis in malignant vs non-malignant T-cell populations.
#'
#' Based on the CTCL single-cell memory manuscript analysis pipeline

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(ggsignif)
  library(MutationalPatterns)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(kableExtra)
  library(future)
  library(conflicted)
})

# Resolve conflicts
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::mutate)

# Set up parallel processing
param <- list(cores = 4, ram = 60, seed = 38)
plan("multicore", workers = param$cores)
set.seed(param$seed)
options(future.globals.maxSize = param$ram * 1024^3)

#' Define cell type mapping for consolidation
#'
#' @return Tibble with cell type mapping information
create_cell_type_mapping <- function() {
  
  cell_type_mapping <- tibble(
    original_type = c(
      # Malignant T cells (key focus for analysis)
      "malignant_Cytotoxic_T", "malignant_TCM1", "malignant_TCM2",
      
      # CD4 T cells (control - expected to have lower TMB)
      "CD4*CTL", "CD4_naive_T", "CD4_T", "CD4_TCM1", "CD4_TCM2", "CD4_TCM3", "CD4_TEM",
      
      # CD8 T cells (control - expected to have lower TMB)
      "CD8_naive_T", "CD8_TEM",
      
      # Other cell types
      "NK", "gdT", "Platelets"
    ),
    consolidated_type = c(
      # Malignant T cells
      rep("malignant_T", 3),
      
      # CD4 T cells
      rep("CD4_T", 7),
      
      # CD8 T cells
      rep("CD8_T", 2),
      
      # Other cell types
      "NK", "gdT", "Platelets"
    ),
    # Define expected malignancy status
    is_malignant = c(
      rep(TRUE, 3),    # malignant T cells
      rep(FALSE, 7),   # CD4 T cells
      rep(FALSE, 2),   # CD8 T cells
      FALSE,           # NK
      FALSE,           # gdT
      FALSE            # Platelets
    ),
    # Define cell type category for analyses
    cell_category = c(
      rep("Malignant", 3),  # malignant T cells
      rep("T cells", 7),    # CD4 T cells
      rep("T cells", 2),    # CD8 T cells
      "NK cells",           # NK
      "T cells",            # gdT
      "Other"               # Platelets
    )
  )
  
  return(cell_type_mapping)
}

#' Load and process callable sites data
#'
#' @param callable_sites_file Path to callable sites TSV file
#' @param cell_type_mapping Cell type mapping tibble
#' @param min_coverage Minimum coverage threshold
#' @return Processed callable sites data
load_callable_sites <- function(callable_sites_file, cell_type_mapping, min_coverage = 10) {
  
  cat("Loading callable sites from:", callable_sites_file, "\n")
  
  # Check if file exists
  if (!file.exists(callable_sites_file)) {
    warning("Callable sites file not found, creating dummy data")
    return(create_dummy_callable_sites(cell_type_mapping))
  }
  
  # Read callable sites data
  all_sites_raw <- read_tsv(callable_sites_file, show_col_types = FALSE) %>%
    # Filter for all relevant cell types
    filter(Cell_types %in% cell_type_mapping$original_type) %>%
    # Filter by minimum coverage threshold
    filter(Cov >= min_coverage) %>%
    # Add the consolidated cell type category
    left_join(cell_type_mapping, by = c("Cell_types" = "original_type"))
  
  # Create summarized data 
  all_sites_summary <- all_sites_raw %>%
    # Group by consolidated category and sample
    group_by(consolidated_type, sample) %>% 
    # Calculate sum of DP and NC for each group
    summarize(
      callable_sites = sum(DP, na.rm = TRUE),
      total_cells = sum(NC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Remove any potential duplicates
    distinct()
  
  # Create final callable sites data frame
  all_sites <- data.frame(
    Cell_types = all_sites_summary$consolidated_type,
    sample = all_sites_summary$sample,
    callable_sites = all_sites_summary$callable_sites,
    total_cells = all_sites_summary$total_cells
  )
  
  return(all_sites)
}

#' Create dummy callable sites data when file is missing
#'
#' @param cell_type_mapping Cell type mapping tibble
#' @return Dummy callable sites data frame
create_dummy_callable_sites <- function(cell_type_mapping) {
  
  # Create dummy data with all combinations
  all_sites <- expand.grid(
    Cell_types = unique(cell_type_mapping$consolidated_type),
    sample = unique(c(paste0("CTCL_", 1:9), paste0("P", c(1,2,4,5,6,7,8)))),
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(
      callable_sites = 1000000,  # Default value
      total_cells = 100          # Default value
    )
  
  return(all_sites)
}

#' Process variant calling results from SComatic step2 output
#'
#' @param file_path Path to step2 TSV file
#' @return Processed variants data frame
process_variant_calling <- function(file_path) {
  
  cat("Processing variant file:", basename(file_path), "\n")
  
  # Read file lines to find header
  lines <- readLines(file_path)
  header_line <- which(grepl("^#CHROM", lines))
  
  if (length(header_line) == 0) {
    warning("Could not find column header line in file:", file_path)
    return(NULL)
  }
  
  # Extract column names and read data
  col_names <- strsplit(sub("^#", "", lines[header_line]), "\t")[[1]]
  variants <- read_tsv(file_path, skip = header_line, col_names = col_names, show_col_types = FALSE)
  
  return(variants)
}

#' Process variants by consolidated cell type categories
#'
#' @param variants_df Variants data frame
#' @param cell_type_mapping Cell type mapping tibble
#' @return Consolidated mutation counts by cell type
process_variants_by_consolidated_category <- function(variants_df, cell_type_mapping) {
  
  if (is.null(variants_df) || nrow(variants_df) == 0) {
    # Return zeros for all cell types
    return(tibble(
      Cell_types = unique(cell_type_mapping$consolidated_type),
      mutations = 0
    ))
  }
  
  # Filter for PASS variants
  variants_pass_filter <- variants_df %>%
    filter(grepl("PASS", FILTER))
  
  if (nrow(variants_pass_filter) == 0) {
    # Return zeros for all cell types
    return(tibble(
      Cell_types = unique(cell_type_mapping$consolidated_type),
      mutations = 0
    ))
  }
  
  # Count variants by their primary cell type assignment
  cell_type_counts <- variants_pass_filter %>%
    count(Cell_types, name = "mutations") %>%
    rename(original_type = Cell_types)
  
  # Join with mapping to get consolidated types
  counts_with_mapping <- cell_type_counts %>%
    left_join(cell_type_mapping, by = "original_type") %>%
    # Handle any cell types not in our mapping
    mutate(consolidated_type = ifelse(is.na(consolidated_type), 
                                    original_type, 
                                    consolidated_type))
  
  # Sum by consolidated type
  consolidated_counts <- counts_with_mapping %>%
    group_by(consolidated_type) %>%
    summarise(mutations = sum(mutations, na.rm = TRUE), .groups = "drop") %>%
    rename(Cell_types = consolidated_type)
  
  # Ensure all expected cell types are present
  all_expected_types <- unique(cell_type_mapping$consolidated_type)
  final_counts <- tibble(Cell_types = all_expected_types) %>%
    left_join(consolidated_counts, by = "Cell_types") %>%
    mutate(mutations = replace_na(mutations, 0))
  
  return(final_counts)
}

#' Process all SComatic variant files
#'
#' @param work_dir Directory containing SComatic results
#' @param sample_mapping Sample mapping data frame
#' @param cell_type_mapping Cell type mapping tibble
#' @return Combined results data frame
process_all_variant_files <- function(work_dir, sample_mapping, cell_type_mapping) {
  
  cat("Searching for variant files in:", work_dir, "\n")
  
  # Find step2 files
  step2_files <- list.files(work_dir, pattern = "*.calling.step2.tsv", 
                           recursive = TRUE, full.names = TRUE)
  
  # Filter for relevant samples
  sample_pattern <- paste0("(", paste(sample_mapping$simple_name, collapse = "|"), ")\\.calling\\.step2\\.tsv")
  final_step2_files <- step2_files[grep(sample_pattern, step2_files)]
  
  cat("Found", length(final_step2_files), "relevant variant files\n")
  
  # Process all samples
  results <- list()
  
  for (file in final_step2_files) {
    simple_name <- sub("\\.calling\\.step2\\.tsv.*", "", basename(file))
    cat("Processing sample:", simple_name, "\n")
    
    variants <- process_variant_calling(file)
    if (!is.null(variants)) {
      cluster_results <- process_variants_by_consolidated_category(variants, cell_type_mapping)
      if (!is.null(cluster_results)) {
        results[[simple_name]] <- cluster_results %>%
          mutate(simple_name = simple_name)
      }
    }
  }
  
  # Combine results
  if (length(results) > 0) {
    all_variants <- bind_rows(results)
  } else {
    warning("No variant results found")
    all_variants <- tibble()
  }
  
  return(all_variants)
}

#' Calculate tumor mutational burden (TMB) and related metrics
#'
#' @param variants_df Variants data frame
#' @param callable_sites_df Callable sites data frame
#' @param sample_mapping Sample mapping data frame
#' @param cell_type_mapping Cell type mapping tibble
#' @param min_callable_sites Minimum callable sites threshold
#' @return Data frame with TMB calculations
calculate_tmb_metrics <- function(variants_df, callable_sites_df, sample_mapping, 
                                 cell_type_mapping, min_callable_sites = 100000) {
  
  cat("Calculating TMB metrics...\n")
  
  if (nrow(variants_df) == 0) {
    warning("No variant data available for TMB calculation")
    return(create_empty_tmb_results(cell_type_mapping))
  }
  
  # Combine variants with sample type and callable sites info
  final_results <- variants_df %>%
    left_join(sample_mapping, by = c("simple_name" = "simple_name")) %>%
    rename(sample = simple_name) %>%
    left_join(callable_sites_df, by = c("Cell_types", "sample")) %>%
    # Filter out data with insufficient callable sites
    filter(!is.na(callable_sites) & callable_sites >= min_callable_sites) %>%
    # Filter out platelets with small number of callable sites
    filter(!(Cell_types == "Platelets" & callable_sites < 50000)) %>%
    # Join cell type metadata
    left_join(
      cell_type_mapping %>% 
        select(consolidated_type, is_malignant, cell_category) %>%
        rename(Cell_types = consolidated_type) %>%
        distinct(),
      by = "Cell_types"
    ) %>%
    distinct()
  
  if (nrow(final_results) == 0) {
    warning("No data remaining after filtering")
    return(create_empty_tmb_results(cell_type_mapping))
  }
  
  # Calculate TMB and related metrics
  final_results <- final_results %>%
    mutate(
      # Standard TMB calculation (mutations per megabase)
      tmb = (mutations * 1e6) / callable_sites,
      
      # Add confidence score based on mutation count
      confidence_score = case_when(
        mutations >= 10 ~ "High",
        mutations >= 3 ~ "Medium", 
        TRUE ~ "Low"
      ),
      
      # Flag outliers
      is_outlier = ifelse(
        nrow(.) > 1 & !all(is.na(tmb)) & sd(tmb, na.rm = TRUE) != 0, 
        tmb > mean(tmb, na.rm = TRUE) + 3 * sd(tmb, na.rm = TRUE), 
        FALSE
      )
    )
  
  # Calculate sample-specific normalization metrics
  final_results <- final_results %>%
    group_by(sample) %>%
    mutate(
      # Within-sample normalization
      sample_norm_tmb = tmb / max(tmb[Cell_types != "Platelets"], na.rm = TRUE),
      
      # Calculate non-malignant baseline in this sample
      non_malignant_tmb = mean(tmb[!is_malignant & Cell_types != "Platelets"], na.rm = TRUE),
      
      # Calculate malignancy enrichment score 
      malignancy_enrichment = ifelse(is_malignant, 
                                    tmb / non_malignant_tmb, 
                                    NA_real_)
    ) %>%
    ungroup()
  
  # Add cell type level metrics
  final_results <- final_results %>%
    group_by(Cell_types, type) %>% 
    mutate(
      # Calculate mean TMB per cell type and sample type
      cell_type_mean_tmb = mean(tmb, na.rm = TRUE),
      
      # Calculate relative TMB compared to cell type average
      cell_type_relative_tmb = tmb / cell_type_mean_tmb
    ) %>%
    ungroup()
  
  return(final_results)
}

#' Create empty TMB results data frame
#'
#' @param cell_type_mapping Cell type mapping tibble
#' @return Empty TMB results data frame with correct structure
create_empty_tmb_results <- function(cell_type_mapping) {
  
  empty_results <- tibble(
    Cell_types = unique(cell_type_mapping$consolidated_type),
    mutations = 0,
    sample = "dummy",
    type = "dummy",
    callable_sites = 1000000,
    total_cells = 100,
    tmb = 0,
    confidence_score = "Low",
    is_outlier = FALSE,
    is_malignant = FALSE,
    cell_category = "Other"
  )
  
  return(empty_results)
}

#' Create TMB boxplot with statistical significance testing
#'
#' @param data TMB results data frame
#' @param test_method Statistical test method
#' @param log_scale Whether to use log scale
#' @param min_tmb Minimum TMB value to include
#' @return ggplot object
plot_tmb_with_significance <- function(data, test_method = "wilcox.test",
                                     log_scale = FALSE, min_tmb = 0) {
  
  # Filter data if needed
  plot_data <- data
  if (min_tmb > 0) {
    plot_data <- data %>% filter(tmb > min_tmb)
  }
  
  # Order cell types
  plot_data <- plot_data %>%
    mutate(Cell_types = factor(Cell_types, 
                              levels = c("malignant_T", "CD4_T", "CD8_T", "NK", "gdT", "Platelets"))) %>%
    filter(Cell_types %in% c("malignant_T", "CD4_T", "CD8_T", "NK", "gdT", "Platelets"))
  
  # Define comparisons
  comparisons_list <- list(
    c("malignant_T", "CD4_T"),
    c("malignant_T", "CD8_T")
  )
  
  # Create the base plot
  p <- ggplot(plot_data, aes(x = Cell_types, y = tmb)) +
    geom_boxplot(aes(fill = is_malignant), 
                outlier.shape = NA, 
                alpha = 0.7,
                width = 0.6) +
    geom_jitter(aes(color = type, shape = sample), 
                position = position_jitter(width = 0.2), 
                size = 2.5, 
                alpha = 0.8) +
    scale_fill_manual(values = c("FALSE" = "#4CAF50", "TRUE" = "#F44336"),
                     name = "Cell Status",
                     labels = c("FALSE" = "Non-malignant", "TRUE" = "Malignant")) +
    scale_color_manual(values = c("Skin" = "#FF9800", "PBMC" = "#2196F3"),
                      name = "Sample Type") +
    scale_shape_manual(values = c(16, 17, 15, 3, 4, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20),
                      name = "Sample") +
    labs(
      title = "Tumor Mutational Burden (TMB) Comparison",
      x = "Cell Type",
      y = "TMB (Mutations per Mb)",
      caption = paste("Statistical test:", test_method)
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 13, face = "bold"),
      plot.title = element_text(size = 15, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  # Apply log scale if requested
  if (log_scale) {
    p <- p + scale_y_log10(labels = scales::comma)
  }
  
  # Add significance testing
  p <- p + 
    geom_signif(
      comparisons = comparisons_list,
      test = test_method,
      map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
      textsize = 4,
      step_increase = 0.1,
      tip_length = 0.01,
      color = "black"
    )
  
  return(p)
}

#' Create comprehensive TMB comparison plots
#'
#' @param data TMB results data frame
#' @return Combined ggplot object
plot_tmb_comparisons <- function(data) {
  
  # Plot 1: TMB by cell type
  p1 <- ggplot(data, aes(x = Cell_types, y = tmb, fill = is_malignant)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(), alpha = 0.6) +
    scale_fill_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "TMB by Cell Type", x = "Cell Type", y = "TMB")
  
  # Plot 2: TMB by sample, colored by type
  p2 <- ggplot(data, aes(x = sample, y = tmb, fill = Cell_types)) +
    geom_col(position = "dodge") +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "TMB by Sample", x = "Sample", y = "TMB")
  
  # Plot 3: Malignant vs Non-malignant TMB ratio
  p3 <- data %>%
    filter(Cell_types != "Platelets") %>%
    group_by(sample, type) %>%
    summarise(
      malignant_tmb = mean(tmb[is_malignant], na.rm = TRUE),
      non_malignant_tmb = mean(tmb[!is_malignant], na.rm = TRUE),
      ratio = malignant_tmb / non_malignant_tmb,
      .groups = "drop"
    ) %>%
    ggplot(aes(x = reorder(sample, ratio), y = ratio, fill = type)) +
    geom_col() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    coord_flip() +
    scale_fill_manual(values = c("Skin" = "lightblue", "PBMC" = "salmon")) +
    theme_bw() +
    labs(title = "Malignant to Non-Malignant TMB Ratio", 
         subtitle = "Values > 1 indicate higher TMB in malignant cells",
         x = "Sample", y = "Ratio")
  
  # Combine plots
  combined_plot <- p1 / p2 / p3 +
    plot_layout(heights = c(2, 2, 3)) +
    plot_annotation(
      title = "Tumor Mutational Burden (TMB) Analysis",
      subtitle = "Comparison across cell types and samples",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  return(combined_plot)
}

#' Perform pairwise statistical tests between cell types
#'
#' @param data TMB results data frame
#' @param test_method Statistical test method
#' @return Data frame with test results
perform_pairwise_tests <- function(data, test_method = "wilcox.test") {
  
  # Filter for relevant cell types
  test_data <- data %>%
    filter(Cell_types %in% c("malignant_T", "CD4_T", "CD8_T", "NK", "gdT", "Platelets")) %>%
    filter(!is.na(tmb))
  
  # Define comparisons
  comparisons <- list(
    c("malignant_T", "CD4_T"),
    c("malignant_T", "CD8_T")
  )
  
  # Perform tests
  test_results <- map_dfr(comparisons, function(comp) {
    group1_data <- test_data %>% filter(Cell_types == comp[1]) %>% pull(tmb)
    group2_data <- test_data %>% filter(Cell_types == comp[2]) %>% pull(tmb)
    
    if (length(group1_data) > 0 && length(group2_data) > 0) {
      if (test_method == "wilcox.test") {
        test_result <- wilcox.test(group1_data, group2_data)
      } else {
        test_result <- t.test(group1_data, group2_data)
      }
      
      tibble(
        comparison = paste(comp[1], "vs", comp[2]),
        group1 = comp[1],
        group2 = comp[2],
        p_value = test_result$p.value,
        p_value_formatted = ifelse(test_result$p.value < 0.001, "< 0.001", 
                                  format(round(test_result$p.value, 4), scientific = FALSE)),
        significance = case_when(
          test_result$p.value < 0.001 ~ "***",
          test_result$p.value < 0.01 ~ "**", 
          test_result$p.value < 0.05 ~ "*",
          TRUE ~ "ns"
        ),
        test_method = test_method,
        n_group1 = length(group1_data),
        n_group2 = length(group2_data)
      )
    } else {
      tibble(
        comparison = paste(comp[1], "vs", comp[2]),
        group1 = comp[1],
        group2 = comp[2],
        p_value = NA_real_,
        p_value_formatted = "Insufficient data",
        significance = "ns",
        test_method = test_method,
        n_group1 = length(group1_data),
        n_group2 = length(group2_data)
      )
    }
  })
  
  return(test_results)
}

#' Main SComatic analysis workflow
#'
#' @param work_dir Directory containing SComatic results
#' @param callable_sites_file Path to callable sites file
#' @param output_dir Output directory for results
#' @param min_coverage Minimum coverage threshold
#' @param min_callable_sites Minimum callable sites threshold
#' @return List containing analysis results
process_scomatic_workflow <- function(work_dir, callable_sites_file,
                                     output_dir = "SComatic_analysis_output",
                                     min_coverage = 10,
                                     min_callable_sites = 100000) {
  
  cat("=== Starting SComatic Analysis Workflow ===\n")
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create cell type mapping
  cell_type_mapping <- create_cell_type_mapping()
  
  # Create sample mapping
  sample_mapping <- data.frame(
    simple_name = c(paste0("CTCL_", 1:9), paste0("P", c(1,2,4,5,6,7,8))),
    type = c(rep("Skin", 9), rep("PBMC", 7))
  )
  
  # Load callable sites
  callable_sites <- load_callable_sites(callable_sites_file, cell_type_mapping, min_coverage)
  
  # Process variant files
  all_variants <- process_all_variant_files(work_dir, sample_mapping, cell_type_mapping)
  
  # Calculate TMB metrics
  tmb_results <- calculate_tmb_metrics(all_variants, callable_sites, sample_mapping, 
                                      cell_type_mapping, min_callable_sites)
  
  # Generate plots
  plots <- list()
  
  if (nrow(tmb_results) > 0 && sum(!is.na(tmb_results$tmb)) > 0) {
    
    # TMB significance plots
    plots$tmb_significance <- plot_tmb_with_significance(
      tmb_results, 
      test_method = "wilcox.test",
      log_scale = FALSE,
      min_tmb = 0
    )
    
    plots$tmb_significance_log <- plot_tmb_with_significance(
      tmb_results,
      test_method = "wilcox.test", 
      log_scale = TRUE,
      min_tmb = 0.1
    )
    
    # Comprehensive comparison plots
    plots$tmb_comparisons <- plot_tmb_comparisons(tmb_results)
    
    # Perform statistical tests
    test_results <- perform_pairwise_tests(tmb_results)
    
  } else {
    cat("No TMB data available for plotting.\n")
    test_results <- tibble()
  }
  
  # Save results
  write_tsv(tmb_results, file.path(output_dir, "tmb_results.tsv"))
  write_tsv(test_results, file.path(output_dir, "statistical_tests.tsv"))
  saveRDS(tmb_results, file.path(output_dir, "tmb_results.rds"))
  
  # Save plots
  for (plot_name in names(plots)) {
    ggsave(
      filename = file.path(output_dir, paste0(plot_name, ".pdf")),
      plot = plots[[plot_name]],
      width = 12, height = 10
    )
  }
  
  cat("=== SComatic Analysis Workflow Completed ===\n")
  cat("Results saved to:", output_dir, "\n")
  
  return(list(
    tmb_results = tmb_results,
    test_results = test_results,
    plots = plots,
    cell_type_mapping = cell_type_mapping
  ))
}

# Example usage (uncomment and modify paths as needed):
# 
# # Run SComatic analysis
# results <- process_scomatic_workflow(
#   work_dir = "/path/to/SComatic/results/final_output/",
#   callable_sites_file = "/path/to/callable_sites/all_callable_sites.tsv",
#   output_dir = "SComatic_analysis_results"
# )

cat("SComatic analysis script loaded successfully!\n")
cat("Use process_scomatic_workflow() to run the complete analysis pipeline.\n")
