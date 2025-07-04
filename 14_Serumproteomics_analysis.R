#!/usr/bin/env Rscript

#' CTCL Serum Proteomic Analysis using Olink Immune Response Panel
#'
#' This script processes serum proteomic data for CTCL patients using the Olink platform.
#' It performs QC, mixed-model-based differential expression analysis, volcano plotting,
#' and identifies significantly deregulated proteins between Sézary syndrome and healthy controls.
#'
#' Dependencies: openxlsx, tidyverse, OlinkAnalyze, EnhancedVolcano

# Load required libraries
suppressPackageStartupMessages({
  library(openxlsx)
  library(tidyverse)
  library(OlinkAnalyze)
  library(EnhancedVolcano)
  library(ggrepel)
})

#' Load and clean NPX data
load_clean_data <- function(npx_path, meta_path) {
  npx <- read_NPX(npx_path)
  meta <- read.xlsx(meta_path)
  data <- merge(npx, meta, by = "SampleID") %>%
    filter(QC_Warning == "Pass") %>%
    filter(Panel == "Olink Immune Response") %>%
    mutate(Disease = recode(Disease, "SS" = "Sézary", "HD" = "Healthy"),
           Disease = factor(Disease, levels = c("Sézary","Healthy")))
  return(data)
}

#' Run linear mixed-effect model
run_lmer_analysis <- function(df) {
  res <- olink_lmer(
    df = df,
    variable = "Disease",
    random = "Patient",
    covariates = c("Age", "Sex")
  )
  return(res)
}

#' Run post-hoc analysis
run_posthoc <- function(df, protein_list) {
  olink_lmer_posthoc(
    df = df,
    olinkid_list = protein_list,
    variable = "Disease",
    random = "Patient",
    covariates = c("Age", "Sex"),
    effect = "Disease"
  )
}

#' Volcano plot with significance thresholds
plot_volcano <- function(res) {
  res$Color <- "NS or estimate < 0.5"
  res$Color[res$Adjusted_pval < 0.05] <- "P < 0.05"
  res$Color[res$estimate > 0.5] <- "estimate > 0.5"
  res$Color[res$estimate < -0.5] <- "estimate < -0.5"
  res$Color[abs(res$Adjusted_pval) > 0.05] <- "NS or estimate < 0.5"
  res$Color <- factor(res$Color, levels = c("estimate > 0.5", "estimate < -0.5", "P < 0.05", "NS or estimate < 0.5"))

  ggplot(res, aes(x = estimate, y = -log10(Adjusted_pval), color = Color, label = Assay)) +
    geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    geom_point() +
    scale_color_manual(values = c("estimate > 0.5" = "dodgerblue",
                                  "P < 0.05" = "#99CC99",
                                  "estimate < -0.5" = "firebrick",
                                  "NS or estimate < 0.5" = "gray"),
                       guide = guide_legend(override.aes = list(size = 4))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    geom_text_repel(
      data = res %>% filter((estimate > 0.5 | estimate < -0.5) & Adjusted_pval < 0.05),
      aes(label = Assay), size = 3, box.padding = 0.5, max.overlaps = Inf
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "right")
}

# Main analysis pipeline
main <- function() {
  cat("Loading and cleaning data...\n")
  data_path <- "../raw/Roggo_Inflammation_Immune Response_NPX_short.xlsx"
  meta_path <- "../data/Serum_Lymphoma_clean_nokey.xlsx"
  data <- load_clean_data(data_path, meta_path)

  cat("Running linear mixed model...\n")
  lmer_results <- run_lmer_analysis(data)

  significant_proteins <- lmer_results %>%
    filter(Threshold == "Significant") %>%
    pull(OlinkID)

  if (length(significant_proteins) == 0) {
    stop("No significant proteins found.")
  }

  cat("Running post-hoc analysis...\n")
  posthoc_res <- run_posthoc(data, significant_proteins)
  rownames(posthoc_res) <- posthoc_res$Assay

  cat("Generating volcano plot...\n")
  print(plot_volcano(posthoc_res))

  cat("Generating per-protein plots...\n")
  for (protein in significant_proteins) {
    plots <- olink_lmer_plot(
      df = data,
      olinkid_list = protein,
      variable = "Disease",
      x_axis_variable = "Disease",
      random = "Patient",
      covariates = c("Age", "Sex"),
      number_of_proteins_per_plot = 1
   )
    print(plots[[1]])
}

}

# Run the workflow
main()
