#!/usr/bin/env Rscript

#' CTCL Malignant T-cell Developmental Potential Scoring (CytoTRACE)
#'
#' This script applies CytoTRACE2 to malignant T cells from CTCL blood and skin
#' samples to infer differentiation potential.
#'
#' Dependencies: Seurat, CytoTRACE2, ggplot2, qs

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(CytoTRACE2)
  library(ggplot2)
  library(qs)
})

#' Load and prepare Seurat object
#' @param path Path to the Seurat .rds object
#' @return Prepared Seurat object
load_data <- function(path) {
  seu <- readRDS(path)
  seu <- SetIdent(seu, value = seu@meta.data$TCM)
  DefaultAssay(seu) <- "RNA"
  return(seu)
}

#' Run CytoTRACE2 analysis
#' @param seu Seurat object
#' @param ncores Number of cores to use
#' @return Seurat object with CytoTRACE2 scores added
run_cytotrace <- function(seu, ncores = 20) {
  cytotrace2_result <- cytotrace2(
    seu,
    is_seurat = TRUE,
    ncores = ncores,
    slot_type = "counts",
    species = "human"
  )
  return(cytotrace2_result)
}

#' Plot CytoTRACE score by TCM subtype and tissue
#' @param seu Seurat object with CytoTRACE2_Score
plot_cytotrace_violin <- function(seu) {
  VlnPlot(
    seu,
    features = "CytoTRACE2_Score",
    pt.size = 0,
    group.by = "TCM",
    split.by = "tissue",
    cols = c("blood" = "#DC143C", "skin" = "#8B4513")
  ) +
    ggtitle("Developmental Potential by Phenotype") +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.2),
      limits = c(0, 1),
      sec.axis = sec_axis(
        trans = ~.,
        breaks = seq(0, 1, length.out = 13),
        labels = c(
          "", "Differentiated", "", "Unipotent", "", "Oligopotent", "",
          "Multipotent", "", "Pluripotent", "", "Totipotent", ""
        )
      )
    ) +
    labs(x = "Phenotype", y = "Potency Score") +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
      axis.ticks.y.right = element_line(color = rep(c("black", NA), 7)),
      axis.ticks.length.y.right = unit(0.3, "cm"),
      aspect.ratio = 0.8
    )
}
