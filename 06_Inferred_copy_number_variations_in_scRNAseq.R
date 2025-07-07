#!/usr/bin/env Rscript

#' CTCL scRNA-seq CNV Inference via inferCNV
#'
#' This script performs copy number variation (CNV) inference on CTCL single-cell RNA-seq data
#' using the `infercnv` package. It processes both blood and skin T-cell samples, generates
#' annotation and count matrices, runs the inferCNV model, and visualizes CNV heatmaps.
#'
#' Dependencies: Seurat, infercnv, ComplexHeatmap, circlize, data.table, qs, SCpubr, dplyr, paletteer

# Load packages
suppressPackageStartupMessages({
  library(Seurat)
  library(infercnv)
  library(SCpubr)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(paletteer)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(dittoSeq)
  library(qs)
})

# Setup
set.seed(42)
options(scipen = 100)

#' Create input files for inferCNV
prepare_infercnv_inputs <- function(seu, sample_type = "blood", output_dir = ".") {
  cat("Preparing input matrices for:", sample_type, "\n")
  DefaultAssay(seu) <- "RNA"
  Idents(seu) <- seu$Tcells
  
  count_path <- file.path(output_dir, paste0("raw_count_matrix_CTCL_T_", sample_type, ".txt"))
  anno_path <- file.path(output_dir, paste0(ucfirst(sample_type), "_Tcells_anno.txt"))
  
  # Write raw counts
  write.table(GetAssayData(seu, slot = "counts"), count_path, sep = "\t", col.names = TRUE, row.names = TRUE)
  
  # Create annotation
  meta <- seu@meta.data
  meta$Barcode <- rownames(meta)
  anno <- meta[, c("Barcode", "Tcells")]
  rownames(anno) <- anno$Barcode
  write.table(anno["Tcells"], anno_path, sep = "\t", col.names = FALSE, row.names = TRUE)
  
  return(list(count_matrix = count_path, annotation = anno_path))
}

#' Run inferCNV
run_infercnv <- function(count_matrix, annotation, gene_order_file, ref_group = "NK", out_dir = "InferCNV_out") {
  cat("Running inferCNV for:", out_dir, "\n")
  infer_obj <- CreateInfercnvObject(
    raw_counts_matrix = count_matrix,
    annotations_file = annotation,
    gene_order_file = gene_order_file,
    ref_group_names = ref_group
  )
  
  infercnv::run(
    infer_obj,
    cutoff = 0.1,
    out_dir = out_dir,
    HMM = TRUE,
    num_threads = 50,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    sd_amplifier = 3,
    noise_logistic = TRUE,
    png_res = 300,
    BayesMaxPNormal = 0.2,
    tumor_subcluster_partition_method = "leiden",
    num_ref_groups = 2,
    z_score_filter = 0.5
  )
}

#' Plot CNV heatmap from inferCNV object
plot_infercnv_heatmap <- function(hmm_path, seu, sample_type = "blood") {
  cat("Plotting CNV heatmap for:", sample_type, "\n")
  
  hmm_obj <- readRDS(hmm_path)
  meta <- seu@meta.data
  meta$barcode <- rownames(meta)
  
  if (sample_type == "blood") {
    # Subsample per patient for clarity
    sample_cells <- function(df) {
      mal <- df %>% filter(Tcells == "Malignant TCM") %>% slice_sample(n = min(200, n()))
      rest <- df %>% filter(Tcells != "Malignant TCM") %>% slice_sample(n = min(100, n()))
      bind_rows(mal, rest)
    }
    
    clusters <- meta %>%
      filter(patient %in% paste0("blood ", 1:8)) %>%
      select(barcode, Tcells, patient) %>%
      group_by(patient) %>%
      group_split() %>%
      lapply(sample_cells) %>%
      bind_rows() %>%
      arrange(Tcells, patient)
    
    celltype_colors <- c(
      "Malignant TCM" = "#000080", "CD4 naive T" = "#1CBE4FFF", "CD4 TCM" = "#90AD1CFF", "CD4 TEM" = "#1C8356FF",
      "Treg" = "#006600FF", "CTL" = "#664466FF", "CD8 naive T" = "#CC99CCFF", "CD8 TEM" = "#FF9966FF",
      "gdT" = "#BB6622FF", "Prolif T" = "#FFCC99FF", "NK TOX+" = "#CC6666FF", "NK" = "#882211FF"
    )
  } else {
    clusters <- meta %>% select(barcode = rownames(meta), Tcells, patient) %>% arrange(Tcells, patient)
    celltype_colors <- c(
      "Malignant TCM" = "#000080", "Malignant IS" = "#0011EEFF", "CD4 TRM" = "#66CC33FF",
      "CD4 TEM" = "#BBAA55FF", "Treg" = "#006600FF", "CD8 T" = "#DD6644FF", "NK" = "#882211FF"
    )
  }
  
  mat <- t(hmm_obj@expr.data[, clusters$barcode])
  gene_order <- hmm_obj@gene_order
  range_vals <- range(mat)
  middle <- as.double(names(which.max(table(mat[, 1]))))
  
  color_fn <- if (range_vals[1] != 1) {
    circlize::colorRamp2(c(range_vals[1], middle, range_vals[2]), c("#00428c", "white", "#9e002a"))
  } else {
    circlize::colorRamp2(range_vals, c("white", "#9e002a"))
  }
  
  chrom_boundaries <- cumsum(table(gene_order$chr))
  chrom_labels <- rep("", ncol(mat))
  chrom_labels[chrom_boundaries] <- unique(gene_order$chr)
  
  col_anno <- HeatmapAnnotation(
    Chromosome = anno_text(chrom_labels, rot = 90, just = "right", gp = gpar(fontsize = 9)),
    annotation_height = unit(1.5, "cm")
  )
  
  row_anno <- rowAnnotation(
    Celltype = clusters$Tcells,
    Patient = clusters$patient,
    col = list(Celltype = celltype_colors, Patient = setNames(brewer.pal(length(unique(clusters$patient)), "Spectral"), unique(clusters$patient)))
  )
  
  Heatmap(
    mat,
    border_gp = gpar(col = "#a6a6a6"),
    use_raster = TRUE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    col = color_fn,
    column_labels = NULL,
    right_annotation = row_anno,
    bottom_annotation = col_anno
  )
}

# Capitalize helper
ucfirst <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))

# Main pipeline
main <- function() {
  gene_order_file <- "/hg38_gencode_v27.txt"
  
  # --- Blood ---
  cat("\n>>> Processing BLOOD T-cells <<<\n")
  blood_seu <- readRDS("/BLOOD_Tcells.rds")
  Idents(blood_seu) <- blood_seu$Tcells
  blood_inputs <- prepare_infercnv_inputs(blood_seu, "blood")
  run_infercnv(blood_inputs$count_matrix, blood_inputs$annotation, gene_order_file, "NK", "/InferCNV_blood_T_NKref")
  blood_hmm_path <- list.files("/InferCNV_blood_T_NKref", pattern = "20_HMM_pred.*infercnv_obj", full.names = TRUE)
  plot_infercnv_heatmap(blood_hmm_path, blood_seu, "blood")
  
  # --- Skin ---
  cat("\n>>> Processing SKIN T-cells <<<\n")
  skin_seu <- readRDS("/SKIN_Tcells.rds")
  Idents(skin_seu) <- skin_seu$Tcells
  skin_inputs <- prepare_infercnv_inputs(skin_seu, "skin")
  run_infercnv(skin_inputs$count_matrix, skin_inputs$annotation, gene_order_file, "NK", "/InferCNV_skin_T_NKref")
  skin_hmm_path <- list.files("/InferCNV_skin_T_NKref", pattern = "20_HMM_pred.*infercnv_obj", full.names = TRUE)
  plot_infercnv_heatmap(skin_hmm_path, skin_seu, "skin")
}

# Run pipeline
main()


