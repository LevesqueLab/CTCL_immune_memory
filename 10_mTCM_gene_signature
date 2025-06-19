#!/usr/bin/env Rscript

#' CTCL mTCM Malignancy Signature Enrichment (KEGG)
#'
#' This script identifies differentially upregulated genes in mTCM cells
#' across CTCL patients and performs KEGG over-representation analysis
#' using enrichR.
#'
#' Dependencies: Seurat, dplyr, enrichR, ggplot2, EnhancedVolcano, qs

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(enrichR)
  library(EnhancedVolcano)
  library(qs)
})

#' Load and prepare Seurat object
load_data <- function(path) {
  seu <- qread(path)
  seu <- subset(seu, subset = TCM %in% c("mTCM", "mTRM"))
  seu <- SetIdent(seu, value = seu$TCM)
  DefaultAssay(seu) <- "RNA"
  return(seu)
}

#' Perform per-patient DE analysis between mTCM and mTRM
run_patient_deg <- function(seu, patients) {
  final_result <- data.frame()
  
  for (pat in patients) {
    seu_pat <- subset(seu, subset = patient == pat)
    if (all(table(seu_pat$TCM) >= 3)) {
      res <- FindMarkers(seu_pat, ident.1 = "mTCM", ident.2 = "mTRM")
      res$gene <- rownames(res)
      res$patient <- pat
      res$regulation <- "non significant"
      res$regulation[res$avg_log2FC > 1.5 & res$p_val_adj < 0.001] <- "relevant up"
      final_result <- bind_rows(final_result, res)
    }
  }
  
  return(final_result)
}

#' Extract genes upregulated in at least 8 patients
get_consensus_genes <- function(deg_df, min_patients = 8) {
  deg_df %>%
    filter(regulation == "relevant up") %>%
    group_by(gene) %>%
    summarise(frequency = n(), .groups = "drop") %>%
    filter(frequency >= min_patients) %>%
    pull(gene)
}

#' Run enrichR KEGG pathway enrichment
run_kegg_enrichment <- function(genes) {
  kegg_db <- "KEGG_2021_Human"
  enriched <- enrichr(genes, databases = kegg_db)
  df <- enriched[[kegg_db]]
  df$Gene_Count <- as.numeric(sapply(strsplit(df$Overlap, "/"), `[`, 1))
  df$Gene_Ratio <- df$Gene_Count / as.numeric(sapply(strsplit(df$Overlap, "/"), `[`, 2))
  df$LCS <- log2(df$Combined.Score)
  return(df)
}

#' Plot top KEGG enriched pathways
plot_kegg_results <- function(df, top_n = 10) {
  df <- df[1:top_n, ]
  ggplot(df, aes(x = LCS, y = reorder(Term, LCS), size = Gene_Ratio, color = Adjusted.P.value)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = "mTCM KEGG Pathway Enrichment",
         x = "log2(Combined Score)",
         y = "KEGG Pathway",
         size = "Gene Ratio",
         color = "Adjusted P-value") +
    theme_minimal()
}

# Main pipeline
main <- function() {
  # Load data
  seu <- load_data("/Malignant_T_blood_skin.qs")
  
  # Patients to include (minimum 3 cells per group)
  patients <- c("blood 1", "blood 2", "blood 4", "blood 5", "blood 6", "blood 7", "blood 8",
                "skin 1", "skin 3", "skin 4", "skin 7", "skin 8")

  # Differential expression
  cat("Running DE analysis per patient...\n")
  deg <- run_patient_deg(seu, patients)
  
  # Consensus upregulated genes
  cat("Identifying consistently upregulated genes in mTCM...\n")
  gene_list <- get_consensus_genes(deg, min_patients = 8)
  cat("Selected", length(gene_list), "genes for enrichment.\n")

  # KEGG enrichment
  cat("Running KEGG enrichment...\n")
  kegg_df <- run_kegg_enrichment(gene_list)

  # Plot top pathways
  print(plot_kegg_results(kegg_df))
}

# Execute
main()
