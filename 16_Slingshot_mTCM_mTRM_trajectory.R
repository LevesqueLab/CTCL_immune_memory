#!/usr/bin/env Rscript

#' CTCL mTCM -> mTRM Slingshot Trajectory Heatmaps
#'
#' This script runs Slingshot + tradeSeq on the malignant T-cell object
#' (CytoTRACE2-annotated) to fit a mTCM -> mTRM pseudotime trajectory, predicts
#' smoothed expression along that trajectory for a curated 19-category gene
#' panel, and renders a heatmap with two identity bars (Malignant_T state and
#' tissue) smoothed along pseudotime.
#'
#' Paper colours: mTCM = darkgreen, mTRM = darkred, mCyto = darkblue.
#'
#' Dependencies: Seurat, qs, slingshot, tradeSeq, ggplot2, dplyr, tidyr,
#'   patchwork, scales, SingleCellExperiment, readxl

suppressPackageStartupMessages({
  library(Seurat); library(qs); library(slingshot); library(tradeSeq)
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)
  library(SingleCellExperiment)
})

setwd("/srv/GT/analysis/p28409-Andrea/Memory_Revision_Slingshot")

# ---------- LOAD ------------------------------------------------------------
cat("Loading seurat object ...\n")
seurat_obj <- qs::qread("Malignant_T_cytotrace.qs", nthreads = 8)

# Map original TCM labels (TCM-NP / TCM-P / TCM-CT) → paper labels (mTCM / mTRM / mCyto)
tcm_map <- c("TCM-NP" = "mTCM", "TCM-P" = "mTRM", "TCM-CT" = "mCyto")
seurat_obj$Malignant_T <- factor(unname(tcm_map[as.character(seurat_obj$TCM)]),
                                 levels = c("mTCM", "mTRM", "mCyto"))

# ---------- SLINGSHOT -------------------------------------------------------
cat("Running slingshot ...\n")
dimred <- Embeddings(seurat_obj, "umap")
Idents(seurat_obj) <- seurat_obj$Malignant_T  # use paper labels for slingshot too
clusters <- Idents(seurat_obj)
set.seed(1)
lineages <- getLineages(data = dimred, clusterLabels = clusters,
                        start.clus = "mTCM", end.clus = "mTRM")
lineages <- as.SlingshotDataSet(lineages)
curves <- getCurves(lineages, thresh = 0.01, stretch = 0.8,
                    allow.breaks = FALSE, shrink = 0.8)
curves <- as.SlingshotDataSet(curves)

# ---------- COUNTS + GENE SELECTION ----------------------------------------
DefaultAssay(seurat_obj) <- "RNA"
counts <- GetAssayData(seurat_obj, layer = "counts")

# Andrea's curated 19-category trajectory gene list (from
# 2026-02-18-FlashSeq-Blood-Integration/gene list for trajectory.xlsx). Each
# column header is a functional category; cells under it are the gene symbols.
suppressPackageStartupMessages(library(readxl))
xlsx_path <- "/srv/GT/analysis/p28409-Andrea/Memory_Revision_Slingshot/gene_list_trajectory.xlsx"
gene_xl <- as.data.frame(read_excel(xlsx_path, sheet = 1))
gene_categories <- unlist(lapply(colnames(gene_xl), function(cat) {
  gs <- unique(stats::na.omit(trimws(as.character(gene_xl[[cat]]))))
  gs <- gs[nchar(gs) > 0]
  setNames(rep(cat, length(gs)), gs)
}))
gene_categories <- gene_categories[!duplicated(names(gene_categories))]
genes_of_interest <- unique(names(gene_categories))
cat("Categories:", length(unique(gene_categories)), "  Total genes:", length(genes_of_interest), "\n")
genes_found <- intersect(genes_of_interest, rownames(counts))
cat("Genes found:", length(genes_found), "/", length(genes_of_interest), "\n")
counts_subset <- counts[genes_found, ]
counts_subset <- counts_subset[!duplicated(rownames(counts_subset)), ]

# ---------- FIT GAM ---------------------------------------------------------
cat("fitGAM ...\n")
sce_subset <- fitGAM(counts = as.matrix(counts_subset), sds = curves, verbose = FALSE,
                     parallel = FALSE)
qs::qsave(sce_subset, "cache_sce_subset.qs", nthreads = 8)
cat("Saved cache_sce_subset.qs\n")

# ---------- PREDICT SMOOTH --------------------------------------------------
cat("predictSmooth ...\n")
smooth_df <- predictSmooth(sce_subset, gene = rownames(counts_subset),
                           nPoints = 100, tidy = TRUE)

# gene_categories is built from the xlsx earlier (one category per spreadsheet column)
smooth_df$category <- gene_categories[as.character(smooth_df$gene)]
smooth_df <- smooth_df |> dplyr::filter(lineage == 1)
smooth_df <- smooth_df |> dplyr::group_by(gene) |>
  dplyr::mutate(expr_scaled = as.numeric(scale(yhat))) |> dplyr::ungroup()

# ---------- CYTOTRACE2 SMOOTH ----------------------------------------------
pt_full <- slingPseudotime(curves)[, 1]
pt_full <- pt_full[!is.na(pt_full)]
cyto <- seurat_obj$CytoTRACE2_Score[names(pt_full)]
loess_fit <- loess(cyto ~ pt_full, span = 0.3)
time_grid <- seq(min(pt_full), max(pt_full), length.out = 100)
cyto_smooth <- predict(loess_fit, newdata = time_grid)
cyto_df <- data.frame(
  gene = "CytoTRACE2", time = time_grid, yhat = cyto_smooth, lineage = 1,
  category = "Stemness score", expr_scaled = as.numeric(scale(cyto_smooth)))

smooth_df$gene <- as.character(smooth_df$gene)
smooth_df_combined <- dplyr::bind_rows(smooth_df, cyto_df)

# ---------- IDENTITY BARS (smoothed along pseudotime) ----------------------
cat("Building identity bars ...\n")
state_pal  <- c(mTCM = "darkgreen", mTRM = "darkred", mCyto = "darkblue")
tissue_pal <- c(blood = "#3C5488", skin = "#E64B35")

states  <- as.character(seurat_obj$Malignant_T[names(pt_full)])
tissues <- as.character(seurat_obj$tissue[names(pt_full)])

# Gaussian-kernel smoothed fraction per time grid point
bandwidth <- diff(range(pt_full)) / 25  # ~4% of pt range
compute_frac <- function(labels, levels_keep) {
  out <- matrix(0, nrow = length(time_grid), ncol = length(levels_keep),
                dimnames = list(NULL, levels_keep))
  for (i in seq_along(time_grid)) {
    w <- stats::dnorm(pt_full, mean = time_grid[i], sd = bandwidth)
    if (sum(w) == 0) next
    w <- w / sum(w)
    for (lv in levels_keep) out[i, lv] <- sum(w[labels == lv])
  }
  as.data.frame(out) |>
    dplyr::mutate(time = time_grid) |>
    tidyr::pivot_longer(-time, names_to = "category", values_to = "frac")
}

state_bar_df  <- compute_frac(states,  c("mTCM", "mTRM", "mCyto"))
tissue_bar_df <- compute_frac(tissues, c("blood", "skin"))

# Time grid step width for geom_col bar width
time_w <- diff(time_grid)[1] * 0.97

p_state_bar <- ggplot(state_bar_df, aes(x = time, y = frac,
                                        fill = factor(category, levels = c("mTCM","mTRM","mCyto")))) +
  geom_col(width = time_w, position = "stack") +
  scale_fill_manual(values = state_pal, name = "State",
                    drop = FALSE, guide = guide_legend(order = 1)) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), limits = c(min(time_grid)-time_w/2, max(time_grid)+time_w/2)) +
  labs(y = "State", x = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5, hjust = 1, margin = margin(r = 4)),
        axis.text.x = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(2, 2, 0, 2))

p_tissue_bar <- ggplot(tissue_bar_df, aes(x = time, y = frac, fill = category)) +
  geom_col(width = time_w, position = "stack") +
  scale_fill_manual(values = tissue_pal, name = "Tissue",
                    guide = guide_legend(order = 2)) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), limits = c(min(time_grid)-time_w/2, max(time_grid)+time_w/2)) +
  labs(y = "Tissue", x = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5, hjust = 1, margin = margin(r = 4)),
        axis.text.x = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0, 2, 0, 2))

# ---------- MAIN HEATMAP (FULL) --------------------------------------------
make_main <- function(df) {
  # Order genes: CytoTRACE2 first, then by xlsx column order
  cat_order <- c("Stemness score", colnames(gene_xl))
  cats_in_data <- intersect(cat_order, unique(df$category))
  df$category <- factor(df$category, levels = cats_in_data)
  ggplot(df, aes(x = time, y = gene, fill = expr_scaled)) +
    geom_tile() +
    scale_fill_gradient2(low = "#00428c", mid = "white", high = "#9e002a",
                         midpoint = 0, limits = c(-2, 2), oob = scales::squish,
                         name = "Scaled\nExpression") +
    facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = "Pseudotime", y = "Gene") +
    theme_bw(base_size = 11) +
    theme(strip.text.y = element_text(angle = 0, hjust = 0),
          axis.text.y = element_text(size = 8),
          panel.spacing = unit(0.2, "lines"),
          plot.margin = margin(0, 2, 2, 2))
}

p_main_full <- make_main(smooth_df_combined)

cat("Saving full heatmap PNG ...\n")
combined_full <- (p_state_bar / p_tissue_bar / p_main_full) +
  patchwork::plot_layout(heights = c(0.6, 0.6, 14), guides = "collect") &
  theme(legend.position = "right")
# Height scales with gene count: ~0.12 inch/gene + 2 inch overhead
height_in <- max(11, ceiling(0.12 * length(unique(smooth_df_combined$gene)) + 2))
ggsave("Trajectory_mTCM_mTRM_heatmap_with_bars.png",
       combined_full, width = 13, height = height_in, units = "in",
       dpi = 600, bg = "white", limitsize = FALSE)

# ---------- SELECTION HEATMAP ----------------------------------------------
genes_to_keep <- c("TOX","PDCD1","KLHL42","TSHZ2","SESN3","DNM3",
                   "SETD2","KDM6A","EZH2","TET2","DNMT3A",
                   "CXCR4","S1PR1","KLF2",
                   "RPS6","RPL10","RPL13",
                   "RORA","ITGA1")
smooth_df_sub <- smooth_df |> dplyr::filter(gene %in% genes_to_keep)
p_main_sel <- make_main(smooth_df_sub)

cat("Saving selection heatmap PNG ...\n")
combined_sel <- (p_state_bar / p_tissue_bar / p_main_sel) +
  patchwork::plot_layout(heights = c(0.6, 0.6, 7), guides = "collect") &
  theme(legend.position = "right")
ggsave("Trajectory_mTCM_mTRM_heatmap_selection_with_bars.png",
       combined_sel, width = 7, height = 7, units = "in",
       dpi = 600, bg = "white")

cat("DONE\n")
