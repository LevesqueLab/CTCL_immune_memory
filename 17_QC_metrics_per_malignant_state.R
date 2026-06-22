#!/usr/bin/env Rscript

#' CTCL QC Metrics per Malignant State (mTCM / mTRM / mCyto)
#'
#' This script summarises quality-control metrics per malignant T-cell state on
#' the CTCL_blood_skin_malignant Seurat object: transcripts per cell
#' (nCount_RNA), genes per cell (nFeature_RNA), fraction of cells with a
#' productive TCR (CTstrict populated), and TCR-chain RNA expression
#' (TRAC, TRBC1/2). Results are tabulated overall, by tissue, and blood-only,
#' with pairwise Wilcoxon tests between states and 4-panel publication figures.
#'
#' Produced for Reviewer Reply 12 of the Cell Reports Medicine revision, to
#' quantify the degree of TCR downregulation in mTRM cells relative to mTCM.
#'
#' Dependencies: Seurat, qs2, dplyr, tidyr, ggplot2, patchwork, rstatix, ggpubr

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(rstatix); library(ggpubr)
})

setwd("/srv/GT/analysis/p28409-Andrea/Revisions_CellRepMed/QC_mTRM_vs_mTCM")

cat("Loading malignant Seurat object ...\n")
obj <- qs_read("/srv/GT/analysis/p28409-Andrea/switchdrive_data/CTCL_blood_skin_malignant.qs2",
               nthreads = 4)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}
cat(sprintf("Cells: %d   Genes: %d\n", ncol(obj), nrow(obj)))
cat("Available metadata cols:\n"); print(colnames(obj@meta.data))

## --- Build the per-cell table -------------------------------------------------
md <- obj@meta.data
md$state  <- factor(md$Malignant_T, levels = c("mTCM", "mTRM", "mCyto"))
md$tissue <- factor(md$tissue, levels = c("blood", "skin"))

# TCR detection - CTstrict populated => productive TCR recovered
md$has_TCR <- !is.na(md$CTstrict) & nchar(as.character(md$CTstrict)) > 0

# TRAC / TRBC expression (log-normalised values)
DefaultAssay(obj) <- "RNA"
tcr_genes <- intersect(c("TRAC", "TRBC1", "TRBC2"), rownames(obj))
cat("TCR genes present:", paste(tcr_genes, collapse = ", "), "\n")
expr <- FetchData(obj, vars = tcr_genes, layer = "data")
md   <- cbind(md, expr)
md$TRBC <- pmax(md$TRBC1, md$TRBC2)   # combine TRBC1/2 to single TRBC channel

## --- Per-state summary --------------------------------------------------------
summarise_qc <- function(df, group_cols) {
  df |>
    dplyr::group_by(across(all_of(group_cols))) |>
    dplyr::summarise(
      n_cells          = dplyr::n(),
      median_nCount    = median(nCount_RNA),
      iqr_nCount       = paste(round(quantile(nCount_RNA, c(.25,.75))), collapse = "-"),
      mean_nCount      = round(mean(nCount_RNA), 0),
      sd_nCount        = round(sd(nCount_RNA), 0),
      median_nFeature  = median(nFeature_RNA),
      iqr_nFeature     = paste(round(quantile(nFeature_RNA, c(.25,.75))), collapse = "-"),
      mean_nFeature    = round(mean(nFeature_RNA), 0),
      sd_nFeature      = round(sd(nFeature_RNA), 0),
      pct_TCR_productive = round(100 * mean(has_TCR), 1),
      mean_TRAC        = round(mean(TRAC), 3),
      mean_TRBC        = round(mean(TRBC), 3),
      pct_TRAC_pos     = round(100 * mean(TRAC > 0), 1),
      pct_TRBC_pos     = round(100 * mean(TRBC > 0), 1),
      .groups = "drop"
    )
}

sum_state         <- summarise_qc(md, "state")
sum_state_tissue  <- summarise_qc(md, c("state", "tissue"))
sum_blood_only    <- summarise_qc(dplyr::filter(md, tissue == "blood"), "state")

write.csv(sum_state,        "QC_metrics_by_state.csv",         row.names = FALSE)
write.csv(sum_state_tissue, "QC_metrics_by_state_tissue.csv",  row.names = FALSE)
write.csv(sum_blood_only,   "QC_metrics_blood_only.csv",       row.names = FALSE)
cat("Summary saved.\n")

## --- Wilcoxon tests -----------------------------------------------------------
wx <- function(df, var) {
  df |>
    rstatix::wilcox_test(as.formula(paste(var, "~ state"))) |>
    rstatix::adjust_pvalue(method = "BH") |>
    dplyr::mutate(metric = var)
}
stats_combined <- dplyr::bind_rows(
  wx(md, "nCount_RNA"),
  wx(md, "nFeature_RNA"),
  wx(md, "TRAC"),
  wx(md, "TRBC")
)
write.csv(stats_combined, "stats_wilcox.csv", row.names = FALSE)

## --- Figures ------------------------------------------------------------------
state_pal <- c(mTCM = "darkgreen", mTRM = "darkred", mCyto = "darkblue")
comps     <- list(c("mTCM", "mTRM"), c("mTCM", "mCyto"), c("mTRM", "mCyto"))

make_violin <- function(df, var, ylab, log_y = FALSE) {
  p <- ggplot(df, aes(x = state, y = .data[[var]], fill = state)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.75) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.7) +
    scale_fill_manual(values = state_pal, guide = "none") +
    stat_compare_means(comparisons = comps, method = "wilcox.test",
                       label = "p.format", size = 3.2) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(face = "bold"))
  if (log_y) p <- p + scale_y_log10()
  p
}

p_n  <- make_violin(md, "nCount_RNA",   "Transcripts / cell (nCount_RNA)", log_y = TRUE)
p_g  <- make_violin(md, "nFeature_RNA", "Genes / cell (nFeature_RNA)")

# %TCR productive: bar chart
tcr_bar <- md |> dplyr::group_by(state) |>
  dplyr::summarise(pct = 100 * mean(has_TCR), n = dplyr::n(), .groups = "drop")
p_t <- ggplot(tcr_bar, aes(x = state, y = pct, fill = state)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct, n)),
            vjust = -0.4, size = 3.4, fontface = "bold") +
  scale_fill_manual(values = state_pal, guide = "none") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "% cells with productive TCR") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(face = "bold"))

# TRAC + TRBC expression
expr_long <- md |>
  dplyr::select(state, TRAC, TRBC) |>
  tidyr::pivot_longer(c(TRAC, TRBC), names_to = "gene", values_to = "expr")
p_e <- ggplot(expr_long, aes(x = state, y = expr, fill = state)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.75) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.7) +
  facet_wrap(~ gene, nrow = 1) +
  scale_fill_manual(values = state_pal, guide = "none") +
  stat_compare_means(comparisons = comps, method = "wilcox.test",
                     label = "p.format", size = 2.8) +
  labs(x = NULL, y = "log-normalised expression") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(face = "bold"),
        strip.text  = element_text(face = "bold"))

combined <- (p_n | p_g) / (p_t | p_e) +
  patchwork::plot_annotation(
    title = "QC metrics per malignant state (CTCL_blood_skin_malignant)",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave("QC_metrics_by_state.png", combined, width = 11, height = 9,
       dpi = 300, bg = "white")
ggsave("QC_metrics_by_state.pdf", combined, width = 11, height = 9,
       device = cairo_pdf, bg = "white")

## Same, split by tissue (blood vs skin) -- shows confound stays under control
md_bt <- md
p_n2 <- make_violin(md_bt, "nCount_RNA", "Transcripts / cell", log_y = TRUE) +
  facet_wrap(~ tissue, nrow = 1) + theme(strip.text = element_text(face = "bold"))
p_g2 <- make_violin(md_bt, "nFeature_RNA", "Genes / cell") +
  facet_wrap(~ tissue, nrow = 1) + theme(strip.text = element_text(face = "bold"))
tissue_combined <- p_n2 / p_g2 +
  patchwork::plot_annotation(
    title = "QC metrics per malignant state, split by tissue",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )
ggsave("QC_metrics_by_state_tissue.png", tissue_combined, width = 11, height = 9,
       dpi = 300, bg = "white")
ggsave("QC_metrics_by_state_tissue.pdf", tissue_combined, width = 11, height = 9,
       device = cairo_pdf, bg = "white")

cat("\nSummary table (per state):\n"); print(sum_state)
cat("\nDONE\n")
