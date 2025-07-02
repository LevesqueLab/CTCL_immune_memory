#!/usr/bin/env Rscript

#' Statistical Analysis Pipeline
#'
#' This script performs statistical testing on 
#' CytoTRACE scores, UCell scores, tumor mutational burden and relative abundance.
#' Includes NMDS, PERMANOVA, SIMPER, mixed-effect modeling, and binomial testing.

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(vegan)
  library(ggpubr)
})

# -----------------------------
# CytoTRACE / UCell Score Analysis
# -----------------------------

# Load and model score data
variable <- CytoTRACE2_Score  # or UCell score
final_results <- final_results  # Output from scoring step

# Linear mixed-effects model
model <- lmer(variable ~ TCM * tissue + (1 | patient), data = final_results)
anova(model)

# Post-hoc pairwise tests with FDR correction
pairwise_results <- emmeans(model, pairwise ~ TCM * tissue, adjust = "fdr")
summary(pairwise_results)

# -----------------------------
# Tumor Mutational Burden (TMB) Analysis
# -----------------------------

model <- lmer(tmb ~ Cell_types + (1 | sample), data = final_results_sub)
anova(model)

pairwise_results <- emmeans(model, pairwise ~ Cell_types, adjust = "fdr")
summary(pairwise_results)

# -----------------------------
# Relative Abundance scRNA-seq Malignant T CTCL blood and skin
# -----------------------------

temp <- table(seu_Mal$patient, seu_Mal$TCM)
temp_2 <- as.data.frame(temp)
df_wide <- pivot_wider(temp_2, names_from = Var2, values_from = Freq)
names(df_wide)[names(df_wide) == "Var1"] <- "Patient"
df_wide$tissue <- c(rep("blood", 7), rep("skin", 9))

# Convert counts to percentages
cell_types <- c("mCyto", "mTCM", "mTRM")
df_percentages <- df_wide %>%
  mutate(Total = rowSums(across(all_of(cell_types)))) %>%
  mutate(across(all_of(cell_types), ~ .x / Total * 100))

data <- df_percentages
rownames(data) <- NULL

# NMDS ordination
md_raw <- data %>%
  column_to_rownames("Patient") %>%
  select(any_of(cell_types)) %>%
  metaMDS(distance = "bray", autotransform = FALSE)

# Compute NMDS centroids
nm <- as.data.frame(md_raw$points) %>%
  rownames_to_column("Patient") %>%
  merge(data) %>%
  unique() %>%
  group_by(tissue) %>%
  mutate(
    cent.1 = mean(MDS1),
    cent.2 = mean(MDS2)
  )

# Plot NMDS
library(ggplot2)
ggplot(nm, aes(x = MDS1, y = MDS2, color = tissue)) +
  geom_point(size = 4) +
  theme_bw() +
  geom_segment(aes(xend = cent.1, yend = cent.2)) +
  scale_color_manual(values = c("blood" = "#DC143C", "skin" = "#8B4513"))

# PERMANOVA
for_perm <- data %>% select(any_of(cell_types))
perm <- adonis2(for_perm ~ tissue, data, method = "bray")
print(perm)

# SIMPER analysis
sim <- simper(for_perm, group = data$tissue)
summary(sim)

# plot per tissue
clr_temp <- clr_pt_long %>%
  filter(cell == "mTCM") %>%
  mutate(fraction = fraction / 100)

a <- ggboxplot(clr_temp, x = "tissue",
               y = "fraction",
               fill = "tissue") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("mTCM") +
  ylab("Proportion") +
  scale_fill_manual(values = c("blood" = "#DC143C", "skin" = "#8B4513"))

clr_temp <- clr_pt_long %>%
  filter(cell == "mTRM") %>%
  mutate(fraction = fraction / 100)

b <- ggboxplot(clr_temp, x = "tissue",
               y = "fraction",
               fill = "tissue") +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("mTRM") +
  ylab("Proportion") +
  scale_fill_manual(values = c("blood" = "#DC143C", "skin" = "#8B4513"))

print(a | b)

# -----------------------------
# Ratio Analysis (mTCM/mTRM)
# -----------------------------

temp <- table(seu_Xe$sample, seu_Xe$rctd_celltype)
temp_2 <- as.data.frame(temp) %>% filter(Var2 %in% c("mTRM", "mTCM"))
df_wide <- pivot_wider(temp_2, names_from = Var2, values_from = Freq)
names(df_wide)[names(df_wide) == "Var1"] <- "donor_id"

df_percentages <- df_wide %>%
  mutate(Total = rowSums(across(c("mTCM", "mTRM")))) %>%
  mutate(
    mTCM = mTCM / Total * 100,
    mTRM = mTRM / Total * 100
  )

# Binomial test: mTRM > mTCM?
successes <- sum(df_percentages$mTRM > df_percentages$mTCM, na.rm = TRUE)
n <- sum(!is.na(df_percentages$mTRM) & !is.na(df_percentages$mTCM))
binom.test(successes, n, p = 0.5, alternative = "two.sided")

# End of script
cat("Statistical analysis completed.\n")
