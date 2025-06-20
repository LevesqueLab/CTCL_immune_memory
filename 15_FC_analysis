#!/usr/bin/env Rscript

#' CTCL TCM Fraction Analysis (Screen + Longitudinal Cohorts)
#'
#' This script loads flow cytometry quantifications of mTRM/mTCM and Malignant/Healthy,
#' visualizes them as stacked barplots and summarizes longitudinal changes across timepoints.
#'
#' Dependencies: ggplot2, dplyr, openxlsx, tidyr

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(openxlsx)
  library(tidyr)
})

# Define plotting function
plot_stacked_bar <- function(data, sample_order, type, colors) {
  data %>%
    select(Sample, Timepoint, Response.Blood, all_of(type)) %>%
    pivot_longer(cols = all_of(type), names_to = "Cell_Type", values_to = "Count") %>%
    group_by(Sample, Timepoint) %>%
    mutate(Percentage = Count / sum(Count) * 100) %>%
    ungroup() %>%
    mutate(
      Sample = factor(Sample, levels = sample_order),
      Response.Blood = factor(Response.Blood, levels = c("yes", "no"))
    ) %>%
    ggplot(aes(x = interaction(Timepoint, Sample), y = Percentage, fill = Cell_Type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = colors) +
    labs(x = "", y = "Percentage", fill = "Cell Type") +
    theme_minimal() +
    theme(text = element_text(size = 14), axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none")
}

# Load screen cohort
screen_df <- read.xlsx("/Screen_typing_CD27_CD28.xlsx", sheet = 1)
screen_df$Sample <- factor(screen_df$Sample, levels = paste("Sézary", 1:15))

# Process and plot screen cohort barplot
screen_long <- screen_df %>%
  select(Sample, `mTCM`, `mTRM`) %>%
  pivot_longer(cols = everything(), names_to = "Category", values_to = "Count") %>%
  group_by(Sample) %>%
  mutate(Proportion = Count / sum(Count))

ggplot(screen_long, aes(x = Sample, y = Proportion, fill = Category)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_trans(y = "sqrt") +
  labs(title = "Screen Cohort: TCM Fractions per Sample",
       x = "", y = "Proportion (sqrt)", fill = "Category") +
  scale_fill_manual(values = c("mTCM" = "darkgreen", "mTRM" = "darkred")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Load longitudinal cohort
long_df <- read.xlsx("/Longitudinal_typing_CD27_CD28.xlsx", sheet = 1)

long_df <- long_df %>%
  mutate(
    Malignancy_load = Malignant / Healthy,
    TCM_load = `mTCM` / `mTRM`,
    log_Malignancy = log(Malignancy_load),
    log_TCM = log(TCM_load)
  )

# Plot ECP stacked barplots
ECP <- long_df %>% filter(Cohort == "ECP") %>% arrange(Timepoint)
plot_stacked_bar(ECP, sample_order = c("ECP1", "ECP3", "ECP6", "ECP2", "ECP4", "ECP5"),
                 type = c("Malignant", "Healthy"),
                 colors = c("Malignant" = "darkblue", "Healthy" = "orange"))
plot_stacked_bar(ECP, sample_order = c("ECP1", "ECP3", "ECP6", "ECP2", "ECP4", "ECP5"),
                 type = c("mTCM", "mTRM"),
                 colors = c("mTCM" = "darkgreen", "mTRM" = "darkred"))

# Plot Moga stacked barplots
Moga <- long_df %>% filter(Cohort == "Moga") %>% arrange(Timepoint)
plot_stacked_bar(Moga, sample_order = c("Moga2", "Moga1", "Moga3", "Moga4", "Moga5", "Moga6"),
                 type = c("Malignant", "Healthy"),
                 colors = c("Malignant" = "darkblue", "Healthy" = "orange"))
plot_stacked_bar(Moga, sample_order = c("Moga2", "Moga1", "Moga3", "Moga4", "Moga5", "Moga6"),
                 type = c("mTCM", "mTRM"),
                 colors = c("mTCM" = "darkgreen", "mTRM" = "darkred"))

# Combine ECP and Moga, normalize per sample baseline
all <- bind_rows(ECP, Moga) %>%
  arrange(Cohort, Sample, Timepoint) %>%
  group_by(Sample) %>%
  mutate(
    norm_Malignancy = Malignancy_load / first(Malignancy_load),
    norm_TCM = TCM_load / first(TCM_load),
    log_Malignancy = log(norm_Malignancy),
    log_TCM = log(norm_TCM)
  ) %>%
  ungroup() %>%
  mutate(Response.Blood = factor(Response.Blood, levels = c("yes", "no")))

# Summary trajectory plot
ggplot(all, aes(x = log_TCM, y = log_Malignancy, color = Response.Blood, shape = Timepoint)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("yes" = "darkblue", "no" = "darkred")) +
  geom_path(aes(group = Sample, color = Response.Blood)) +
  labs(title = "Normalized log(Malignancy) vs log(TCM) Load",
       x = "log(mTCM / mTRM) normalized to BL",
       y = "log(Malignant / Healthy) normalized to BL",
       color = "Response",
       shape = "Timepoint") +
  facet_wrap(~Cohort) +
  theme_minimal()
