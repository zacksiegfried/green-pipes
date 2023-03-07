library(tidyverse)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(this.path)
library(arsenal)

### This script should be placed in a directory with standardized csv data sets and ran to:
### GENERATE A NEW SUBDIRECTORY CONTAINING P-VALUE TABLES (ANOVA and Tukey)

##### LOADING DATA
# sets working directory as where this script is, constructs a list of file path names matching *std.csv in directory
parent_dir <- this.path::here(..=0)
setwd(paste0(parent_dir))
file_path_list <- list.files(pattern="*std.csv")

# reads in csv files and constructs a list of data frames
csv_df_list <- list()
for (i in 1:length(file_path_list)){
  csv_df_list[[i]] <- read.csv(file_path_list[[i]])
}

# creating directories to save tables
dir.create(paste0(parent_dir,"/1tables"))


##### CALCULATING STATISTICAL TEST P-VALUES

### comparison of biomarker mean between every group for a fixed timepoint
p_val_group_list <- anova_group_list <- list()

for (df in 1:length(csv_df_list)){

  # selecting working df, extracting column names and biomarker name
  temp_df <- as.data.frame(csv_df_list[df])
  temp_df$group <- as.factor(temp_df$group)
  col_names <- colnames(temp_df)
  biomarker_name <- col_names[4]
  timepoints <- unique(temp_df[,'time'])
  biomarker_formula <- as.formula(paste0(biomarker_name, " ~ group"))


  # STATISTICAL TESTS
  group_mods <- temp_df %>%
    group_by(time) %>%
    nest() %>%
    mutate(fit = map(data, ~ glm(biomarker_formula, data = .x)))

  group_mods_res <- group_mods %>%
    mutate(# Get Pairwise Comparisons across groups within each time point
      means= map(fit, ~pairs(emmeans::emmeans(.x,~group,data=.x$data))),
      # Get overall ANOVA effect of group
      aov=map(fit, ~car::Anova(.x, 3))) %>%
    mutate(emmeans_tidy = map(means, tidy)) %>%
    mutate(aov_tidy = map(aov, tidy))


  p_val_group_list[[df]] <- group_mods_res %>% unnest(emmeans_tidy) %>%
    dplyr::select(time, term, contrast, estimate, std.error, adj.p.value) %>%
    mutate(p.adj = ifelse(adj.p.value < 0.001, "p<0.001", paste0("", format(round(adj.p.value, 3), nsmall=3))), biomarker = biomarker_name)

  anova_group_list[[df]] <- group_mods_res %>% unnest(aov_tidy) %>%
    dplyr::select(time, term, statistic, df, p.value) %>%
    mutate(p.adj = ifelse(p.value < 0.001, "p<0.001", paste0("", format(round(p.value, 3), nsmall=3))), biomarker = biomarker_name)
}

group_p_merged_df <- p_val_group_list %>%
  reduce(full_join) %>%
  dplyr::select(-c("term"))
write.csv(group_p_merged_df, paste0(parent_dir, "/1tables/group_pvals.csv"), row.names = FALSE)

group_p_merged_df <- anova_group_list %>%
  reduce(full_join)
write.csv(group_p_merged_df, paste0(parent_dir, "/1tables/group_anova.csv"), row.names = FALSE)
