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

# % change
### comparing biomarker percentage change from baseline timepoint to timepoint (X) between every group
p_val_pct_list <- anova_pct_list <- list()

for (df in 1:length(csv_df_list)){

  # selecting working df, extracting column names and biomarker name. Creating list of unique timepoints
  temp_df <- as.data.frame(csv_df_list[df])
  col_names <- colnames(temp_df)
  biomarker_name <- col_names[4]
  timepoints <- unique(temp_df[,'time'])
  biomarker_formula <- as.formula(paste0(biomarker_name, " ~ group"))

  # converting df to wide format
  change_wide <- temp_df %>%
    pivot_wider(id_cols = c('subject', 'group'), names_from = 'time', values_from = as.character(biomarker_name))

  # selecting baseline time as 0 or closest -ve number to 0
  baseline_time <- -10000
  if (0 %in% timepoints){
    baseline_time <- 0
  } else {
    for (i in 1:length(timepoints)){
      if (timepoints[i] < 0 & timepoints[i] > baseline_time){
        baseline_time <- timepoints[i]
      }
    }
  }

  # calculating pct change from baseline to every timepoint (data generated into non baseline column as a proportion)
  timepoints_no_bl <- timepoints[(baseline_time!=timepoints)]
  for (t in 1:length(timepoints_no_bl)){
    if (timepoints_no_bl[t] > baseline_time){
      change_wide[,paste0(timepoints_no_bl[t])] <- ((change_wide[,as.character(timepoints_no_bl[t])] - change_wide[,as.character(baseline_time)])/change_wide[,as.character(baseline_time)])
    } else {
      change_wide[,paste0(timepoints_no_bl[t])] <- ((change_wide[,as.character(baseline_time)] - change_wide[,as.character(timepoints_no_bl[t])])/change_wide[,as.character(timepoints_no_bl[t])])
    }
  }

  # converting back to long format and grouping by time
  pct_change_set <- change_wide %>%
    dplyr::select(-c(as.character(baseline_time))) %>%
    pivot_longer(cols = num_range('', -1000000:1000000), names_to = 'time', values_to = biomarker_name) %>%
    group_by(time) %>%
    mutate(group = as.factor(group), baseline = baseline_time)
  pct_change_set[sapply(pct_change_set, is.infinite)] <- NA

  # STATISTICAL TESTS
  group_mods <- pct_change_set %>%
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


  # amending lists with relevant data
  p_val_pct_list[[df]] <- group_mods_res %>% unnest(emmeans_tidy) %>%
    dplyr::select(time, term, contrast, estimate, std.error, adj.p.value) %>%
    mutate(p.adj = ifelse(adj.p.value < 0.001, "p<0.001", paste0("", format(round(adj.p.value, 3), nsmall=3)))) %>%
    mutate(sig = ifelse(adj.p.value < 0.05, "*", '')) %>%
    dplyr::select(-c('adj.p.value')) %>%
    mutate(time = ifelse(baseline_time > time, paste0(time, ' -> ', baseline_time), paste0(baseline_time, ' -> ', time)), biomarker = biomarker_name)

  anova_pct_list[[df]] <- group_mods_res %>% unnest(aov_tidy) %>%
    dplyr::select(time, term, statistic, df, p.value) %>%
    mutate(p.adj = ifelse(p.value < 0.001, "p<0.001", paste0("", format(round(p.value, 3), nsmall=3)))) %>%
    mutate(sig = ifelse(p.value < 0.05, "*", '')) %>%
    dplyr::select(-c('p.value')) %>%
    mutate(time = ifelse(baseline_time > time, paste0(time, ' -> ', baseline_time), paste0(baseline_time, ' -> ', time)), biomarker = biomarker_name)

}

# merging and exporting data and p-value tables
pct_p_merged_df <- p_val_pct_list %>%
  reduce(full_join) %>%
  dplyr::select(-c("term"))
write.csv(pct_p_merged_df, paste0(parent_dir, "/1tables/pct_pvals.csv"), row.names = FALSE)

pct_p_merged_df <- anova_pct_list %>%
  reduce(full_join)
write.csv(pct_p_merged_df, paste0(parent_dir, "/1tables/pct_anova.csv"), row.names = FALSE)


# absoulte change
### comparing biomarker percentage change from baseline timepoint to timepoint (X) between every group
p_val_abs_list <- anova_abs_list <- list()

for (df in 1:length(csv_df_list)){
  
  # selecting working df, extracting column names and biomarker name. Creating list of unique timepoints
  temp_df <- as.data.frame(csv_df_list[df])
  col_names <- colnames(temp_df)
  biomarker_name <- col_names[4]
  timepoints <- unique(temp_df[,'time'])
  biomarker_formula <- as.formula(paste0(biomarker_name, " ~ group"))
  
  # converting df to wide format
  change_wide <- temp_df %>%
    pivot_wider(id_cols = c('subject', 'group'), names_from = 'time', values_from = as.character(biomarker_name))
  
  # selecting baseline time as 0 or closest -ve number to 0
  baseline_time <- -10000
  if (0 %in% timepoints){
    baseline_time <- 0
  } else {
    for (i in 1:length(timepoints)){
      if (timepoints[i] < 0 & timepoints[i] > baseline_time){
        baseline_time <- timepoints[i]
      }
    }
  }
  
  # calculating pct change from baseline to every timepoint (data generated into non baseline column as a proportion)
  timepoints_no_bl <- timepoints[(baseline_time!=timepoints)]
  for (t in 1:length(timepoints_no_bl)){
    if (timepoints_no_bl[t] > baseline_time){
      change_wide[,paste0(timepoints_no_bl[t])] <- (change_wide[,as.character(timepoints_no_bl[t])] - change_wide[,as.character(baseline_time)])
    } else {
      change_wide[,paste0(timepoints_no_bl[t])] <- (change_wide[,as.character(baseline_time)] - change_wide[,as.character(timepoints_no_bl[t])])
    }
  }
  
  # converting back to long format and grouping by time
  abs_change_set <- change_wide %>%
    dplyr::select(-c(as.character(baseline_time))) %>%
    pivot_longer(cols = num_range('', -1000000:1000000), names_to = 'time', values_to = biomarker_name) %>%
    group_by(time) %>%
    mutate(group = as.factor(group), baseline = baseline_time)
  abs_change_set[sapply(abs_change_set, is.infinite)] <- NA
  
  # STATISTICAL TESTS
  group_mods <- abs_change_set %>%
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
  
  
  # amending lists with relevant data
  p_val_abs_list[[df]] <- group_mods_res %>% unnest(emmeans_tidy) %>%
    dplyr::select(time, term, contrast, estimate, std.error, adj.p.value) %>%
    mutate(p.adj = ifelse(adj.p.value < 0.001, "p<0.001", paste0("", format(round(adj.p.value, 3), nsmall=3)))) %>%
    mutate(sig = ifelse(adj.p.value < 0.05, "*", '')) %>%
    dplyr::select(-c('adj.p.value')) %>%
    mutate(time = ifelse(baseline_time > time, paste0(time, ' -> ', baseline_time), paste0(baseline_time, ' -> ', time)), biomarker = biomarker_name)
  
  anova_abs_list[[df]] <- group_mods_res %>% unnest(aov_tidy) %>%
    dplyr::select(time, term, statistic, df, p.value) %>%
    mutate(p.adj = ifelse(p.value < 0.001, "p<0.001", paste0("", format(round(p.value, 3), nsmall=3)))) %>%
    mutate(sig = ifelse(p.value < 0.05, "*", '')) %>%
    dplyr::select(-c('p.value')) %>%
    mutate(time = ifelse(baseline_time > time, paste0(time, ' -> ', baseline_time), paste0(baseline_time, ' -> ', time)), biomarker = biomarker_name)
  
}

# merging and exporting data and p-value tables
abs_p_merged_df <- p_val_abs_list %>%
  reduce(full_join) %>%
  dplyr::select(-c("term"))
write.csv(abs_p_merged_df, paste0(parent_dir, "/1tables/abs_pvals.csv"), row.names = FALSE)

abs_p_merged_df <- anova_abs_list %>%
  reduce(full_join)
write.csv(abs_p_merged_df, paste0(parent_dir, "/1tables/abs_anova.csv"), row.names = FALSE)