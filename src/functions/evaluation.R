# src/functions/evaluation.R
# ==============================================================================
# EVALUATION MODULE
# Description: Computes statistical metrics to evaluate batch correction efficacy
#              and biological variance preservation.
# Dependencies: SingleCellExperiment, dplyr, tidyr, stats
# ==============================================================================

library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(stats)

#' @title Evaluate Batch Correction
#' @description Computes ANOVA F-statistics for batch removal and Spearman 
#'              correlation for biological signal preservation.
#' @param sce_pre Uncorrected SCE object.
#' @param sce_post Corrected SCE object.
#' @param markers Character vector of markers to evaluate.
#' @return A list containing three data frames (F_Statistics, Biological_Preservation, Summary).
evaluate_batch_correction <- function(sce_pre, sce_post, markers) {
  
  message("[Evaluation] Initializing metric calculation...")
  
  # Ensure strict ordering of cells just in case
  if (ncol(sce_pre) != ncol(sce_post)) {
    stop("[Evaluation Fatal] Pre and Post SCE objects have different cell counts.")
  }
  
  # 1. Extract Data
  meta_df <- as.data.frame(colData(sce_pre))
  exprs_pre  <- t(assay(sce_pre, "exprs"))
  exprs_post <- t(assay(sce_post, "exprs"))
  
  # Combine for aggregation
  df_pre  <- bind_cols(as.data.frame(exprs_pre[, markers, drop = FALSE], check.names = FALSE), meta_df)
  df_post <- bind_cols(as.data.frame(exprs_post[, markers, drop = FALSE], check.names = FALSE), meta_df)
  
  # ============================================================================
  # METRIC 1: BATCH EFFECT REMOVAL (ANOVA F-Statistic)
  # ============================================================================
  message("[Evaluation] Computing ANOVA F-Statistics for batch variance...")
  
  # Downsample if matrix is massive (>100k cells) to save time/memory for aov()
  n_cells <- nrow(df_pre)
  max_cells_stat <- 100000
  
  if (n_cells > max_cells_stat) {
    message(sprintf("             Subsampling %d out of %d cells for F-stat computation...", max_cells_stat, n_cells))
    set.seed(42)
    sample_idx <- sample(seq_len(n_cells), max_cells_stat)
    stat_pre <- df_pre[sample_idx, ]
    stat_post <- df_post[sample_idx, ]
  } else {
    stat_pre <- df_pre
    stat_post <- df_post
  }
  
  f_stats_df <- data.frame(Marker = markers, F_Stat_Pre = NA, F_Stat_Post = NA)
  
  for (i in seq_along(markers)) {
    m <- markers[i]
    
    # Pre-correction F-value
    fit_pre <- aov(stat_pre[[m]] ~ stat_pre$batch)
    f_pre <- summary(fit_pre)[[1]][["F value"]][1]
    
    # Post-correction F-value
    fit_post <- aov(stat_post[[m]] ~ stat_post$batch)
    f_post <- summary(fit_post)[[1]][["F value"]][1]
    
    f_stats_df$F_Stat_Pre[i]  <- round(f_pre, 2)
    f_stats_df$F_Stat_Post[i] <- round(f_post, 2)
  }
  
  # Calculate Percentage Drop (Higher is better)
  f_stats_df <- f_stats_df %>%
    mutate(Percent_Reduction = round((1 - (F_Stat_Post / F_Stat_Pre)) * 100, 2)) %>%
    arrange(desc(Percent_Reduction))
  
  # ============================================================================
  # METRIC 2: BIOLOGICAL PRESERVATION (Sample-Level Spearman Correlation)
  # ============================================================================
  message("[Evaluation] Computing Biological Preservation (Spearman correlation of medians)...")
  
  # Calculate sample-level medians
  medians_pre <- df_pre %>%
    group_by(sample_id) %>%
    summarise(across(all_of(markers), median, .names = "{.col}")) %>%
    arrange(sample_id)
  
  medians_post <- df_post %>%
    group_by(sample_id) %>%
    summarise(across(all_of(markers), median, .names = "{.col}")) %>%
    arrange(sample_id)
  
  preservation_df <- data.frame(Marker = markers, Spearman_Rho = NA)
  
  for (i in seq_along(markers)) {
    m <- markers[i]
    # Suppress warnings for zero variance markers in specific edge cases
    rho <- suppressWarnings(cor(medians_pre[[m]], medians_post[[m]], method = "spearman"))
    preservation_df$Spearman_Rho[i] <- round(rho, 4)
  }
  
  preservation_df <- preservation_df %>% arrange(Spearman_Rho)
  
  # ============================================================================
  # 3. GLOBAL SUMMARY GENERATION
  # ============================================================================
  mean_reduction <- mean(f_stats_df$Percent_Reduction, na.rm = TRUE)
  mean_preservation <- mean(preservation_df$Spearman_Rho, na.rm = TRUE)
  
  summary_df <- data.frame(
    Metric = c("Mean Batch F-Stat Reduction (%)", "Mean Bio Preservation (Spearman)"),
    Value = c(round(mean_reduction, 2), round(mean_preservation, 4)),
    Target = c("> 80.00%", "> 0.8500"),
    Status = c(
      ifelse(mean_reduction > 80, "PASS", "WARNING"),
      ifelse(mean_preservation > 0.85, "PASS", "WARNING")
    )
  )
  
  message(sprintf("[Evaluation] Finished. Mean batch reduction: %.2f%% | Mean correlation: %.4f", 
                  mean_reduction, mean_preservation))
  
  return(list(
    Summary = summary_df,
    Batch_Removal = f_stats_df,
    Biological_Preservation = preservation_df
  ))
}