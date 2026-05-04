# src/legacy/functions/evaluation.R
# ==============================================================================
# [LEGACY] EVALUATION MODULE — MFI CORRECTION QUALITY METRICS
# Description: ANOVA F-statistics for batch removal and Spearman correlation
#              for biological preservation. Used only in the legacy MFI
#              correction path (run_per_batch = false).
# Dependencies: SingleCellExperiment, dplyr, tidyr, stats
# ==============================================================================

library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(stats)

#' @title Evaluate Batch Correction
#' @description Computes ANOVA F-statistics for batch removal and Spearman
#'   correlation for biological signal preservation.
#' @param sce_pre Uncorrected SCE object.
#' @param sce_post Corrected SCE object.
#' @param markers Character vector of markers to evaluate.
#' @return Named list: $Summary, $Batch_Removal, $Biological_Preservation.
evaluate_batch_correction <- function(sce_pre, sce_post, markers) {

  message("[Evaluation] Initializing metric calculation...")

  if (ncol(sce_pre) != ncol(sce_post)) {
    stop("[Evaluation Fatal] Pre and Post SCE objects have different cell counts.")
  }

  meta_df    <- as.data.frame(colData(sce_pre))
  exprs_pre  <- t(assay(sce_pre,  "exprs"))
  exprs_post <- t(assay(sce_post, "exprs"))

  df_pre  <- dplyr::bind_cols(
    as.data.frame(exprs_pre[,  markers, drop = FALSE], check.names = FALSE), meta_df)
  df_post <- dplyr::bind_cols(
    as.data.frame(exprs_post[, markers, drop = FALSE], check.names = FALSE), meta_df)

  # --------------------------------------------------------------------------
  # METRIC 1: Batch effect removal — ANOVA F-statistic
  # --------------------------------------------------------------------------
  message("[Evaluation] Computing ANOVA F-Statistics for batch variance...")

  n_cells <- nrow(df_pre)
  if (n_cells > 100000L) {
    message(sprintf("[Evaluation] Subsampling 100000 / %d cells for F-stat computation...", n_cells))
    set.seed(42L)
    idx      <- sample(seq_len(n_cells), 100000L)
    stat_pre  <- df_pre[idx, ]
    stat_post <- df_post[idx, ]
  } else {
    stat_pre  <- df_pre
    stat_post <- df_post
  }

  f_stats_df <- data.frame(Marker = markers, F_Stat_Pre = NA_real_, F_Stat_Post = NA_real_)
  for (i in seq_along(markers)) {
    m <- markers[i]
    f_stats_df$F_Stat_Pre[i]  <- summary(aov(stat_pre[[m]]  ~ stat_pre$batch))[[1]][["F value"]][1L]
    f_stats_df$F_Stat_Post[i] <- summary(aov(stat_post[[m]] ~ stat_post$batch))[[1]][["F value"]][1L]
  }
  f_stats_df <- f_stats_df %>%
    dplyr::mutate(Percent_Reduction = round((1 - F_Stat_Post / F_Stat_Pre) * 100, 2)) %>%
    dplyr::arrange(dplyr::desc(Percent_Reduction))

  # --------------------------------------------------------------------------
  # METRIC 2: Biological preservation — Spearman ρ of sample medians
  # --------------------------------------------------------------------------
  message("[Evaluation] Computing Biological Preservation (Spearman ρ of sample medians)...")

  medians_pre <- df_pre %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(markers), stats::median), .groups = "drop") %>%
    dplyr::arrange(sample_id)

  medians_post <- df_post %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(markers), stats::median), .groups = "drop") %>%
    dplyr::arrange(sample_id)

  preservation_df <- data.frame(Marker = markers, Spearman_Rho = NA_real_)
  for (i in seq_along(markers)) {
    m <- markers[i]
    preservation_df$Spearman_Rho[i] <- round(
      suppressWarnings(stats::cor(medians_pre[[m]], medians_post[[m]], method = "spearman")), 4L)
  }
  preservation_df <- dplyr::arrange(preservation_df, Spearman_Rho)

  # --------------------------------------------------------------------------
  # Summary
  # --------------------------------------------------------------------------
  mean_reduction   <- mean(f_stats_df$Percent_Reduction,  na.rm = TRUE)
  mean_preservation <- mean(preservation_df$Spearman_Rho, na.rm = TRUE)

  summary_df <- data.frame(
    Metric = c("Mean Batch F-Stat Reduction (%)", "Mean Bio Preservation (Spearman)"),
    Value  = c(round(mean_reduction, 2), round(mean_preservation, 4)),
    Target = c("> 80.00%", "> 0.8500"),
    Status = c(
      ifelse(mean_reduction    > 80,   "PASS", "WARNING"),
      ifelse(mean_preservation > 0.85, "PASS", "WARNING")
    )
  )

  message(sprintf("[Evaluation] Mean batch reduction: %.2f%% | Mean Spearman ρ: %.4f",
                  mean_reduction, mean_preservation))

  list(Summary = summary_df, Batch_Removal = f_stats_df, Biological_Preservation = preservation_df)
}
