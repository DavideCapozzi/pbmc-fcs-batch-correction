#!/usr/bin/env Rscript
# src/steps/05_integrate.R
# ==============================================================================
# STEP 05: STRATIFIED ECOLOGICAL BASELINE
# Description: Computes independent reference distributions (Mean, SD, 95%
#              non-parametric Reference Intervals) for each matched population,
#              stratified by batch. Batches are treated as isolated, independent
#              experimental conditions — no cross-batch pooling is performed.
#
# Outputs:
#   stratified_baseline_dictionary.rds  — reference dict (consumed by Step 06)
#   stratified_baseline_report.xlsx     — two sheets:
#     StratifiedBaseline : mean/SD/RI per (batch, population), optional labels
#     Centroids          : marker centroids per population per batch (for
#                          biological annotation — fill population_labels in config)
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(tidyr)
  library(dplyr)
  library(writexl)
})

source("src/functions/integration.R")

message("\n=== PIPELINE STEP 5: STRATIFIED ECOLOGICAL BASELINE ===")

config      <- read_yaml("config/global_params.yml")
int_cfg     <- config$integration
out_dir     <- config$directories$processed
int_out_dir <- config$directories$integration

if (is.null(int_out_dir)) int_out_dir <- "results/integration"
if (!dir.exists(int_out_dir)) dir.create(int_out_dir, recursive = TRUE)

required_files <- c(
  freq_file      = file.path(out_dir, "frequency_matrix.rds"),
  match_file     = file.path(out_dir, "cluster_match_table.rds"),
  centroids_file = file.path(out_dir, "cluster_centroids.rds")
)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("[Step 05] Missing files: ", paste(names(missing), collapse = ", "),
       ". Run Steps 03-04 first.")
}

tryCatch({

  freq_result    <- readRDS(required_files["freq_file"])
  match_table    <- readRDS(required_files["match_file"])
  centroids_list <- readRDS(required_files["centroids_file"])
  fm             <- freq_result$freq_matrix

  # Minimum samples guard — SD and non-parametric RI require replication
  min_samp <- if (isTRUE(config$testing$enabled)) 2L else
              as.integer(int_cfg$min_samples_per_batch %||% 3L)

  n_per_batch <- table(fm$batch)
  message("[Step 05] Samples per batch:")
  for (bat in names(n_per_batch)) {
    message(sprintf("           %s: %d", bat, n_per_batch[[bat]]))
  }

  if (any(n_per_batch < min_samp)) {
    stop(sprintf("[Step 05] Insufficient samples. Require >= %d per batch (current min: %d).",
                 min_samp, min(n_per_batch)))
  }

  message("[Step 05] Computing stratified reference distributions...")

  baseline_dict <- compute_stratified_baseline(fm)

  # Optional biological annotation from config:
  # fill population_labels in global_params.yml after reviewing the
  # Centroids sheet in stratified_baseline_report.xlsx
  pop_labels    <- config$population_labels
  if (!is.null(pop_labels)) {
    message(sprintf("[Step 05] Applying %d population labels from config.", length(pop_labels)))
  }

  baseline_report <- format_stratified_report(baseline_dict, match_table, pop_labels)

  for (bat in sort(unique(baseline_dict$batch))) {
    bd_bat <- baseline_dict[baseline_dict$batch == bat, ]
    message(sprintf("[Step 05] [%s] %d populations | mean freq range: [%.3f, %.3f]",
                    bat,
                    nrow(bd_bat),
                    min(bd_bat$mean_freq, na.rm = TRUE),
                    max(bd_bat$mean_freq, na.rm = TRUE)))
  }

  markers         <- colnames(centroids_list[[names(centroids_list)[1L]]])
  centroid_report <- format_centroid_report(centroids_list, match_table, markers)

  saveRDS(baseline_dict, file.path(out_dir, "stratified_baseline_dictionary.rds"))

  writexl::write_xlsx(
    list(
      StratifiedBaseline = baseline_report,
      Centroids          = centroid_report
    ),
    path = file.path(int_out_dir, "stratified_baseline_report.xlsx")
  )

  message(sprintf("[Step 05] Baseline saved: %d strata (%d batches x %d populations)",
                  nrow(baseline_dict),
                  length(unique(baseline_dict$batch)),
                  length(unique(baseline_dict$population))))
  message(sprintf("[Step 05] Centroid sheet: %d rows (%d populations x 2 batches)",
                  nrow(centroid_report),
                  length(unique(centroid_report$population))))

  message("\n=== STEP 5 COMPLETE ===\n")

}, error = function(e) {
  message("\n[Error] Step 05 failed:")
  message(e$message)
  stop(paste("[Step 05 Fatal]", e$message), call. = FALSE)
})
