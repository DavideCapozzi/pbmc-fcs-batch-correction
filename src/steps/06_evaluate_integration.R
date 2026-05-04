#!/usr/bin/env Rscript
# src/steps/06_evaluate_integration.R
# ==============================================================================
# STEP 06: INTEGRATION VALIDATION — INTRA-BATCH STABILITY METRICS
# Description: Validates the stratified ecological baseline with two metrics:
#   (1) JSD between matched cluster centroid profiles — validates centroid
#       matching quality from Step 03.
#   (2) Intra-Batch Coefficient of Variation (CV = SD / Mean) — flags
#       populations with extreme within-batch frequency variability, indicating
#       reference intervals that may be unreliable for patient Z-scoring.
# Dependencies: writexl, yaml, dplyr
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(writexl)
})

source("src/functions/cluster_matching.R")

message("\n=== PIPELINE STEP 6: INTEGRATION VALIDATION ===")

config      <- read_yaml("config/global_params.yml")
out_dir     <- config$directories$processed
int_out_dir <- config$directories$integration %||% "results/integration"
if (!dir.exists(int_out_dir)) dir.create(int_out_dir, recursive = TRUE)

required_files <- c(
  baseline_file  = file.path(out_dir, "stratified_baseline_dictionary.rds"),
  match_file     = file.path(out_dir, "cluster_match_table.rds"),
  centroids_file = file.path(out_dir, "cluster_centroids.rds")
)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("[Step 06] Missing files: ", paste(names(missing), collapse = ", "),
       ". Run Steps 03-05 first.")
}

tryCatch({

  baseline_dict  <- readRDS(required_files["baseline_file"])
  match_table    <- readRDS(required_files["match_file"])
  centroids_list <- readRDS(required_files["centroids_file"])

  batches <- names(centroids_list)
  markers <- colnames(centroids_list[[batches[1L]]])

  # --------------------------------------------------------------------------
  # METRIC 1: JSD between matched cluster centroid profiles
  # --------------------------------------------------------------------------
  message("[Step 06] Computing Jensen-Shannon Divergence per matched pair...")
  validation <- validate_cluster_matches(match_table, centroids_list, markers)

  cluster_match_report <- validation$summary
  if (nrow(cluster_match_report) > 0) {
    n_pass <- sum(cluster_match_report$jsd_pass, na.rm = TRUE)
    message(sprintf("[Step 06] JSD < 0.1: %d / %d matched populations",
                    n_pass, nrow(cluster_match_report)))
  }

  # --------------------------------------------------------------------------
  # METRIC 2: Intra-Batch Coefficient of Variation (CV)
  # --------------------------------------------------------------------------
  message("[Step 06] Computing intra-batch stability (CV per population per batch)...")

  cv_threshold <- as.numeric(config$integration$cv_flag_threshold %||% 1.5)

  cv_table         <- baseline_dict
  cv_table$cv      <- cv_table$sd_freq / cv_table$mean_freq
  cv_table$cv_flag <- !is.na(cv_table$cv) & cv_table$cv > cv_threshold

  zero_mean_n <- sum(!is.finite(cv_table$cv), na.rm = TRUE)
  if (zero_mean_n > 0) {
    message(sprintf("[Step 06] %d strata with mean_freq = 0 (CV undefined): excluded from CV flagging.",
                    zero_mean_n))
  }
  cv_table$cv[!is.finite(cv_table$cv)]  <- NA_real_
  cv_table$cv_flag[is.na(cv_table$cv)]  <- NA

  n_flagged <- sum(cv_table$cv_flag, na.rm = TRUE)
  message(sprintf("[Step 06] Strata with CV > %.1f (high intra-batch variability): %d / %d",
                  cv_threshold, n_flagged, sum(!is.na(cv_table$cv_flag))))
  if (n_flagged > 0) {
    flagged_ids <- cv_table[!is.na(cv_table$cv_flag) & cv_table$cv_flag,
                             c("batch", "population", "cv")]
    for (i in seq_len(nrow(flagged_ids))) {
      message(sprintf("           %s / %s: CV = %.3f",
                      flagged_ids$batch[i], flagged_ids$population[i], flagged_ids$cv[i]))
    }
  }

  # --------------------------------------------------------------------------
  # EXPORT
  # --------------------------------------------------------------------------
  writexl::write_xlsx(
    list(
      ClusterMatch        = cluster_match_report,
      IntraBatchStability = cv_table
    ),
    path = file.path(int_out_dir, "integration_validation_report.xlsx")
  )

  message(sprintf("[Step 06] Validation report saved to: %s",
                  file.path(int_out_dir, "integration_validation_report.xlsx")))

  message("\n=== STEP 6 COMPLETE ===\n")

}, error = function(e) {
  message("\n[Error] Step 06 failed:")
  message(e$message)
  stop(paste("[Step 06 Fatal]", e$message), call. = FALSE)
})
