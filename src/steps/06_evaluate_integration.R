#!/usr/bin/env Rscript
# src/steps/06_evaluate_integration.R
# ==============================================================================
# STEP 06: INTEGRATION VALIDATION — INTRA-BATCH STABILITY METRICS
# (1) JSD between matched cluster centroid profiles — validates matching quality.
# (2) Intra-Batch CV (SD/Mean) — flags populations with high within-batch
#     frequency variability that would produce unreliable Z-score baselines.
# Output: integration_validation_report.xlsx (ClusterMatch + IntraBatchStability)
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(writexl)
})

source("src/functions/utils.R")
source("src/functions/cluster_matching.R")
source("src/functions/step_logger.R")

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

log_obj <- init_step_log(
  step_name   = "06_evaluate_integration",
  step_number = 6L,
  input_files = unname(required_files)
)

tryCatch({

  baseline_dict  <- readRDS(required_files["baseline_file"])
  match_table    <- readRDS(required_files["match_file"])
  centroids_list <- readRDS(required_files["centroids_file"])

  batches <- names(centroids_list)
  markers <- colnames(centroids_list[[batches[1L]]])

  # Metric 1: Per-marker delta (primary) + JSD (secondary) per matched pair
  delta_threshold <- as.numeric(config$integration$max_marker_delta_arcsinh %||% 2.0)
  jsd_threshold   <- as.numeric(config$integration$jsd_pass_threshold        %||% 0.1)
  message("[Step 06] Computing per-marker delta and JSD per matched pair...")
  validation           <- validate_cluster_matches(match_table, centroids_list, markers,
                                                    delta_threshold = delta_threshold,
                                                    jsd_threshold   = jsd_threshold)
  cluster_match_report <- validation$summary

  n_jsd_pass        <- 0L
  n_delta_pass      <- 0L
  mean_delta_global <- NA_real_
  if (nrow(cluster_match_report) > 0L) {
    n_jsd_pass        <- sum(cluster_match_report$jsd_pass,  na.rm = TRUE)
    n_delta_pass      <- sum(cluster_match_report$delta_pass, na.rm = TRUE)
    mean_delta_global <- round(mean(cluster_match_report$mean_marker_delta, na.rm = TRUE), 4L)
    message(sprintf("[Step 06] JSD < 0.1: %d / %d matched populations",
                    n_jsd_pass, nrow(cluster_match_report)))
    message(sprintf("[Step 06] MarkerDelta < %.1f arcsinh: %d / %d matched populations (mean=%.3f)",
                    delta_threshold, n_delta_pass, nrow(cluster_match_report), mean_delta_global))
  }

  # Metric 2: Intra-batch CV
  message("[Step 06] Computing intra-batch CV per (batch, population)...")
  cv_threshold <- as.numeric(config$integration$cv_flag_threshold %||% 1.5)

  cv_table         <- baseline_dict
  cv_table$cv      <- safe_div(cv_table$sd_freq, cv_table$mean_freq, default = NA_real_)
  cv_table$cv_flag <- !is.na(cv_table$cv) & cv_table$cv > cv_threshold

  zero_mean_n <- sum(is.na(cv_table$cv), na.rm = TRUE)
  if (zero_mean_n > 0) {
    message(sprintf("[Step 06] %d strata with mean_freq ~ 0 (CV undefined): set to NA.",
                    zero_mean_n))
  }
  cv_table$cv_flag[is.na(cv_table$cv)] <- NA

  n_flagged <- sum(cv_table$cv_flag, na.rm = TRUE)
  message(sprintf("[Step 06] Strata with CV > %.1f: %d / %d",
                  cv_threshold, n_flagged, sum(!is.na(cv_table$cv_flag))))

  if (n_flagged > 0) {
    flagged_ids <- cv_table[!is.na(cv_table$cv_flag) & cv_table$cv_flag,
                             c("batch", "population", "cv")]
    for (i in seq_len(nrow(flagged_ids))) {
      message(sprintf("           %s / %s: CV = %.3f",
                      flagged_ids$batch[i], flagged_ids$population[i], flagged_ids$cv[i]))
    }
  }

  out_xlsx <- file.path(int_out_dir, "integration_validation_report.xlsx")
  writexl::write_xlsx(
    list(ClusterMatch = cluster_match_report, IntraBatchStability = cv_table),
    path = out_xlsx
  )
  message(sprintf("[Step 06] Validation report saved: %s", out_xlsx))

  add_metric(log_obj, "n_matched_pairs",       nrow(cluster_match_report))
  add_metric(log_obj, "n_jsd_pass",            n_jsd_pass)
  add_metric(log_obj, "n_delta_pass",          n_delta_pass)
  add_metric(log_obj, "mean_marker_delta_global", mean_delta_global)
  add_metric(log_obj, "delta_threshold",       delta_threshold)
  add_metric(log_obj, "n_cv_flagged",          n_flagged)
  add_metric(log_obj, "cv_threshold",          cv_threshold)

  finalize_step_log(log_obj, output_files = out_xlsx, status = "SUCCESS")
  write_step_json(log_obj, config$directories$logs)

  message("\n=== STEP 6 COMPLETE ===\n")

}, error = function(e) {
  add_metric(log_obj, "error_message", e$message)
  finalize_step_log(log_obj, output_files = character(0L), status = "FAILURE")
  write_step_json(log_obj, config$directories$logs)
  message("\n[Error] Step 06 failed: ", e$message)
  stop(paste("[Step 06 Fatal]", e$message), call. = FALSE)
})
