#!/usr/bin/env Rscript
# src/steps/07_patient_scoring.R
# ==============================================================================
# STEP 07: CLINICAL INFERENCE (Z-SCORE)
# Evaluates patient samples against the stratified baseline dictionary and
# generates a clinical flagging report (EXPANDED / DEPLETED / NORMAL).
# Output: patient_zscores.rds, clinical_patient_report.xlsx
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(writexl)
})

source("src/functions/scoring.R")
source("src/functions/step_logger.R")

message("\n=== PIPELINE STEP 7: CLINICAL PATIENT SCORING ===")

config      <- read_yaml("config/global_params.yml")
out_dir     <- config$directories$processed
int_out_dir <- config$directories$integration %||% "results/integration"

if (!isTRUE(config$scoring$run_patient_scoring)) {
  message("[Step 07] scoring.run_patient_scoring is FALSE — skipping.")
  quit(save = "no", status = 0)
}

req_files <- c(
  baseline = file.path(out_dir, "stratified_baseline_dictionary.rds"),
  match    = file.path(out_dir, "cluster_match_table.rds"),
  patients = file.path(out_dir, "frequency_matrix.rds")
)
missing <- req_files[!file.exists(req_files)]
if (length(missing) > 0) {
  stop("[Step 07 Fatal] Missing dependency files: ", paste(names(missing), collapse = ", "))
}

log_obj <- init_step_log(
  step_name   = "07_patient_scoring",
  step_number = 7L,
  input_files = unname(req_files)
)

tryCatch({

  baseline_dict <- readRDS(req_files["baseline"])
  match_table   <- readRDS(req_files["match"])
  patient_freq  <- readRDS(req_files["patients"])$freq_matrix

  z_thresh <- as.numeric(config$scoring$z_score_threshold %||% 2.0)
  message(sprintf("[Step 07] Z-threshold: +/- %.1f", z_thresh))

  scored_data      <- compute_patient_zscores(patient_freq, baseline_dict, z_threshold = z_thresh)
  auto_labels_file <- file.path(out_dir, "auto_population_labels.rds")
  pop_labels <- if (!is.null(config$population_labels) && length(config$population_labels) > 0L) {
    config$population_labels
  } else if (file.exists(auto_labels_file)) {
    message("[Step 07] Loading auto-generated population labels from Step 05.")
    readRDS(auto_labels_file)
  } else {
    NULL
  }
  clinical_report <- format_clinical_report(scored_data, match_table, pop_labels)

  n_expanded     <- sum(clinical_report$clinical_flag == "EXPANDED",      na.rm = TRUE)
  n_depleted     <- sum(clinical_report$clinical_flag == "DEPLETED",      na.rm = TRUE)
  n_undetermined <- sum(clinical_report$clinical_flag == "UNDETERMINED",  na.rm = TRUE)

  message(sprintf("[Step 07] Scored %d strata — EXPANDED: %d, DEPLETED: %d, UNDETERMINED: %d",
                  nrow(clinical_report), n_expanded, n_depleted, n_undetermined))

  out_rds  <- file.path(out_dir, "patient_zscores.rds")
  out_xlsx <- file.path(int_out_dir, "clinical_patient_report.xlsx")

  saveRDS(clinical_report, out_rds)
  writexl::write_xlsx(list(ClinicalScoring = clinical_report), path = out_xlsx)
  message(sprintf("[Step 07] Clinical report exported: %s", out_xlsx))

  add_metric(log_obj, "n_strata_scored", nrow(clinical_report))
  add_metric(log_obj, "n_expanded",      n_expanded)
  add_metric(log_obj, "n_depleted",      n_depleted)
  add_metric(log_obj, "n_undetermined",  n_undetermined)
  add_metric(log_obj, "z_threshold",     z_thresh)

  finalize_step_log(log_obj, output_files = c(out_rds, out_xlsx), status = "SUCCESS")
  write_step_json(log_obj, config$directories$logs)

  message("\n=== STEP 7 COMPLETE ===\n")

}, error = function(e) {
  add_metric(log_obj, "error_message", e$message)
  finalize_step_log(log_obj, output_files = character(0L), status = "FAILURE")
  write_step_json(log_obj, config$directories$logs)
  stop(paste("[Step 07 Fatal]", e$message), call. = FALSE)
})
