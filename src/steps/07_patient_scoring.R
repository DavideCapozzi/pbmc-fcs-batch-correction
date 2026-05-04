#!/usr/bin/env Rscript
# src/steps/07_patient_scoring.R
# ==============================================================================
# STEP 07: CLINICAL INFERENCE (Z-SCORE)
# Description: Evaluates patient samples against the stratified baseline dictionary
#              and generates a clinical flagging report.
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(writexl)
})

source("src/functions/scoring.R")

message("\n=== PIPELINE STEP 7: CLINICAL PATIENT SCORING ===")

config      <- read_yaml("config/global_params.yml")
out_dir     <- config$directories$processed
int_out_dir <- config$directories$integration %||% "results/integration"

# Check if scoring is enabled
if (!isTRUE(config$scoring$run_patient_scoring)) {
  message("[Step 07] scoring.run_patient_scoring is FALSE. Skipping inference.")
  quit(save = "no", status = 0)
}

req_files <- c(
  baseline = file.path(out_dir, "stratified_baseline_dictionary.rds"),
  match    = file.path(out_dir, "cluster_match_table.rds"),
  # NOTE: Assuming patient frequencies are generated and saved here 
  # (in a real scenario, this might come from a distinct patient ingestion path)
  patients = file.path(out_dir, "frequency_matrix.rds") 
)

missing <- req_files[!file.exists(req_files)]
if (length(missing) > 0) {
  stop("[Step 07 Fatal] Missing dependency files: ", paste(names(missing), collapse = ", "))
}

tryCatch({
  baseline_dict <- readRDS(req_files["baseline"])
  match_table   <- readRDS(req_files["match"])
  # For the initial implementation, we score the existing frequency matrix
  # Later, this points to patient_frequency_matrix.rds
  patient_freq  <- readRDS(req_files["patients"])$freq_matrix 
  
  z_thresh <- as.numeric(config$scoring$z_score_threshold %||% 2.0)
  message(sprintf("[Step 07] Initiating clinical inference (Z-Threshold: +/- %.1f)...", z_thresh))
  
  scored_data <- compute_patient_zscores(patient_freq, baseline_dict, z_threshold = z_thresh)
  
  pop_labels <- config$population_labels
  clinical_report <- format_clinical_report(scored_data, match_table, pop_labels)
  
  # Summary diagnostics
  expanded_count <- sum(clinical_report$clinical_flag == "EXPANDED", na.rm = TRUE)
  depleted_count <- sum(clinical_report$clinical_flag == "DEPLETED", na.rm = TRUE)
  
  message(sprintf("[Step 07] Scoring completed for %d patient/sample strata.", nrow(clinical_report)))
  message(sprintf("[Step 07] Detected %d Expanded events and %d Depleted events.", expanded_count, depleted_count))
  
  # Export Output
  saveRDS(clinical_report, file.path(out_dir, "patient_zscores.rds"))
  writexl::write_xlsx(
    list(ClinicalScoring = clinical_report),
    path = file.path(int_out_dir, "clinical_patient_report.xlsx")
  )
  
  message(sprintf("[Step 07] Clinical report exported to: %s", file.path(int_out_dir, "clinical_patient_report.xlsx")))
  message("\n=== STEP 7 COMPLETE ===\n")
  
}, error = function(e) {
  stop(paste("\n[Step 07 Fatal]", e$message), call. = FALSE)
})