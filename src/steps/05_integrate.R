#!/usr/bin/env Rscript
# src/steps/05_integrate.R
# ==============================================================================
# STEP 05: STRATIFIED ECOLOGICAL BASELINE
# Computes independent reference distributions (Mean, SD, 95% non-parametric
# Reference Intervals) for each matched population, stratified by batch.
# Batches are treated as isolated experimental conditions — no cross-batch
# pooling is performed.
#
# Output:
#   stratified_baseline_dictionary.rds
#   stratified_baseline_report.xlsx (StratifiedBaseline + Centroids sheets)
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(tidyr)
  library(dplyr)
  library(writexl)
})

source("src/functions/integration.R")
source("src/functions/step_logger.R")

message("\n=== PIPELINE STEP 5: STRATIFIED ECOLOGICAL BASELINE ===")

config      <- read_yaml("config/global_params.yml")
int_cfg     <- config$integration
out_dir     <- config$directories$processed
int_out_dir <- config$directories$integration %||% "results/integration"

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

log_obj <- init_step_log(
  step_name   = "05_integrate",
  step_number = 5L,
  input_files = unname(required_files)
)

tryCatch({

  freq_result    <- readRDS(required_files["freq_file"])
  match_table    <- readRDS(required_files["match_file"])
  centroids_list <- readRDS(required_files["centroids_file"])

  fm_check <- freq_result$freq_matrix
  if (is.null(fm_check) || nrow(fm_check) == 0L) {
    stop("[Step 05] frequency_matrix.rds is empty. Run Step 04 first.")
  }
  if (!all(c("sample_id", "batch") %in% colnames(fm_check))) {
    stop("[Step 05] frequency_matrix is missing required columns: sample_id, batch.")
  }

  # --- AUTO-PHENOTYPING ---
  source("src/functions/auto_phenotyping.R")
  message("[Step 05] Running auto-phenotyping on matched population centroids...")
  pheno_result <- tryCatch(
    run_auto_phenotyping(
      centroids_list            = centroids_list,
      match_table               = match_table,
      markers                   = colnames(centroids_list[[names(centroids_list)[1L]]]),
      uninformative_gap         = as.numeric(config$auto_phenotyping$uninformative_gap %||% 0.5),
      low_reliability_threshold = as.numeric(config$auto_phenotyping$low_reliability_threshold %||% 0.30),
      kmeans_nstart             = as.integer(config$auto_phenotyping$kmeans$nstart  %||% 10L),
      kmeans_iter_max           = as.integer(config$auto_phenotyping$kmeans$iter_max %||% 50L),
      marker_min_thresholds     = config$auto_phenotyping$marker_min_thresholds %||% list()
    ),
    error = function(e) {
      warning(sprintf("[Step 05] Auto-phenotyping failed: %s — continuing without labels.",
                      e$message), call. = FALSE)
      NULL
    }
  )
  pheno_xlsx <- NULL
  if (!is.null(pheno_result)) {
    pheno_xlsx <- file.path(int_out_dir, "phenotype_annotation_report.xlsx")
    writexl::write_xlsx(list(PhenotypeAnnotation = pheno_result$report), path = pheno_xlsx)
    message(sprintf("[Step 05] Phenotype annotation saved: %s  (%d rows)",
                    pheno_xlsx, nrow(pheno_result$report)))
    saveRDS(as.list(pheno_result$population_labels),
            file.path(out_dir, "auto_population_labels.rds"))
  }
  # --- END AUTO-PHENOTYPING ---

  fm             <- freq_result$freq_matrix

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

  # Flag populations with suspiciously high frequency — may indicate clustering collapse
  max_freq_warn  <- as.numeric(int_cfg$max_population_freq_warning %||% 0.40)
  dominant_flags <- baseline_dict[baseline_dict$mean_freq > max_freq_warn, ]
  if (nrow(dominant_flags) > 0L) {
    for (i in seq_len(nrow(dominant_flags))) {
      warning(sprintf(
        "[Step 05] Dominant population: %s in batch '%s' (mean_freq=%.3f > %.2f) — verify in UMAP",
        dominant_flags$population[i], dominant_flags$batch[i],
        dominant_flags$mean_freq[i], max_freq_warn
      ), call. = FALSE)
    }
    add_metric(log_obj, "dominant_populations_flagged",
               lapply(seq_len(nrow(dominant_flags)), function(i) list(
                 batch      = dominant_flags$batch[i],
                 population = dominant_flags$population[i],
                 mean_freq  = round(dominant_flags$mean_freq[i], 4L)
               )))
  }
  baseline_dict$freq_warning <- baseline_dict$mean_freq > max_freq_warn

  pop_labels <- config$population_labels
  if (!is.null(pop_labels) && length(pop_labels) > 0L) {
    message(sprintf("[Step 05] Applying %d population labels from config.", length(pop_labels)))
  } else if (!is.null(pheno_result)) {
    pop_labels <- as.list(pheno_result$population_labels)
    message(sprintf("[Step 05] Using %d auto-generated population labels.", length(pop_labels)))
  }

  baseline_report <- format_stratified_report(baseline_dict, match_table, pop_labels)

  freq_ranges <- list()
  for (bat in sort(unique(baseline_dict$batch))) {
    bd_bat <- baseline_dict[baseline_dict$batch == bat, ]
    min_f  <- round(min(bd_bat$mean_freq, na.rm = TRUE), 4)
    max_f  <- round(max(bd_bat$mean_freq, na.rm = TRUE), 4)
    message(sprintf("[Step 05] [%s] %d populations | mean freq range: [%.3f, %.3f]",
                    bat, nrow(bd_bat), min_f, max_f))
    freq_ranges[[bat]] <- list(min = min_f, max = max_f)
  }

  markers         <- colnames(centroids_list[[names(centroids_list)[1L]]])
  centroid_report <- format_centroid_report(centroids_list, match_table, markers)

  out_rds  <- file.path(out_dir, "stratified_baseline_dictionary.rds")
  out_xlsx <- file.path(int_out_dir, "stratified_baseline_report.xlsx")

  saveRDS(baseline_dict, out_rds)
  writexl::write_xlsx(
    list(StratifiedBaseline = baseline_report, Centroids = centroid_report),
    path = out_xlsx
  )

  n_strata <- nrow(baseline_dict)
  message(sprintf("[Step 05] Baseline saved: %d strata (%d batches x %d populations)",
                  n_strata,
                  length(unique(baseline_dict$batch)),
                  length(unique(baseline_dict$population))))

  add_metric(log_obj, "n_strata",     n_strata)
  add_metric(log_obj, "n_batches",    length(unique(baseline_dict$batch)))
  add_metric(log_obj, "n_populations",length(unique(baseline_dict$population)))
  add_metric(log_obj, "freq_ranges_per_batch", freq_ranges)
  if (!is.null(pheno_result)) {
    add_metric(log_obj, "n_auto_phenotyped", nrow(pheno_result$report))
    add_metric(log_obj, "phenotype_labels",  as.list(pheno_result$population_labels))
  }

  pheno_out <- if (!is.null(pheno_xlsx)) pheno_xlsx else character(0L)
  finalize_step_log(log_obj, output_files = c(out_rds, out_xlsx, pheno_out), status = "SUCCESS")
  write_step_json(log_obj, config$directories$logs)

  message("\n=== STEP 5 COMPLETE ===\n")

}, error = function(e) {
  add_metric(log_obj, "error_message", e$message)
  finalize_step_log(log_obj, output_files = character(0L), status = "FAILURE")
  write_step_json(log_obj, config$directories$logs)
  message("\n[Error] Step 05 failed: ", e$message)
  stop(paste("[Step 05 Fatal]", e$message), call. = FALSE)
})
