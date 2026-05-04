#!/usr/bin/env Rscript
# src/steps/05_integrate.R
# ==============================================================================
# STEP 05: DERSIMONIANLAIRD RANDOM EFFECTS META-ANALYSIS — HEALTHY BASELINE
# Description: Pools per-sample population frequencies across cohorts using
#              the DerSimonian-Laird estimator after arcsine-sqrt transformation.
#              Produces the pooled healthy immune baseline with 95% CIs,
#              heterogeneity statistics, and bootstrap validation.
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(writexl)
})

source("src/functions/integration.R")

message("\n=== PIPELINE STEP 5: RANDOM EFFECTS INTEGRATION (HEALTHY BASELINE) ===")

config      <- read_yaml("config/global_params.yml")
int_cfg     <- config$integration
out_dir     <- config$directories$processed
int_out_dir <- config$directories$integration

if (is.null(int_out_dir)) int_out_dir <- "results/integration"
if (!dir.exists(int_out_dir)) dir.create(int_out_dir, recursive = TRUE)

freq_file  <- file.path(out_dir, "frequency_matrix.rds")
match_file <- file.path(out_dir, "cluster_match_table.rds")

if (!file.exists(freq_file))  stop("[Step 05] frequency_matrix.rds not found. Run Step 04 first.")
if (!file.exists(match_file)) stop("[Step 05] cluster_match_table.rds not found. Run Step 03 first.")

tryCatch({

  freq_result <- readRDS(freq_file)
  match_table <- readRDS(match_file)
  fm          <- freq_result$freq_matrix
  cm          <- freq_result$count_matrix

  # Minimum samples guard (relaxed in testing mode)
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

  n_boot <- as.integer(int_cfg$n_bootstrap %||% 1000L)
  seed   <- as.integer(config$random_seed %||% 42L)

  message(sprintf("[Step 05] Running DL meta-analysis (%d bootstrap iterations)...", n_boot))

  integration_result <- run_random_effects_integration(fm, cm, n_boot, seed)

  baseline <- integration_result$baseline
  high_het  <- baseline[baseline$flag_heterogeneity, "population"]
  message(sprintf("[Step 05] Populations with I² > 75%% (high heterogeneity): %d",
                  length(high_het)))
  if (length(high_het) > 0) {
    message("           ", paste(high_het, collapse = ", "))
  }

  baseline_report <- format_baseline_report(integration_result, match_table)

  saveRDS(integration_result, file.path(out_dir, "healthy_baseline.rds"))

  writexl::write_xlsx(
    list(
      Baseline   = baseline_report,
      PerSample  = integration_result$per_sample,
      Bootstrap  = integration_result$bootstrap
    ),
    path = file.path(int_out_dir, "healthy_baseline_report.xlsx")
  )

  message(sprintf("[Step 05] Baseline saved: %d populations", nrow(baseline_report)))
  message(sprintf("[Step 05] Frequency range: [%.3f, %.3f]",
                  min(baseline_report$freq_pooled, na.rm = TRUE),
                  max(baseline_report$freq_pooled, na.rm = TRUE)))

  message("\n=== STEP 5 COMPLETE ===\n")

}, error = function(e) {
  message("\n[Error] Step 05 failed:")
  message(e$message)
  stop(paste("[Step 05 Fatal]", e$message), call. = FALSE)
})
