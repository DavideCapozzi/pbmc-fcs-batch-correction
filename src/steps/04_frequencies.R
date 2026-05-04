#!/usr/bin/env Rscript
# src/steps/04_frequencies.R
# ==============================================================================
# STEP 04: PER-SAMPLE POPULATION FREQUENCY COMPUTATION
# Computes relative cell frequencies for each matched metacluster population.
# Denominator = all cells (matched + unmatched) to preserve true proportions.
# Output: frequency_matrix.rds, frequency_matrix.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
  library(dplyr)
})

source("src/functions/frequency.R")
source("src/functions/step_logger.R")

message("\n=== PIPELINE STEP 4: POPULATION FREQUENCY COMPUTATION ===")

config  <- read_yaml("config/global_params.yml")
out_dir <- config$directories$processed
batches <- names(config$directories$raw)

match_file <- file.path(out_dir, "cluster_match_table.rds")
sce_files  <- file.path(out_dir, paste0("sce_batch_", batches, ".rds"))
missing    <- sce_files[!file.exists(sce_files)]

if (!file.exists(match_file)) {
  stop("[Step 04] cluster_match_table.rds not found. Run Step 03 first.")
}
if (length(missing) > 0) {
  stop("[Step 04] Missing per-batch SCE file(s): ", paste(missing, collapse = ", "))
}

log_obj <- init_step_log(
  step_name   = "04_frequencies",
  step_number = 4L,
  input_files = c(match_file, sce_files)
)

tryCatch({

  match_table <- readRDS(match_file)
  sce_list    <- lapply(setNames(sce_files, batches), readRDS)

  n_matched_pops <- sum(!is.na(match_table$matched_population_id) & match_table$is_matched)
  message(sprintf("[Step 04] Matched populations to quantify: %d", n_matched_pops))

  freq_result <- compute_population_frequencies(sce_list, match_table)

  fm       <- freq_result$freq_matrix
  cm       <- freq_result$count_matrix
  pop_cols <- grep("^pop_", colnames(cm), value = TRUE)

  matched_cells <- rowSums(cm[, pop_cols, drop = FALSE])
  pct_assigned  <- mean(matched_cells / cm$total_cells, na.rm = TRUE) * 100

  message(sprintf("[Step 04] Samples processed: %d", nrow(fm)))
  message(sprintf("[Step 04] Batches: %s", paste(unique(fm$batch), collapse = ", ")))
  message(sprintf("[Step 04] Mean fraction of cells in matched populations: %.1f%%", pct_assigned))

  if (pct_assigned < 50) {
    warning("[Step 04] Less than 50% of cells assigned to matched populations. ",
            "Consider lowering integration.similarity_threshold.", call. = FALSE)
  }

  out_rds <- file.path(out_dir, "frequency_matrix.rds")
  out_csv <- file.path(out_dir, "frequency_matrix.csv")
  saveRDS(freq_result, out_rds)
  write.csv(fm, out_csv, row.names = FALSE)

  message(sprintf("[Step 04] Frequency matrix saved: %d samples x %d populations",
                  nrow(fm), length(pop_cols)))

  add_metric(log_obj, "n_samples",         nrow(fm))
  add_metric(log_obj, "n_populations",     length(pop_cols))
  add_metric(log_obj, "mean_pct_assigned", round(pct_assigned, 2))
  add_metric(log_obj, "batches",           unique(fm$batch))

  finalize_step_log(log_obj, output_files = c(out_rds, out_csv), status = "SUCCESS")
  write_step_json(log_obj, config$directories$logs)

  message("\n=== STEP 4 COMPLETE ===\n")

}, error = function(e) {
  add_metric(log_obj, "error_message", e$message)
  finalize_step_log(log_obj, output_files = character(0L), status = "FAILURE")
  write_step_json(log_obj, config$directories$logs)
  message("\n[Error] Step 04 failed: ", e$message)
  stop(paste("[Step 04 Fatal]", e$message), call. = FALSE)
})
