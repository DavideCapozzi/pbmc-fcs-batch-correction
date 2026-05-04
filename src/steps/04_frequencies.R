#!/usr/bin/env Rscript
# src/steps/04_frequencies.R
# ==============================================================================
# STEP 04: PER-SAMPLE POPULATION FREQUENCY COMPUTATION
# Description: Computes relative cell frequencies for each matched metacluster
#              population across all samples. Produces frequency_matrix.rds
#              (proportions) and count_matrix (raw counts) for meta-analysis.
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
  library(dplyr)
})

source("src/functions/frequency.R")

message("\n=== PIPELINE STEP 4: POPULATION FREQUENCY COMPUTATION ===")

config  <- read_yaml("config/global_params.yml")
out_dir <- config$directories$processed

batches <- names(config$directories$raw)

match_file <- file.path(out_dir, "cluster_match_table.rds")
if (!file.exists(match_file)) {
  stop("[Step 04] cluster_match_table.rds not found. Run Step 03 (per-batch) first.")
}

sce_files <- file.path(out_dir, paste0("sce_batch_", batches, ".rds"))
missing   <- sce_files[!file.exists(sce_files)]
if (length(missing) > 0) {
  stop("[Step 04] Missing per-batch SCE file(s): ", paste(missing, collapse = ", "))
}

tryCatch({

  match_table <- readRDS(match_file)

  sce_list <- lapply(setNames(sce_files, batches), readRDS)

  n_matched_pops <- sum(!is.na(match_table$matched_population_id) & match_table$is_matched)
  message(sprintf("[Step 04] Matched populations to quantify: %d", n_matched_pops))

  freq_result <- compute_population_frequencies(sce_list, match_table)

  fm <- freq_result$freq_matrix
  cm <- freq_result$count_matrix

  message(sprintf("[Step 04] Samples processed: %d", nrow(fm)))
  message(sprintf("[Step 04] Batches: %s", paste(unique(fm$batch), collapse = ", ")))

  # Diagnostic: fraction of cells assigned to matched populations
  pop_cols     <- grep("^pop_", colnames(cm), value = TRUE)
  matched_cells <- rowSums(cm[, pop_cols, drop = FALSE])
  total_cells   <- cm$total_cells
  pct_assigned  <- mean(matched_cells / total_cells, na.rm = TRUE) * 100
  message(sprintf("[Step 04] Mean fraction of cells in matched populations: %.1f%%", pct_assigned))

  if (pct_assigned < 50) {
    warning("[Step 04] Less than 50% of cells assigned to matched populations. ",
            "Consider lowering integration.similarity_threshold.", call. = FALSE)
  }

  saveRDS(freq_result, file.path(out_dir, "frequency_matrix.rds"))
  write.csv(fm, file.path(out_dir, "frequency_matrix.csv"), row.names = FALSE)
  message(sprintf("[Step 04] Frequency matrix saved: %d samples x %d populations",
                  nrow(fm), length(pop_cols)))

  message("\n=== STEP 4 COMPLETE ===\n")

}, error = function(e) {
  message("\n[Error] Step 04 failed:")
  message(e$message)
  stop(paste("[Step 04 Fatal]", e$message), call. = FALSE)
})
