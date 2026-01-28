#!/usr/bin/env Rscript

# src/steps/02_correct.R
# ==============================================================================
# STEP 02: BATCH CORRECTION
# Description: Loads uncorrected SCE, performs cyCombine correction, saves result.
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
  library(cyCombine)
  library(tidyverse)
})

source("src/functions/correction.R")

message("\n=== PIPELINE STEP 2: BATCH CORRECTION ===")

# 1. Load Config
config <- read_yaml("config/global_params.yml")
in_file <- file.path(config$directories$intermediate, "uncorrected_sce.rds")
out_file <- file.path(config$directories$processed, "corrected_sce.rds")

# 2. Checks
if (!file.exists(in_file)) stop("[Fatal] Step 1 output not found: ", in_file)
if (!dir.exists(dirname(out_file))) dir.create(dirname(out_file), recursive = TRUE)

tryCatch({
  # 3. Load Data
  message("[IO] Loading uncorrected data...")
  sce <- readRDS(in_file)
  
  # Detect markers (using the columns present in the SCE assays)
  # Note: Step 1 already filtered for shared markers.
  markers_to_correct <- rownames(sce)
  
  message(sprintf("[Config] Correcting %d markers across %d batches.", 
                  length(markers_to_correct), 
                  length(unique(sce$batch))))
  
  # 4. Run Correction Module
  sce_corrected <- run_batch_correction(
    sce = sce,
    markers = markers_to_correct,
    batch_col = "batch",
    config = config$correction
  )
  
  # 5. Quick Validation
  # Compare variance before/after to ensure we didn't flatten everything to 0
  var_pre <- mean(apply(assay(sce, "exprs"), 1, var))
  var_post <- mean(apply(assay(sce_corrected, "exprs"), 1, var))
  
  message(sprintf("[QC] Mean Variance Change: %.3f -> %.3f", var_pre, var_post))
  if(var_post < 1e-4) warning("[QC Warning] Post-correction variance is extremely low!")
  
  # 6. Save
  message(sprintf("[IO] Saving corrected data to %s", out_file))
  saveRDS(sce_corrected, out_file)
  
  message("=== STEP 2 COMPLETE ===\n")
  
}, error = function(e) {
  message("\n[Error] Step 2 Failed:")
  message(e$message)
  stop(paste("[Step 2 Fatal]", e$message), call. = FALSE)
})
