#!/usr/bin/env Rscript
# src/legacy/steps/02_correct.R
# ==============================================================================
# [LEGACY] STEP 02: CYCOMBINE BATCH CORRECTION — MFI ALIGNMENT PATH
# Description: Loads uncorrected SCE, performs cyCombine EMD correction, saves
#              result. Only executed when clustering.run_per_batch = false.
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
  library(cyCombine)
  library(tidyverse)
})

source("src/legacy/functions/correction.R")

message("\n=== [LEGACY] PIPELINE STEP 2: BATCH CORRECTION ===")

config   <- read_yaml("config/global_params.yml")
in_file  <- file.path(config$directories$intermediate, "uncorrected_sce.rds")
out_file <- file.path(config$directories$processed,    "corrected_sce.rds")

if (!file.exists(in_file)) stop("[Fatal] Step 1 output not found: ", in_file)
if (!dir.exists(dirname(out_file))) dir.create(dirname(out_file), recursive = TRUE)

tryCatch({

  message("[IO] Loading uncorrected data...")
  sce              <- readRDS(in_file)
  markers_to_correct <- rownames(sce)

  message(sprintf("[Config] Correcting %d markers across %d batches.",
                  length(markers_to_correct), length(unique(sce$batch))))

  sce_corrected <- run_batch_correction(
    sce       = sce,
    markers   = markers_to_correct,
    batch_col = "batch",
    config    = config$correction
  )

  var_pre  <- mean(apply(assay(sce,           "exprs"), 1L, var))
  var_post <- mean(apply(assay(sce_corrected, "exprs"), 1L, var))
  message(sprintf("[QC] Mean variance: %.3f -> %.3f", var_pre, var_post))
  if (var_post < 1e-4) warning("[QC] Post-correction variance is extremely low!", call. = FALSE)

  message(sprintf("[IO] Saving corrected data to %s", out_file))
  saveRDS(sce_corrected, out_file)

  message("=== [LEGACY] STEP 2 COMPLETE ===\n")

}, error = function(e) {
  message("\n[Error] Step 2 failed:")
  message(e$message)
  stop(paste("[Step 02 Fatal]", e$message), call. = FALSE)
})
