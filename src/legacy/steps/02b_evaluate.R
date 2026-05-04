#!/usr/bin/env Rscript
# src/legacy/steps/02b_evaluate.R
# ==============================================================================
# [LEGACY] STEP 02b: BATCH CORRECTION EVALUATION — MFI ALIGNMENT PATH
# Description: Quantitative report on batch effect removal vs biological
#              preservation (ANOVA F-stats + Spearman ρ). Only executed when
#              clustering.run_per_batch = false.
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
  library(writexl)
})

source("src/legacy/functions/evaluation.R")

message("\n=== [LEGACY] PIPELINE STEP 2b: CORRECTION EVALUATION ===")

config     <- read_yaml("config/global_params.yml")
in_pre     <- file.path(config$directories$intermediate, "uncorrected_sce.rds")
in_post    <- file.path(config$directories$processed,    "corrected_sce.rds")
out_report <- file.path(config$directories$processed,    "batch_correction_metrics.xlsx")

if (!file.exists(in_pre))  stop("[Fatal] Step 1 output not found: ", in_pre)
if (!file.exists(in_post)) stop("[Fatal] Step 2 output not found: ", in_post)

tryCatch({

  message("[IO] Loading PRE and POST correction data...")
  sce_pre  <- readRDS(in_pre)
  sce_post <- readRDS(in_post)
  markers  <- rownames(sce_post)

  metrics_list <- evaluate_batch_correction(sce_pre = sce_pre, sce_post = sce_post,
                                            markers  = markers)

  message(sprintf("[IO] Writing quantitative report to: %s", out_report))
  writexl::write_xlsx(metrics_list, path = out_report)

  message("=== [LEGACY] STEP 2b COMPLETE ===\n")

}, error = function(e) {
  message("\n[Error] Step 2b failed:")
  message(e$message)
  stop(paste("[Step 02b Fatal]", e$message), call. = FALSE)
})
