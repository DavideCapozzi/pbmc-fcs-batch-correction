#!/usr/bin/env Rscript

# src/steps/02b_evaluate.R
# ==============================================================================
# STEP 02b: BATCH CORRECTION EVALUATION
# Description: Generates quantitative reports on batch effect removal vs 
#              biological preservation.
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
  library(writexl)
})

source("src/functions/evaluation.R")

message("\n=== PIPELINE STEP 2b: CORRECTION EVALUATION ===")

# 1. Load Config
config <- read_yaml("config/global_params.yml")
in_pre   <- file.path(config$directories$intermediate, "uncorrected_sce.rds")
in_post  <- file.path(config$directories$processed, "corrected_sce.rds")

# We save the report in the processed directory (or figures, depending on preference)
out_report <- file.path(config$directories$processed, "batch_correction_metrics.xlsx")

# 2. Checks
if (!file.exists(in_pre)) stop("[Fatal] Step 1 output not found: ", in_pre)
if (!file.exists(in_post)) stop("[Fatal] Step 2 output not found: ", in_post)

tryCatch({
  # 3. Load Data
  message("[IO] Loading PRE and POST correction data...")
  sce_pre  <- readRDS(in_pre)
  sce_post <- readRDS(in_post)
  
  # Ensure markers are aligned
  markers <- rownames(sce_post)
  
  # 4. Run Evaluation Module
  metrics_list <- evaluate_batch_correction(
    sce_pre = sce_pre, 
    sce_post = sce_post, 
    markers = markers
  )
  
  # 5. Export Report
  message(sprintf("[IO] Writing quantitative report to: %s", out_report))
  writexl::write_xlsx(metrics_list, path = out_report)
  
  message("=== STEP 2b COMPLETE ===\n")
  
}, error = function(e) {
  message("\n[Error] Step 2b Failed:")
  message(e$message)
  stop(paste("[Step 2b Fatal]", e$message), call. = FALSE)
})