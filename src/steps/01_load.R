#!/usr/bin/env Rscript

# src/steps/01_load.R
# ==============================================================================
# STEP 01: DATA INGESTION & TRANSFORMATION
# Description: Loads raw FCS files, dynamically detects shared markers, 
#              applies arcsinh transformation, and creates the intermediate SCE.
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
  library(flowCore)
  library(SingleCellExperiment)
  library(cyCombine) 
})

# Load local functions
source("src/functions/fcs_io.R")

message("\n=== PIPELINE STEP 1: DATA INGESTION ===")

# 1. Load Configuration
# ------------------------------------------------------------------------------
config_path <- "config/global_params.yml"
if (!file.exists(config_path)) stop("[Fatal] Config file not found: ", config_path)
config <- read_yaml(config_path)

raw_dir  <- config$directories$raw
out_dir  <- config$directories$intermediate

# Verify Input Directories
for (folder_name in names(raw_dir)) {
  if (!dir.exists(raw_dir[[folder_name]])) stop("[Fatal] Raw data directory not found for: ", folder_name)
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 2. Execution
# ------------------------------------------------------------------------------

cofactor_val <- config$markers$transform_cofactor
  
message(sprintf("[Config] Transformation cofactor: %s", cofactor_val))

# Sanity Check for FACS methodology
if (!is.na(cofactor_val) && cofactor_val < 10) {
  message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  message("[QC WARNING] Cofactor is set to < 10 (Current: ", cofactor_val, ").")
  message("             This is typical for Mass Cytometry (CyTOF).")
  message("             For standard FACS, consider using 150 to avoid")
  message("             artificial 'peak splitting' of negative populations.")
  message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
}

tryCatch({
  
  # Check test mode defaults to avoid null errors
  is_test <- if (!is.null(config$testing$enabled)) config$testing$enabled else FALSE
  test_lim <- if (!is.null(config$testing$max_files_per_panel)) config$testing$max_files_per_panel else 5
  
  # Orchestrate Data Loading (Dynamic Marker Detection internal)
  uncorrected_tbl <- read_and_prep_data(
    data_dirs = raw_dir,
    cofactor = config$markers$transform_cofactor,
    exclude  = config$markers$exclude_channels,
    test_mode = is_test,
    test_limit = test_lim
  )
  
  # 3. Quality Check 
  # ----------------------------------------------------------------------------
  message("[QC] Performing basic integrity checks...")
  
  # Check Batch Distribution
  batch_counts <- table(uncorrected_tbl$batch)
  message("    -> Cell counts per batch (folder):")
  print(batch_counts)
  
  # Print Specific Samples Processed per Batch
  message("    -> Samples successfully mapped per batch:")
  for (b in unique(uncorrected_tbl$batch)) {
    samples_in_batch <- unique(uncorrected_tbl$sample_id[uncorrected_tbl$batch == b])
    message(sprintf("       [Batch: %s] -> %s", b, paste(samples_in_batch, collapse = ", ")))
  }
  
  # Check Dimensions
  n_cells <- nrow(uncorrected_tbl)
  n_markers <- ncol(uncorrected_tbl) - 2 # excluding sample_id and batch
  message(sprintf("    -> Total Dimensions: %d cells x %d biological markers", n_cells, n_markers))
  
  if (n_markers < 3) {
    stop("[Fatal] Too few shared markers detected between folders. Check 'exclude_channels' in config.")
  }

  # 4. Save Artifacts
  # ----------------------------------------------------------------------------
  message("[Export] Converting to SingleCellExperiment...")
  
  # Split metadata and expression for SCE construction
  meta_cols <- c("sample_id", "batch")
  
  # Extract Expression Matrix (Cells x Markers)
  exprs_data <- uncorrected_tbl %>% 
    select(-all_of(meta_cols)) %>% 
    as.matrix()
  
  # Extract Metadata
  meta_data <- uncorrected_tbl %>% 
    select(all_of(meta_cols))
  
  # Create SCE (Assays expect Features x Cells, so we transpose)
  sce <- SingleCellExperiment(
    assays = list(exprs = t(exprs_data)),
    colData = meta_data
  )
  
  out_file <- file.path(out_dir, "uncorrected_sce.rds")
  saveRDS(sce, out_file)
  
  message(sprintf("   -> Artifact saved to: %s", out_file))
  message("=== STEP 1 COMPLETE ===\n")
  
}, error = function(e) {
  message("\n[Error] Pipeline Step 1 Failed:")
  message(e$message)
  stop(paste("[Step 1 Error]", e$message), call. = FALSE)
})
