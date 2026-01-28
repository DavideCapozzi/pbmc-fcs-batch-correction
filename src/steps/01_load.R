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
  library(cyCombine) # Validates installation immediately
})

# Load local functions
source("src/functions/fcs_io.R")

message("\n=== PIPELINE STEP 1: DATA INGESTION ===")

# 1. Load Configuration
# ------------------------------------------------------------------------------
config_path <- "config/params.yaml"
if (!file.exists(config_path)) stop("[Fatal] Config file not found: ", config_path)
config <- read_yaml(config_path)

raw_dir  <- config$directories$raw
out_dir  <- config$directories$intermediate

# Verify Input Directory
if (!dir.exists(raw_dir)) stop("[Fatal] Raw data directory not found: ", raw_dir)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 2. Execution
# ------------------------------------------------------------------------------
message(sprintf("[Config] Transformation cofactor: %s", config$markers$transform_cofactor))

tryCatch({
  
  # Orchestrate Data Loading (Dynamic Marker Detection internal)
  uncorrected_tbl <- read_and_prep_data(
    data_dir = raw_dir,
    patterns = config$batches,
    cofactor = config$markers$transform_cofactor,
    exclude  = config$markers$exclude_channels
  )
  
  # 3. Quality Check (Basic)
  # ----------------------------------------------------------------------------
  message("[QC] Performing basic integrity checks...")
  
  # Check Batch Distribution
  batch_counts <- table(uncorrected_tbl$batch)
  message("    -> Cell counts per batch:")
  print(batch_counts)
  
  # Check Dimensions
  n_cells <- nrow(uncorrected_tbl)
  n_markers <- ncol(uncorrected_tbl) - 2 # excluding sample_id and batch
  message(sprintf("    -> Total Dimensions: %d cells x %d biological markers", n_cells, n_markers))

  if (any(names(batch_counts) == "Unknown")) {
    stop("[Fatal] 'Unknown' batches detected. Please check filename patterns in config.")
  }

  if (n_markers < 3) {
    stop("[Fatal] Too few markers detected. Check 'exclude_channels' in config.")
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
  quit(status = 1)
})
