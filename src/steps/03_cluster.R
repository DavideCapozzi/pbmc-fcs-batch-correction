#!/usr/bin/env Rscript

# src/steps/03_cluster.R
# ==============================================================================
# STEP 03: CLUSTERING
# Description: Performs FlowSOM clustering on corrected data.
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
  library(FlowSOM)
})

source("src/functions/clustering.R")

message("\n=== PIPELINE STEP 3: CLUSTERING ===")

# 1. Load Config
config <- read_yaml("config/global_params.yml")
in_file <- file.path(config$directories$processed, "corrected_sce.rds")
out_file <- file.path(config$directories$processed, "clustered_sce.rds")

# 2. Checks
if (!file.exists(in_file)) stop("[Fatal] Step 2 output not found. Run correction first.")

tryCatch({
  # 3. Load Data
  message("[IO] Loading corrected data...")
  sce <- readRDS(in_file)
  
  # Markers for clustering (Same as correction, usually)
  markers_clustering <- rownames(sce)
  
  # 4. Run Clustering Module
  sce_clustered <- run_flowsom_clustering(
    sce = sce,
    markers = markers_clustering,
    config = config$clustering
  )
  
  # 5. Quick Validation
  n_clus <- length(unique(sce_clustered$metacluster_id))
  message(sprintf("[QC] Generated %d metaclusters (Target: %d)", 
                  n_clus, config$clustering$n_metaclusters))
  
  # Print Frequency Table
  freqs <- table(sce_clustered$metacluster_id)
  message("[QC] Cluster sizes:")
  print(freqs)
  
  if(min(freqs) < 10) warning("[QC Warning] Some clusters have very few cells (<10).")
  
  # 6. Save
  message(sprintf("[IO] Saving clustered data to %s", out_file))
  saveRDS(sce_clustered, out_file)
  
  message("=== STEP 3 COMPLETE ===\n")
  
}, error = function(e) {
  message("\n[Error] Step 3 Failed:")
  message(e$message)
  stop(paste("[Step 3 Fatal]", e$message), call. = FALSE)
})
