#!/usr/bin/env Rscript
# src/steps/04_dimred.R
# ==============================================================================
# STEP 04: DIMENSIONALITY REDUCTION (UMAP)
# Description: Generates UMAP embeddings and plots for QC.
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
  library(scater)
  library(ggplot2) # Explicitly load ggplot2 for theme elements if needed
})

source("src/functions/dimred.R")

message("\n=== PIPELINE STEP 4: DIMENSIONALITY REDUCTION ===")

# 1. Load Config
config <- read_yaml("config/global_params.yml")

# Define Paths
in_file <- file.path(config$directories$processed, "clustered_sce.rds")

# FIXED: Use 'figures' directory as defined in the YAML structure provided.
# Previously 'results' was used, causing a NULL path error.
out_dir <- config$directories$figures 

# 2. Safety Checks
if (is.null(out_dir)) {
  stop("[Fatal] 'directories$figures' is not defined in global_params.yml")
}

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

if (!file.exists(in_file)) {
  stop("[Fatal] Input file not found: ", in_file, "\nHint: Did Step 3 (Clustering) complete successfully?")
}

tryCatch({
  
  # 3. Load Data
  message("[IO] Loading clustered data...")
  sce <- readRDS(in_file)
  
  # 4. Compute UMAP
  # Check for dimred config, set defaults if missing to prevent crashes
  umap_cfg <- config$dimred
  if(is.null(umap_cfg)) {
    message("[Warning] 'dimred' config missing. Using defaults.")
    umap_cfg <- list(n_neighbors=15, min_dist=0.2, seed=1234)
  }
  
  # Run UMAP wrapper
  sce <- run_umap_generation(sce, rownames(sce), umap_cfg)
  
  # 5. Plotting
  message(sprintf("[Output] Saving plots to: %s", out_dir))
  
  # Plot by Batch (Technical Check)
  plot_umap(sce, "batch", file.path(out_dir, "umap_batch.png"))
  
  # Plot by Metacluster (Biological Check)
  # Check if metacluster_id exists (it should from Step 3)
  if ("metacluster_id" %in% names(colData(sce))) {
    plot_umap(sce, "metacluster_id", file.path(out_dir, "umap_clusters.png"))
  } else {
    warning("[QC Warning] 'metacluster_id' not found in SCE. Skipping cluster plot.")
  }
  
  # 6. Save Updated SCE (Optional: if you want to keep the UMAP coords for later)
  # For now, we usually don't overwrite the Step 3 file unless explicitly desired.
  # We just output the figures.
  
  message("=== STEP 4 COMPLETE ===\n")
  
}, error = function(e) {
  message("\n[Error] Step 4 Failed:")
  message(e$message)
  stop(paste("[Step 4 Fatal]", e$message), call. = FALSE)
})