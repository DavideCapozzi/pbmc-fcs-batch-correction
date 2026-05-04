#!/usr/bin/env Rscript
# src/legacy/steps/03_cluster.R
# ==============================================================================
# [LEGACY] STEP 03: FLOWSOM CLUSTERING ON COMBINED CORRECTED DATA
# Description: Performs FlowSOM clustering on cyCombine-corrected data. Only
#              executed when clustering.run_per_batch = false.
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
  library(FlowSOM)
})

source("src/functions/clustering.R")

message("\n=== [LEGACY] PIPELINE STEP 3: CLUSTERING ===")

config   <- read_yaml("config/global_params.yml")
in_file  <- file.path(config$directories$processed, "corrected_sce.rds")
out_file <- file.path(config$directories$processed, "clustered_sce.rds")

if (!file.exists(in_file)) stop("[Fatal] Step 2 output not found. Run correction first.")

tryCatch({

  message("[IO] Loading corrected data...")
  sce                <- readRDS(in_file)
  markers_clustering <- rownames(sce)

  sce_clustered <- run_flowsom_clustering(
    sce     = sce,
    markers = markers_clustering,
    config  = config$clustering
  )

  n_clus <- length(unique(sce_clustered$metacluster_id))
  message(sprintf("[QC] Generated %d metaclusters (Target: %d)",
                  n_clus, config$clustering$n_metaclusters))
  freqs <- table(sce_clustered$metacluster_id)
  message("[QC] Cluster sizes:")
  print(freqs)
  if (min(freqs) < 10L) warning("[QC] Some clusters have very few cells (< 10).", call. = FALSE)

  message(sprintf("[IO] Saving clustered data to %s", out_file))
  saveRDS(sce_clustered, out_file)

  message("=== [LEGACY] STEP 3 COMPLETE ===\n")

}, error = function(e) {
  message("\n[Error] Step 3 failed:")
  message(e$message)
  stop(paste("[Step 03 Fatal]", e$message), call. = FALSE)
})
