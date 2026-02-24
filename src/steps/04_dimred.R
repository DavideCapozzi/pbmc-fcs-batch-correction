#!/usr/bin/env Rscript
# src/steps/04_dimred.R
# ==============================================================================
# STEP 04: DIMENSIONALITY REDUCTION (UMAP QC)
# Description: Generates UMAP embeddings PRE and POST batch correction.
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
  library(scater)
  library(ggplot2)
})

source("src/functions/dimred.R")

message("\n=== PIPELINE STEP 4: DIMENSIONALITY REDUCTION (QC) ===")

# 1. Load Config
config <- read_yaml("config/global_params.yml")

in_pre  <- file.path(config$directories$intermediate, "uncorrected_sce.rds")
in_post <- file.path(config$directories$processed, "clustered_sce.rds")
out_dir <- config$directories$figures 

# 2. Safety Checks
if (is.null(out_dir)) stop("[Fatal] 'directories$figures' is not defined in global_params.yml")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!file.exists(in_pre) || !file.exists(in_post)) stop("[Fatal] Missing input SCE files. Run steps 1-3 first.")

tryCatch({
  
  # Extract UMAP config
  umap_cfg <- config$dimred
  if(is.null(umap_cfg)) {
    message("[Warning] 'dimred' config missing. Using defaults.")
    umap_cfg <- list(n_neighbors=15, min_dist=0.2, seed=1234, max_cells=100000)
  }
  
  # ----------------------------------------------------------------------------
  # PART A: PRE-CORRECTION UMAP
  # ----------------------------------------------------------------------------
  message("\n[Phase A] Processing PRE-correction data...")
  sce_pre <- readRDS(in_pre)
  
  # Run UMAP (using common markers available in the pre-object)
  sce_pre <- run_umap_generation(sce_pre, rownames(sce_pre), umap_cfg)
  
  # Save Pre-correction plots (PDF)
  plot_umap(sce_pre, "batch", file.path(out_dir, "umap_01_PRE_correction_by_batch.pdf"))
  
  # ----------------------------------------------------------------------------
  # PART B: POST-CORRECTION UMAP
  # ----------------------------------------------------------------------------
  message("\n[Phase B] Processing POST-correction data...")
  sce_post <- readRDS(in_post)
  
  # Run UMAP
  sce_post <- run_umap_generation(sce_post, rownames(sce_post), umap_cfg)
  
  # Save Post-correction plots (PDF)
  # 1. Basic batch overlap
  plot_umap(sce_post, "batch", file.path(out_dir, "umap_02_POST_correction_by_batch.pdf"))
  
  # 2. Split view (Faceted by batch, colored by sample_id to see condition variance)
  plot_umap_split(
    sce = sce_post, 
    color_by = "sample_id", 
    split_by = "batch", 
    out_path = file.path(out_dir, "umap_03_POST_correction_split_by_batch.pdf")
  )
  
  # 3. Biological Clusters
  if ("metacluster_id" %in% names(colData(sce_post))) {
    plot_umap(sce_post, "metacluster_id", file.path(out_dir, "umap_04_POST_correction_clusters.pdf"))
  }
  
  # ----------------------------------------------------------------------------
  # PART C: SAVE ARTIFACT
  # ----------------------------------------------------------------------------
  out_file_sce <- file.path(config$directories$processed, "dimred_sce.rds")
  message(sprintf("\n[IO] Saving post-correction UMAP coordinates to %s", out_file_sce))
  saveRDS(sce_post, out_file_sce)
  
  message("=== STEP 4 COMPLETE ===\n")
  
}, error = function(e) {
  message("\n[Error] Step 4 Failed:")
  message(e$message)
  stop(paste("[Step 4 Fatal]", e$message), call. = FALSE)
})