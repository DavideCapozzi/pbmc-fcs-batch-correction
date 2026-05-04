#!/usr/bin/env Rscript
# src/steps/04_dimred.R
# ==============================================================================
# STEP 04 (DIMRED): DIMENSIONALITY REDUCTION (UMAP QC)
# Description: Generates UMAP embeddings for QC visualization.
#   - Per-batch mode (run_per_batch=TRUE): independent UMAPs per batch,
#     matched-population coloring, side-by-side comparison.
#   - Legacy mode (run_per_batch=FALSE): pre/post correction UMAPs as before.
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
  library(scater)
  library(ggplot2)
  library(patchwork)
})

source("src/functions/dimred.R")

message("\n=== PIPELINE STEP DIMRED: DIMENSIONALITY REDUCTION (QC) ===")

config        <- read_yaml("config/global_params.yml")
run_per_batch <- isTRUE(config$clustering$run_per_batch)
out_dir       <- config$directories$processed
fig_dir       <- config$directories$figures
umap_cfg      <- config$dimred

if (is.null(umap_cfg)) umap_cfg <- list(n_neighbors = 15, min_dist = 0.1, seed = 1234)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

tryCatch({

  if (run_per_batch) {

    # --------------------------------------------------------------------------
    # PER-BATCH MODE: independent UMAP per batch + matched population coloring
    # --------------------------------------------------------------------------
    message("[DimRed] Mode: per-batch (new integration paradigm)")

    batches   <- names(config$directories$raw)
    sce_files <- file.path(out_dir, paste0("sce_batch_", batches, ".rds"))
    missing   <- sce_files[!file.exists(sce_files)]
    if (length(missing) > 0) {
      stop("[DimRed] Missing per-batch SCE files: ", paste(missing, collapse = ", "),
           ". Run Step 03 first.")
    }

    sce_list <- setNames(lapply(sce_files, readRDS), batches)

    # Also generate pre-correction UMAP on the full uncorrected SCE for reference
    in_pre <- file.path(config$directories$intermediate, "uncorrected_sce.rds")
    if (file.exists(in_pre)) {
      sce_pre <- readRDS(in_pre)
      sce_pre <- run_umap_generation(sce_pre, rownames(sce_pre), umap_cfg)
      plot_umap(sce_pre, "batch",
                file.path(fig_dir, "umap_00_PRE_correction_by_batch.pdf"))
    }

    # Per-batch UMAPs
    for (b in batches) {
      sce_b   <- sce_list[[b]]
      markers <- rownames(sce_b)

      sce_b_umap <- run_umap_generation(sce_b, markers, umap_cfg)
      sce_list[[b]] <- sce_b_umap

      plot_umap(sce_b_umap, "batch",
                file.path(fig_dir, paste0("umap_01_", b, "_by_batch.pdf")))

      plot_umap(sce_b_umap, "metacluster_id",
                file.path(fig_dir, paste0("umap_02_", b, "_metaclusters.pdf")))

      if ("matched_population_id" %in% names(colData(sce_b_umap))) {
        plot_umap_matched_populations(
          sce_b_umap, b,
          file.path(fig_dir, paste0("umap_03_", b, "_matched_populations.pdf"))
        )
      }
    }

    # Side-by-side comparison with consistent color palette
    if (all(sapply(sce_list, function(s) "matched_population_id" %in% names(colData(s))))) {
      plot_umap_batch_comparison(
        sce_list,
        file.path(fig_dir, "umap_04_batch_comparison.pdf")
      )
    }

    # Save per-batch SCEs with UMAP coordinates
    for (b in batches) {
      saveRDS(sce_list[[b]],
              file.path(out_dir, paste0("dimred_sce_batch_", b, ".rds")))
    }

  } else {

    # --------------------------------------------------------------------------
    # LEGACY MODE: pre/post correction UMAPs
    # --------------------------------------------------------------------------
    message("[DimRed] Mode: legacy (MFI correction path)")

    in_pre  <- file.path(config$directories$intermediate, "uncorrected_sce.rds")
    in_post <- file.path(out_dir, "clustered_sce.rds")

    if (!file.exists(in_pre) || !file.exists(in_post)) {
      stop("[DimRed] Missing input SCE files. Run steps 1-3 (legacy) first.")
    }

    sce_pre <- readRDS(in_pre)
    sce_pre <- run_umap_generation(sce_pre, rownames(sce_pre), umap_cfg)
    plot_umap(sce_pre, "batch", file.path(fig_dir, "umap_01_PRE_correction_by_batch.pdf"))

    sce_post <- readRDS(in_post)
    sce_post <- run_umap_generation(sce_post, rownames(sce_post), umap_cfg)
    plot_umap(sce_post, "batch", file.path(fig_dir, "umap_02_POST_correction_by_batch.pdf"))
    plot_umap_split(sce_post, "sample_id", "batch",
                    file.path(fig_dir, "umap_03_POST_correction_split_by_batch.pdf"))

    if ("metacluster_id" %in% names(colData(sce_post))) {
      plot_umap(sce_post, "metacluster_id",
                file.path(fig_dir, "umap_04_POST_correction_clusters.pdf"))
    }

    saveRDS(sce_post, file.path(out_dir, "dimred_sce.rds"))
  }

  message("=== DIMRED STEP COMPLETE ===\n")

}, error = function(e) {
  message("\n[Error] DimRed step failed:")
  message(e$message)
  stop(paste("[DimRed Fatal]", e$message), call. = FALSE)
})
