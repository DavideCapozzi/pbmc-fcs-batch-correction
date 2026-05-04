#!/usr/bin/env Rscript
# src/steps/08_dimred.R
# ==============================================================================
# STEP 08: DIMENSIONALITY REDUCTION (UMAP QC VISUALIZATION)
# Generates UMAP embeddings for visual QC.
#   - Per-batch mode: independent UMAPs per batch, matched-population coloring,
#     side-by-side batch comparison.
#   - Legacy mode: pre/post MFI-correction UMAPs.
# Output: PDF figures in results/figures/
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
  library(scater)
  library(ggplot2)
  library(patchwork)
})

source("src/functions/dimred.R")
source("src/functions/step_logger.R")

message("\n=== PIPELINE STEP 8: DIMENSIONALITY REDUCTION (UMAP QC) ===")

config        <- read_yaml("config/global_params.yml")
run_per_batch <- isTRUE(config$clustering$run_per_batch)
out_dir       <- config$directories$processed
fig_dir       <- config$directories$figures
umap_cfg      <- config$dimred %||% list(n_neighbors = 15, min_dist = 0.1, seed = 1234)

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

log_obj <- init_step_log(
  step_name   = "08_dimred",
  step_number = 8L,
  input_files = character(0L)
)

umap_files_written <- character(0L)

tryCatch({

  if (run_per_batch) {

    message("[Step 08] Mode: per-batch")

    batches   <- names(config$directories$raw)
    sce_files <- file.path(out_dir, paste0("sce_batch_", batches, ".rds"))
    missing   <- sce_files[!file.exists(sce_files)]
    if (length(missing) > 0) {
      stop("[Step 08] Missing per-batch SCE files: ", paste(missing, collapse = ", "),
           ". Run Step 03 first.")
    }

    add_metric(log_obj, "inputs", sce_files)

    sce_list <- setNames(lapply(sce_files, readRDS), batches)

    # Pre-clustering UMAP on full filtered SCE for reference
    in_pre <- file.path(config$directories$intermediate, "filtered_sce.rds")
    if (file.exists(in_pre)) {
      sce_pre <- readRDS(in_pre)
      sce_pre <- run_umap_generation(sce_pre, rownames(sce_pre), umap_cfg)
      out_f   <- file.path(fig_dir, "umap_00_pre_qc_by_batch.pdf")
      plot_umap(sce_pre, "batch", out_f)
      umap_files_written <- c(umap_files_written, out_f)
    }

    # Per-batch UMAPs
    for (b in batches) {
      sce_b      <- sce_list[[b]]
      markers    <- rownames(sce_b)
      sce_b_umap <- run_umap_generation(sce_b, markers, umap_cfg)
      sce_list[[b]] <- sce_b_umap

      f1 <- file.path(fig_dir, paste0("umap_01_", b, "_by_batch.pdf"))
      f2 <- file.path(fig_dir, paste0("umap_02_", b, "_metaclusters.pdf"))
      plot_umap(sce_b_umap, "batch",         f1)
      plot_umap(sce_b_umap, "metacluster_id", f2)
      umap_files_written <- c(umap_files_written, f1, f2)

      if ("matched_population_id" %in% names(colData(sce_b_umap))) {
        f3 <- file.path(fig_dir, paste0("umap_03_", b, "_matched_populations.pdf"))
        plot_umap_matched_populations(sce_b_umap, b, f3)
        umap_files_written <- c(umap_files_written, f3)
      }
    }

    if (all(sapply(sce_list, function(s) "matched_population_id" %in% names(colData(s))))) {
      f4 <- file.path(fig_dir, "umap_04_batch_comparison.pdf")
      plot_umap_batch_comparison(sce_list, f4)
      umap_files_written <- c(umap_files_written, f4)
    }

    # Save per-batch SCEs with UMAP coordinates
    dimred_files <- character(length(batches))
    for (i in seq_along(batches)) {
      b       <- batches[i]
      out_b   <- file.path(out_dir, paste0("dimred_sce_batch_", b, ".rds"))
      saveRDS(sce_list[[b]], out_b)
      dimred_files[i] <- out_b
    }

  } else {

    message("[Step 08] Mode: legacy (MFI correction path)")

    in_pre  <- file.path(config$directories$intermediate, "filtered_sce.rds")
    in_post <- file.path(out_dir, "clustered_sce.rds")

    if (!file.exists(in_pre) || !file.exists(in_post)) {
      stop("[Step 08] Missing input SCE files. Run Steps 01-03 (legacy) first.")
    }

    add_metric(log_obj, "inputs", c(in_pre, in_post))

    sce_pre <- readRDS(in_pre)
    sce_pre <- run_umap_generation(sce_pre, rownames(sce_pre), umap_cfg)
    f1      <- file.path(fig_dir, "umap_01_pre_correction_by_batch.pdf")
    plot_umap(sce_pre, "batch", f1)
    umap_files_written <- c(umap_files_written, f1)

    sce_post <- readRDS(in_post)
    sce_post <- run_umap_generation(sce_post, rownames(sce_post), umap_cfg)
    f2       <- file.path(fig_dir, "umap_02_post_correction_by_batch.pdf")
    f3       <- file.path(fig_dir, "umap_03_post_correction_split_by_batch.pdf")
    plot_umap(sce_post, "batch", f2)
    plot_umap_split(sce_post, "sample_id", "batch", f3)
    umap_files_written <- c(umap_files_written, f2, f3)

    if ("metacluster_id" %in% names(colData(sce_post))) {
      f4 <- file.path(fig_dir, "umap_04_post_correction_clusters.pdf")
      plot_umap(sce_post, "metacluster_id", f4)
      umap_files_written <- c(umap_files_written, f4)
    }

    saveRDS(sce_post, file.path(out_dir, "dimred_sce.rds"))
    dimred_files <- file.path(out_dir, "dimred_sce.rds")
  }

  message(sprintf("[Step 08] UMAP figures written: %d", length(umap_files_written)))

  add_metric(log_obj, "n_umap_files",    length(umap_files_written))
  add_metric(log_obj, "umap_files",      umap_files_written)

  all_outputs <- c(umap_files_written, dimred_files)
  finalize_step_log(log_obj, output_files = all_outputs, status = "SUCCESS")
  write_step_json(log_obj, config$directories$logs)

  message("=== STEP 8 COMPLETE ===\n")

}, error = function(e) {
  add_metric(log_obj, "error_message", e$message)
  finalize_step_log(log_obj, output_files = character(0L), status = "FAILURE")
  write_step_json(log_obj, config$directories$logs)
  message("\n[Error] Step 08 failed: ", e$message)
  stop(paste("[Step 08 Fatal]", e$message), call. = FALSE)
})
