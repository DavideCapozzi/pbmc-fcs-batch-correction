#!/usr/bin/env Rscript

# src/steps/01_load.R
# ==============================================================================
# STEP 01: DATA INGESTION & TRANSFORMATION
# Loads raw FCS files, detects shared biological markers, applies arcsinh
# transformation, and builds the intermediate SCE. If Step 02 produced
# qc_cell_filters.rds, those per-sample cell index filters are applied before
# the SCE is constructed, removing doublets/debris identified by QC gating.
# Output: results/intermediate/filtered_sce.rds
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
  library(flowCore)
  library(SingleCellExperiment)
})

source("src/functions/fcs_io.R")
source("src/functions/step_logger.R")

message("\n=== PIPELINE STEP 1: DATA INGESTION ===")

config_path <- "config/global_params.yml"
if (!file.exists(config_path)) stop("[Fatal] Config file not found: ", config_path)
config <- read_yaml(config_path)

raw_dir  <- config$directories$raw
out_dir  <- config$directories$intermediate
proc_dir <- config$directories$processed

for (folder_name in names(raw_dir)) {
  if (!dir.exists(raw_dir[[folder_name]])) stop("[Fatal] Raw data directory not found for: ", folder_name)
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cofactor_val <- config$markers$transform_cofactor
message(sprintf("[Config] Transformation cofactor: %s", cofactor_val))

if (!is.na(cofactor_val) && cofactor_val < 10) {
  message("[WARNING] Cofactor < 10 is typical for CyTOF. For FACS use 150.")
}

log_obj <- init_step_log(
  step_name   = "01_load",
  step_number = 1L,
  input_files = unlist(lapply(raw_dir, function(d) list.files(d, "\\.fcs$", full.names = TRUE)))
)

tryCatch({

  is_test  <- isTRUE(config$testing$enabled)
  test_lim <- as.integer(config$testing$max_files_per_panel %||% 5L)

  cell_tbl <- read_and_prep_data(
    data_dirs  = raw_dir,
    cofactor   = cofactor_val,
    exclude    = config$markers$exclude_channels,
    test_mode  = is_test,
    test_limit = test_lim
  )

  n_cells_raw <- nrow(cell_tbl)
  add_metric(log_obj, "n_cells_before_qc_filter", n_cells_raw)

  # Apply QC cell filters produced by Step 02 (if present)
  qc_filter_file <- file.path(proc_dir, "qc_cell_filters.rds")
  if (file.exists(qc_filter_file)) {
    message("[Step 01] QC cell filters detected — applying per-sample filters.")
    qc_filters   <- readRDS(qc_filter_file)
    filtered_rows <- integer(0L)

    for (sid in unique(cell_tbl$sample_id)) {
      f <- qc_filters[[sid]]
      if (is.null(f)) {
        message(sprintf("   [QC] Sample '%s': EXCLUDED (QC step removed it).", sid))
        next
      }
      sample_rows     <- which(cell_tbl$sample_id == sid)
      valid_in_sample <- f$valid_indices[f$valid_indices <= length(sample_rows)]
      filtered_rows   <- c(filtered_rows, sample_rows[valid_in_sample])
    }

    cell_tbl  <- cell_tbl[filtered_rows, ]
    n_after   <- nrow(cell_tbl)
    pct_removed <- 100 * (1 - n_after / n_cells_raw)
    message(sprintf("   [QC] Cells: %d → %d  (%.1f%% removed by QC filters)",
                    n_cells_raw, n_after, pct_removed))
    add_metric(log_obj, "n_cells_after_qc_filter", n_after)
    add_metric(log_obj, "pct_qc_removed",           round(pct_removed, 2))
  } else {
    message("[Step 01] No QC filters found — loading all events (run Step 02 first for QC).")
    add_metric(log_obj, "n_cells_after_qc_filter", n_cells_raw)
  }

  # Integrity checks
  batch_counts <- table(cell_tbl$batch)
  message("[Step 01] Cell counts per batch:")
  print(batch_counts)

  for (b in unique(cell_tbl$batch)) {
    samples_in_batch <- unique(cell_tbl$sample_id[cell_tbl$batch == b])
    message(sprintf("   [Batch: %s] %s", b, paste(samples_in_batch, collapse = ", ")))
  }

  n_markers <- ncol(cell_tbl) - 2L
  message(sprintf("[Step 01] Dimensions: %d cells x %d biological markers", nrow(cell_tbl), n_markers))

  if (n_markers < 3L) {
    stop("[Fatal] Too few shared markers detected. Check 'exclude_channels' in config.")
  }

  add_metric(log_obj, "n_markers",        n_markers)
  add_metric(log_obj, "n_cells_final",    nrow(cell_tbl))
  add_metric(log_obj, "samples_per_batch",
             lapply(split(cell_tbl$sample_id, cell_tbl$batch),
                    function(x) length(unique(x))))

  # Build SCE (assay: Features x Cells)
  meta_cols  <- c("sample_id", "batch")
  exprs_data <- as.matrix(cell_tbl[, setdiff(colnames(cell_tbl), meta_cols)])
  meta_data  <- cell_tbl[, meta_cols]

  sce <- SingleCellExperiment(
    assays  = list(exprs = t(exprs_data)),
    colData = meta_data
  )

  out_file <- file.path(out_dir, "filtered_sce.rds")
  saveRDS(sce, out_file)
  message(sprintf("[Step 01] Artifact saved: %s", out_file))

  finalize_step_log(log_obj, output_files = out_file, status = "SUCCESS")
  write_step_json(log_obj, config$directories$logs)

  message("=== STEP 1 COMPLETE ===\n")

}, error = function(e) {
  add_metric(log_obj, "error_message", e$message)
  finalize_step_log(log_obj, output_files = character(0L), status = "FAILURE")
  write_step_json(log_obj, config$directories$logs)
  message("\n[Error] Step 01 failed: ", e$message)
  stop(paste("[Step 01 Fatal]", e$message), call. = FALSE)
})
