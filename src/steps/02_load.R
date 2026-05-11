#!/usr/bin/env Rscript

# src/steps/02_load.R
# ==============================================================================
# STEP 02: DATA INGESTION & TRANSFORMATION
# Loads raw FCS files, detects shared biological markers, applies arcsinh
# transformation, and builds the intermediate SCE. If Step 01 produced
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

source("src/functions/utils.R")
source("src/functions/parallel_utils.R")
source("src/functions/fcs_io.R")
source("src/functions/step_logger.R")

message("\n=== PIPELINE STEP 2: DATA INGESTION ===")

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
  step_name   = "02_load",
  step_number = 2L,
  input_files = unlist(lapply(raw_dir, function(d) list.files(d, "\\.fcs$", full.names = TRUE)))
)

tryCatch({

  is_test  <- isTRUE(config$testing$enabled)
  test_lim <- as.integer(config$testing$max_files_per_panel %||% 5L)
  n_workers <- get_n_workers(config)

  cell_tbl <- read_and_prep_data(
    data_dirs    = raw_dir,
    cofactor     = cofactor_val,
    exclude      = config$markers$exclude_channels,
    test_mode    = is_test,
    test_limit   = test_lim,
    raw_filters  = config$directories$raw_filters  %||% list(),
    batch_labels = config$batch_labels             %||% list(),
    n_workers    = n_workers
  )

  n_cells_raw <- nrow(cell_tbl)
  add_metric(log_obj, "n_cells_before_qc_filter", n_cells_raw)

  # Apply QC cell filters produced by Step 01 (if present)
  qc_filter_file <- file.path(proc_dir, "qc_cell_filters.rds")
  if (file.exists(qc_filter_file)) {
    message("[Step 02] QC cell filters detected — applying per-sample filters.")
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
    warning(
      "[Step 02] qc_cell_filters.rds not found — loading ALL events without QC filtering.\n",
      "          Run Step 01 (01_qc.R) first to apply per-sample cell filters.",
      call. = FALSE
    )
    add_metric(log_obj, "n_cells_after_qc_filter", n_cells_raw)
  }

  # Integrity checks
  batch_counts <- table(cell_tbl$batch)
  message("[Step 02] Cell counts per batch:")
  print(batch_counts)

  for (b in unique(cell_tbl$batch)) {
    samples_in_batch <- unique(cell_tbl$sample_id[cell_tbl$batch == b])
    message(sprintf("   [Batch: %s] %s", b, paste(samples_in_batch, collapse = ", ")))
  }

  n_markers <- ncol(cell_tbl) - 2L
  message(sprintf("[Step 02] Dimensions: %d cells x %d biological markers", nrow(cell_tbl), n_markers))

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

  # CD3+ gate (optional): restrict SCE to CD3+ T cells before clustering.
  # Frequencies after this gate are expressed as % of CD3+ cells, not total PBMCs —
  # consistent with standard clinical immunology reporting (EuroFlow, NCI BioBank).
  if (isTRUE(config$cd3_gate$enabled) && "CD3" %in% rownames(sce)) {
    message("[Step 02] CD3+ gate enabled — restricting to CD3+ T cells.")

    cd3_gap_thr  <- as.numeric(config$cd3_gate$uninformative_gap %||% 1.0)
    cd3_fallback <- as.numeric(config$cd3_gate$fallback_threshold %||% 4.0)
    min_cd3_pct  <- as.numeric(config$cd3_gate$min_cd3_pct        %||% 10.0)

    cd3_expr   <- assay(sce, "exprs")["CD3", ]
    sample_ids <- sce$sample_id
    keep_cells <- logical(ncol(sce))
    cd3_metrics <- list()

    for (sid in unique(sample_ids)) {
      idx      <- which(sample_ids == sid)
      vals     <- cd3_expr[idx]
      n_before <- length(vals)

      thr <- tryCatch({
        km   <- kmeans(vals, centers = 2L, nstart = 10L, iter.max = 50L)
        ctrs <- sort(km$centers[, 1L])
        gap  <- ctrs[2L] - ctrs[1L]
        if (gap < cd3_gap_thr) {
          message(sprintf("   [CD3 gate] '%s': k-means gap=%.2f < %.2f — fallback thr=%.2f",
                          sid, gap, cd3_gap_thr, cd3_fallback))
          cd3_fallback
        } else {
          mean(ctrs)
        }
      }, error = function(e) {
        warning(sprintf("[CD3 gate] k-means failed for '%s': %s — fallback thr=%.2f",
                        sid, e$message, cd3_fallback), call. = FALSE)
        cd3_fallback
      })

      pos_idx <- idx[vals > thr]
      n_after <- length(pos_idx)
      pct_cd3 <- 100 * n_after / n_before

      if (pct_cd3 < min_cd3_pct) {
        warning(sprintf("[CD3 gate] '%s': only %.1f%% CD3+ (thr=%.2f) — verify acquisition.",
                        sid, pct_cd3, thr), call. = FALSE)
      }
      min_cells_ok <- as.integer(config$qc$min_cells_final %||% 5000L)
      if (n_after < min_cells_ok) {
        warning(sprintf("[CD3 gate] '%s': only %d CD3+ cells remain (min=%d) — too few for clustering.",
                        sid, n_after, min_cells_ok), call. = FALSE)
      }

      keep_cells[pos_idx] <- TRUE
      cd3_metrics[[sid]]  <- list(
        n_before  = n_before,
        n_after   = n_after,
        pct_cd3   = round(pct_cd3, 1L),
        threshold = round(thr, 3L)
      )
      message(sprintf("   [CD3 gate] %-35s %d → %d (%.1f%% CD3+, thr=%.2f)",
                      sid, n_before, n_after, pct_cd3, thr))
    }

    n_before_gate <- ncol(sce)
    sce           <- sce[, keep_cells]
    n_after_gate  <- ncol(sce)

    message(sprintf("[Step 02] CD3 gate: %d → %d cells (%.1f%% retained)",
                    n_before_gate, n_after_gate,
                    100 * n_after_gate / n_before_gate))

    add_metric(log_obj, "cd3_gate_applied",      TRUE)
    add_metric(log_obj, "cd3_gate_n_before",     n_before_gate)
    add_metric(log_obj, "cd3_gate_n_after",      n_after_gate)
    add_metric(log_obj, "cd3_gate_pct_retained", round(100 * n_after_gate / n_before_gate, 1L))
    add_metric(log_obj, "cd3_gate_per_sample",   cd3_metrics)
    add_metric(log_obj, "n_cells_final",         n_after_gate)
  } else {
    add_metric(log_obj, "cd3_gate_applied", FALSE)
  }

  # CD8+ gate (optional): further restrict to CD8+ T cells after the CD3 gate.
  # Aligns the denominator with DuraClone TEX panel design (CD8 exhaustion panel).
  # Frequencies after this gate are expressed as % of CD8+ T cells.
  if (isTRUE(config$cd8_gate$enabled) && "CD8" %in% rownames(sce)) {
    message("[Step 02] CD8+ gate enabled — restricting to CD8+ T cells.")

    cd8_gap_thr  <- as.numeric(config$cd8_gate$uninformative_gap %||% 1.0)
    cd8_fallback <- as.numeric(config$cd8_gate$fallback_threshold %||% 3.5)
    min_cd8_pct  <- as.numeric(config$cd8_gate$min_cd8_pct        %||% 10.0)

    cd8_expr   <- assay(sce, "exprs")["CD8", ]
    sample_ids <- sce$sample_id
    keep_cells <- logical(ncol(sce))
    cd8_metrics <- list()

    for (sid in unique(sample_ids)) {
      idx      <- which(sample_ids == sid)
      vals     <- cd8_expr[idx]
      n_before <- length(vals)

      thr <- tryCatch({
        km   <- kmeans(vals, centers = 2L, nstart = 10L, iter.max = 50L)
        ctrs <- sort(km$centers[, 1L])
        gap  <- ctrs[2L] - ctrs[1L]
        if (gap < cd8_gap_thr) {
          message(sprintf("   [CD8 gate] '%s': k-means gap=%.2f < %.2f — fallback thr=%.2f",
                          sid, gap, cd8_gap_thr, cd8_fallback))
          cd8_fallback
        } else {
          mean(ctrs)
        }
      }, error = function(e) {
        warning(sprintf("[CD8 gate] k-means failed for '%s': %s — fallback thr=%.2f",
                        sid, e$message, cd8_fallback), call. = FALSE)
        cd8_fallback
      })

      pos_idx <- idx[vals > thr]
      n_after <- length(pos_idx)
      pct_cd8 <- 100 * n_after / n_before

      if (pct_cd8 < min_cd8_pct) {
        warning(sprintf("[CD8 gate] '%s': only %.1f%% CD8+ (thr=%.2f) — verify acquisition.",
                        sid, pct_cd8, thr), call. = FALSE)
      }
      min_cells_ok <- as.integer(config$qc$min_cells_final %||% 5000L)
      if (n_after < min_cells_ok) {
        warning(sprintf("[CD8 gate] '%s': only %d CD8+ cells remain (min=%d) — too few for clustering.",
                        sid, n_after, min_cells_ok), call. = FALSE)
      }

      keep_cells[pos_idx] <- TRUE
      cd8_metrics[[sid]]  <- list(
        n_before  = n_before,
        n_after   = n_after,
        pct_cd8   = round(pct_cd8, 1L),
        threshold = round(thr, 3L)
      )
      message(sprintf("   [CD8 gate] %-35s %d → %d (%.1f%% CD8+, thr=%.2f)",
                      sid, n_before, n_after, pct_cd8, thr))
    }

    n_before_gate <- ncol(sce)
    sce           <- sce[, keep_cells]
    n_after_gate  <- ncol(sce)

    message(sprintf("[Step 02] CD8 gate: %d → %d cells (%.1f%% retained)",
                    n_before_gate, n_after_gate,
                    100 * n_after_gate / n_before_gate))

    add_metric(log_obj, "cd8_gate_applied",      TRUE)
    add_metric(log_obj, "cd8_gate_n_before",     n_before_gate)
    add_metric(log_obj, "cd8_gate_n_after",      n_after_gate)
    add_metric(log_obj, "cd8_gate_pct_retained", round(100 * n_after_gate / n_before_gate, 1L))
    add_metric(log_obj, "cd8_gate_per_sample",   cd8_metrics)
    add_metric(log_obj, "n_cells_final",         n_after_gate)
  } else {
    add_metric(log_obj, "cd8_gate_applied", FALSE)
  }

  out_file <- file.path(out_dir, "filtered_sce.rds")
  saveRDS(sce, out_file)
  message(sprintf("[Step 02] Artifact saved: %s", out_file))

  finalize_step_log(log_obj, output_files = out_file, status = "SUCCESS")
  write_step_json(log_obj, config$directories$logs)

  message("=== STEP 2 COMPLETE ===\n")

}, error = function(e) {
  add_metric(log_obj, "error_message", e$message)
  finalize_step_log(log_obj, output_files = character(0L), status = "FAILURE")
  write_step_json(log_obj, config$directories$logs)
  message("\n[Error] Step 02 failed: ", e$message)
  stop(paste("[Step 02 Fatal]", e$message), call. = FALSE)
})
