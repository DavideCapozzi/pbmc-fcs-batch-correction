#!/usr/bin/env Rscript
# src/steps/03_cluster_per_batch.R
# ==============================================================================
# STEP 03: INDEPENDENT PER-BATCH FLOWSOM CLUSTERING + CLUSTER MATCHING
# Runs FlowSOM independently on each batch (no cross-batch MFI alignment).
# Matches equivalent metaclusters via cosine similarity of marker centroids.
# Annotates matched_population_id in each per-batch SCE.
# Input:  results/intermediate/filtered_sce.rds
# Output: cluster_centroids.rds, cluster_match_table.rds, sce_batch_<name>.rds
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
  library(FlowSOM)
  library(ConsensusClusterPlus)
  library(dplyr)
})

source("src/functions/clustering.R")
source("src/functions/cluster_matching.R")
source("src/functions/step_logger.R")

message("\n=== PIPELINE STEP 3: PER-BATCH CLUSTERING + CLUSTER MATCHING ===")

config  <- read_yaml("config/global_params.yml")
in_file <- file.path(config$directories$intermediate, "filtered_sce.rds")
out_dir <- config$directories$processed

if (!isTRUE(config$clustering$run_per_batch)) {
  stop("[Step 03] clustering.run_per_batch is not TRUE. Use legacy 03_cluster.R for the MFI path.",
       call. = FALSE)
}

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!file.exists(in_file)) {
  stop("[Step 03] filtered_sce.rds not found. Run Steps 01 and 02 first.")
}

log_obj <- init_step_log(
  step_name   = "03_cluster_per_batch",
  step_number = 3L,
  input_files = in_file
)

tryCatch({

  sce_full <- readRDS(in_file)

  if (ncol(sce_full) == 0L) stop("[Step 03] filtered_sce.rds contains zero cells.")
  batches <- sort(unique(as.character(sce_full$batch)))
  if (length(batches) < 2L) {
    stop(sprintf(
      "[Step 03] At least 2 batches required for cluster matching; found: %s",
      paste(batches, collapse = ", ")
    ))
  }

  markers  <- rownames(sce_full)
  message(sprintf("[Step 03] Detected %d batch(es): %s", length(batches), paste(batches, collapse = ", ")))
  message(sprintf("[Step 03] Markers: %s", paste(markers, collapse = ", ")))

  add_metric(log_obj, "n_batches", length(batches))
  add_metric(log_obj, "batches",   batches)
  add_metric(log_obj, "markers",   markers)
  add_metric(log_obj, "n_cells_total", ncol(sce_full))

  sce_batch_list   <- list()
  centroids_list   <- list()
  n_meta_per_batch <- list()

  for (b in batches) {
    message(sprintf("\n[Step 03] --- Clustering batch: %s ---", b))

    sce_b     <- sce_full[, sce_full$batch == b]
    n_cells_b <- ncol(sce_b)
    message(sprintf("[Step 03]   Cells in batch: %d", n_cells_b))

    sce_b <- run_flowsom_clustering(sce_b, markers, config$clustering)

    centroids_list[[b]] <- compute_metacluster_centroids(sce_b, markers)

    n_meta <- length(unique(sce_b$metacluster_id))
    message(sprintf("[Step 03]   Metaclusters found: %d", n_meta))
    n_meta_per_batch[[b]] <- n_meta

    sce_batch_list[[b]] <- sce_b
  }

  add_metric(log_obj, "n_metaclusters_per_batch", n_meta_per_batch)

  saveRDS(centroids_list, file.path(out_dir, "cluster_centroids.rds"))
  message("\n[Step 03] Centroids saved.")

  sim_threshold <- config$integration$similarity_threshold %||% 0.70
  message(sprintf("[Step 03] Running cosine similarity matching (threshold = %.2f)...", sim_threshold))
  match_table <- match_metaclusters(centroids_list, markers, sim_threshold)

  n_matched   <- sum(match_table$is_matched, na.rm = TRUE)
  n_pairs     <- nrow(match_table[!is.na(match_table$cluster_a) & !is.na(match_table$cluster_b), ])
  pct_matched <- if (n_pairs > 0) n_matched / n_pairs else 0

  message(sprintf("[Step 03] Matched: %d / %d pairs (%.0f%%)", n_matched, n_pairs, pct_matched * 100))

  add_metric(log_obj, "n_matched_pairs", n_matched)
  add_metric(log_obj, "n_total_pairs",   n_pairs)
  add_metric(log_obj, "pct_matched",     round(pct_matched * 100, 1))

  if (pct_matched < 0.5) {
    warning(sprintf("[Step 03] Only %.0f%% of pairs matched. Consider lowering similarity_threshold.",
                    pct_matched * 100), call. = FALSE)
  }

  validation <- validate_cluster_matches(match_table, centroids_list, markers)
  n_jsd_pass <- 0L
  if (!is.null(validation$summary) && nrow(validation$summary) > 0) {
    n_jsd_pass <- sum(validation$summary$jsd_pass, na.rm = TRUE)
    message(sprintf("[Step 03] JSD < 0.1 (well-matched): %d / %d populations",
                    n_jsd_pass, nrow(validation$summary)))
  }
  add_metric(log_obj, "n_jsd_pass", n_jsd_pass)

  saveRDS(match_table, file.path(out_dir, "cluster_match_table.rds"))
  message("[Step 03] Match table saved.")

  out_batch_files <- character(length(batches))
  names(out_batch_files) <- batches

  for (b in batches) {
    sce_b      <- sce_batch_list[[b]]
    is_batch_a <- (b == batches[1L])

    cluster_col <- if (is_batch_a) "cluster_a" else "cluster_b"
    batch_col   <- if (is_batch_a) "batch_a"   else "batch_b"

    lookup <- match_table[!is.na(match_table[[batch_col]]) & match_table[[batch_col]] == b, ]

    pop_ids <- lookup$matched_population_id[
      match(as.character(sce_b$metacluster_id), lookup[[cluster_col]])
    ]

    sce_b$matched_population_id <- as.factor(pop_ids)

    out_b <- file.path(out_dir, paste0("sce_batch_", b, ".rds"))
    saveRDS(sce_b, out_b)
    out_batch_files[[b]] <- out_b
    message(sprintf("[Step 03] Saved annotated SCE for batch '%s': %s", b, out_b))
  }

  all_outputs <- c(
    file.path(out_dir, "cluster_centroids.rds"),
    file.path(out_dir, "cluster_match_table.rds"),
    unname(out_batch_files)
  )
  finalize_step_log(log_obj, output_files = all_outputs, status = "SUCCESS")
  write_step_json(log_obj, config$directories$logs)

  message("\n=== STEP 3 COMPLETE ===\n")

}, error = function(e) {
  add_metric(log_obj, "error_message", e$message)
  finalize_step_log(log_obj, output_files = character(0L), status = "FAILURE")
  write_step_json(log_obj, config$directories$logs)
  message("\n[Error] Step 03 failed: ", e$message)
  stop(paste("[Step 03 Fatal]", e$message), call. = FALSE)
})
