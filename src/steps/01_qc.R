#!/usr/bin/env Rscript
# src/steps/01_qc.R
# ==============================================================================
# STEP 01: PER-SAMPLE QC FILTERING
# Reads raw FCS files with scatter channels, applies PeacoQC for temporal drift
# and margin-event removal, singlet gating (FSC-H/FSC-A slope), debris removal
# (FSC-A vs SSC-A), and viability dye exclusion (ViaKrome: median + 3×MAD).
# Produces a per-sample valid cell index list that Step 02 applies when building
# the SCE.
#
# Output artifacts:
#   results/processed/qc_cell_filters.rds   — named list, one entry per sample
#   results/logs/steps/01_qc_{ts}.json      — step audit JSON
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(flowCore)
  library(PeacoQC)
})

source("src/functions/qc.R")
source("src/functions/step_logger.R")

message("\n=== PIPELINE STEP 1: QC FILTERING ===")

config <- read_yaml("config/global_params.yml")

if (!isTRUE(config$qc$enabled)) {
  message("[Step 01] qc.enabled = FALSE in config — skipping QC step.")
  message("[Step 01] Step 02 will load all events without cell filtering.")
} else {

raw_dir  <- config$directories$raw
out_dir  <- config$directories$processed
logs_dir <- config$directories$logs
qc_cfg   <- config$qc

if (!dir.exists(out_dir))  dir.create(out_dir,  recursive = TRUE)
if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)

# Enumerate FCS files per batch, respecting test mode
all_files_df <- do.call(rbind, lapply(names(raw_dir), function(batch_name) {
  dir_path <- raw_dir[[batch_name]]
  if (!dir.exists(dir_path)) {
    warning(sprintf("[Step 01] Directory not found for batch '%s': %s", batch_name, dir_path),
            call. = FALSE)
    return(NULL)
  }
  files <- list.files(dir_path, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0L) {
    warning(sprintf("[Step 01] No FCS files found in: %s", dir_path), call. = FALSE)
    return(NULL)
  }
  is_test  <- isTRUE(config$testing$enabled)
  test_lim <- as.integer(config$testing$max_files_per_panel %||% 5L)
  if (is_test) files <- head(files, test_lim)
  data.frame(file_path = files, batch_id = batch_name, stringsAsFactors = FALSE)
}))

if (is.null(all_files_df) || nrow(all_files_df) == 0L) {
  stop("[Step 01 Fatal] No FCS files found across any batch directory.")
}

message(sprintf("[Step 01] Files to process: %d across %d batch(es)",
                nrow(all_files_df), length(unique(all_files_df$batch_id))))

log_obj <- init_step_log(
  step_name   = "01_qc",
  step_number = 1L,
  input_files = all_files_df$file_path
)

# Per-sample QC loop
qc_filters      <- vector("list", nrow(all_files_df))
names(qc_filters) <- basename(all_files_df$file_path)
qc_summary_list <- list()

for (i in seq_len(nrow(all_files_df))) {
  fp  <- all_files_df$file_path[i]
  sid <- basename(fp)
  bid <- all_files_df$batch_id[i]

  message(sprintf("[Step 01] [%d/%d] [%s] %s", i, nrow(all_files_df), bid, sid))

  result <- tryCatch(
    run_sample_qc(fp, qc_config = qc_cfg, sample_id = sid),
    error = function(e) {
      warning(sprintf("[Step 01] QC failed for '%s': %s — sample EXCLUDED.", sid, e$message),
              call. = FALSE)
      NULL
    }
  )

  if (is.null(result)) {
    qc_filters[[sid]] <- NULL
    qc_summary_list[[sid]] <- list(status = "EXCLUDED", n_raw = NA, n_final = NA)
    next
  }

  qc_filters[[sid]] <- list(
    valid_indices = result$valid_indices,
    batch_id      = bid,
    file_path     = fp
  )
  qc_summary_list[[sid]] <- result$qc_metrics

  message(sprintf("         n_raw=%d  n_final=%d  (%.1f%% removed)",
                  result$qc_metrics$n_raw,
                  result$qc_metrics$n_final,
                  result$qc_metrics$pct_total_removed))
}

n_processed      <- length(qc_filters)
n_excluded       <- nrow(all_files_df) - n_processed
passing_sids     <- names(qc_filters)
mean_pct_removed <- mean(
  sapply(qc_summary_list[passing_sids], function(x) x$pct_total_removed),
  na.rm = TRUE
)

message(sprintf("[Step 01] Samples processed: %d  |  Excluded: %d  |  Mean removal: %.1f%%",
                n_processed, n_excluded, mean_pct_removed))

out_file <- file.path(out_dir, "qc_cell_filters.rds")
saveRDS(qc_filters, out_file)
message(sprintf("[Step 01] Cell filters saved: %s", out_file))

add_metric(log_obj, "n_samples_total",     nrow(all_files_df))
add_metric(log_obj, "n_samples_processed", n_processed)
add_metric(log_obj, "n_samples_excluded",  n_excluded)
add_metric(log_obj, "mean_pct_removed",    round(mean_pct_removed, 2))
add_metric(log_obj, "per_sample_qc",       qc_summary_list)

finalize_step_log(log_obj, output_files = out_file, status = "SUCCESS")
write_step_json(log_obj, logs_dir)

message("\n=== STEP 1 COMPLETE ===\n")
} # end if (isTRUE(config$qc$enabled))
