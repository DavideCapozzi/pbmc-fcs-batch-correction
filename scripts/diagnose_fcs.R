#!/usr/bin/env Rscript
# scripts/diagnose_fcs.R
# ==============================================================================
# FCS DIAGNOSTIC SCRIPT
# Analyses all FCS files in the configured raw directories and produces a
# comprehensive JSON report to diagnose QC exclusion patterns.
#
# For each file the script runs five analysis blocks:
#   [A] FCS header keywords (instrument, acquisition timing, channel count)
#   [B] Per-channel statistics on raw events (mean, SD, saturation, negatives)
#   [C] Temporal analysis — event flow stability and per-channel drift scores
#   [D] PeacoQC IT_limit sweep — tests IT_limit 0.3 → 1.0 to find the minimum
#       value that keeps removal below the exclusion threshold
#   [E] Full gate cascade (PeacoQC → singlet → debris → viability) with
#       per-step cell counts and the final exclusion reason
#
# A comparison_summary section aggregates excluded vs passing files so you can
# immediately identify the discriminating metric.
#
# Usage:
#   Rscript scripts/diagnose_fcs.R
#   Rscript scripts/diagnose_fcs.R results/diagnostics/custom_name.json
#
# Output:
#   results/diagnostics/fcs_diagnostic_<timestamp>.json
# ==============================================================================

suppressPackageStartupMessages({
  library(flowCore)
  library(yaml)
  library(jsonlite)
  library(parallel)
})

source("src/functions/qc.R")

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

parse_acq_time <- function(t_str) {
  # Parses FCS $BTIM / $ETIM strings of the form "HH:MM:SS[.ss]" to seconds.
  if (is.null(t_str) || is.na(t_str) || !nzchar(t_str)) return(NA_real_)
  parts <- strsplit(trimws(t_str), ":")[[1]]
  if (length(parts) < 3) return(NA_real_)
  as.numeric(parts[1]) * 3600 + as.numeric(parts[2]) * 60 + as.numeric(parts[3])
}

safe_kw <- function(kw, key) {
  val <- kw[[key]]
  if (is.null(val) || length(val) == 0) NA_character_ else as.character(val)
}

# Computes mean over a numeric vector, guarding for all-NA input.
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# ------------------------------------------------------------------------------
# [A] Extract FCS header keywords
# ------------------------------------------------------------------------------

extract_keywords <- function(ff) {
  kw   <- flowCore::keyword(ff)
  btim <- parse_acq_time(safe_kw(kw, "$BTIM"))
  etim <- parse_acq_time(safe_kw(kw, "$ETIM"))
  dur  <- if (!is.na(btim) && !is.na(etim) && etim > btim) etim - btim else NA_real_

  n_par <- suppressWarnings(as.integer(safe_kw(kw, "$PAR")))
  n_tot <- suppressWarnings(as.integer(safe_kw(kw, "$TOT")))

  list(
    date                   = safe_kw(kw, "$DATE"),
    btim                   = safe_kw(kw, "$BTIM"),
    etim                   = safe_kw(kw, "$ETIM"),
    acquisition_duration_sec = if (is.na(dur)) NULL else round(dur, 1),
    cytometer              = safe_kw(kw, "$CYT"),
    tube_name              = safe_kw(kw, "$TUBE NAME"),
    src                    = safe_kw(kw, "$SRC"),
    vol_ul                 = safe_kw(kw, "$VOL"),
    n_parameters           = if (is.na(n_par)) NULL else n_par,
    n_events_declared      = if (is.na(n_tot)) NULL else n_tot
  )
}

# ------------------------------------------------------------------------------
# [B] Per-channel statistics on raw events
# ------------------------------------------------------------------------------

channel_stats <- function(ff) {
  mat     <- flowCore::exprs(ff)
  params  <- flowCore::pData(flowCore::parameters(ff))
  results <- list()

  for (ch in colnames(mat)) {
    vals     <- mat[, ch]
    rng_row  <- which(params$name == ch)
    rng_max  <- if (length(rng_row) > 0) suppressWarnings(as.numeric(params$range[rng_row[1]])) else NA_real_
    n_total  <- length(vals)
    n_finite <- sum(is.finite(vals))

    pct_above <- if (!is.na(rng_max) && rng_max > 0)
      round(100 * sum(vals >= rng_max, na.rm = TRUE) / n_total, 3) else NA_real_
    pct_neg   <- round(100 * sum(vals < 0, na.rm = TRUE) / n_total, 3)

    results[[ch]] <- list(
      mean           = round(safe_mean(vals), 2),
      sd             = round(sd(vals, na.rm = TRUE), 2),
      median         = round(median(vals, na.rm = TRUE), 2),
      iqr            = round(IQR(vals, na.rm = TRUE), 2),
      min            = round(min(vals, na.rm = TRUE), 2),
      max            = round(max(vals, na.rm = TRUE), 2),
      range_max      = if (is.na(rng_max)) NULL else rng_max,
      pct_above_range = if (is.na(pct_above)) NULL else pct_above,
      pct_negative   = pct_neg,
      n_finite       = n_finite
    )
  }
  results
}

# ------------------------------------------------------------------------------
# [C] Temporal analysis — event flow + per-channel drift
# ------------------------------------------------------------------------------

temporal_analysis <- function(ff, n_bins = 20L) {
  mat      <- flowCore::exprs(ff)
  n_events <- nrow(mat)

  # Locate the TIME channel (case-insensitive)
  time_col <- colnames(mat)[grepl("^time$", colnames(mat), ignore.case = TRUE)][1]

  if (!is.na(time_col) && !is.null(time_col)) {
    time_vals <- mat[, time_col]
    # Some cytometers store time in 10ms ticks; normalise to rank if the range
    # is tiny (< 1) so bin boundaries still make sense.
    if (all(is.na(time_vals)) || max(time_vals, na.rm = TRUE) < 2) time_vals <- seq_len(n_events)
  } else {
    time_vals <- seq_len(n_events)
  }

  breaks   <- seq(min(time_vals, na.rm = TRUE), max(time_vals, na.rm = TRUE),
                  length.out = n_bins + 1L)
  bin_idx  <- cut(time_vals, breaks = breaks, labels = FALSE, include.lowest = TRUE)

  # Event count per bin
  counts <- tabulate(bin_idx, nbins = n_bins)
  cv_events <- if (mean(counts) > 0) sd(counts) / mean(counts) else NA_real_
  max_fold  <- if (min(counts) > 0) max(counts) / min(counts) else NA_real_

  # Drift per biological (non-scatter) channel
  scatter_re <- "^FSC|^SSC|^Time$|Width|Length|Event"
  bio_cols   <- colnames(mat)[!grepl(scatter_re, colnames(mat), ignore.case = TRUE)]

  drift_scores <- list()
  for (ch in bio_cols) {
    vals    <- mat[, ch]
    sd_glob <- sd(vals, na.rm = TRUE)
    if (is.na(sd_glob) || sd_glob < 1e-9) next
    bin_means <- vapply(seq_len(n_bins), function(b) {
      safe_mean(vals[!is.na(bin_idx) & bin_idx == b])
    }, numeric(1))
    ref <- bin_means[1]
    score <- if (!is.na(ref)) max(abs(bin_means - ref), na.rm = TRUE) / sd_glob else NA_real_
    drift_scores[[ch]] <- round(score, 4)
  }

  valid_drift <- Filter(Negate(is.na), drift_scores)
  top_ch      <- if (length(valid_drift) > 0) names(which.max(unlist(valid_drift))) else NA_character_
  top_score   <- if (length(valid_drift) > 0) max(unlist(valid_drift), na.rm = TRUE) else NA_real_

  list(
    n_bins                = n_bins,
    events_per_bin        = as.integer(counts),
    cv_events_per_bin     = round(cv_events, 4),
    max_fold_change_events = round(max_fold, 3),
    channel_drift_scores  = drift_scores,
    highest_drift_channel = top_ch,
    highest_drift_score   = round(top_score, 4)
  )
}

# ------------------------------------------------------------------------------
# [D] PeacoQC IT_limit sweep
# ------------------------------------------------------------------------------

peacoqc_sweep <- function(ff, base_cfg, max_pct_removed = 85) {
  it_values   <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  results     <- list()
  min_passing <- NA_real_

  for (it in it_values) {
    cfg_i <- modifyList(base_cfg, list(IT_limit = it))
    res   <- tryCatch(
      apply_peacoqc_filter(ff, channels_to_use = NULL, peacoqc_params = cfg_i),
      error = function(e) list(pct_removed = NA_real_, n_after = NA_integer_)
    )
    pct    <- round(res$pct_removed, 3)
    surv   <- if (!is.null(res$n_after)) as.integer(res$n_after) else NA_integer_
    passes <- !is.na(pct) && pct < max_pct_removed

    key  <- paste0("it_", gsub("\\.", "_", as.character(it)))
    results[[key]] <- list(
      IT_limit    = it,
      pct_removed = pct,
      n_surviving = surv,
      passes      = passes
    )

    if (passes && is.na(min_passing)) min_passing <- it
  }

  valid_pcts <- sapply(results, function(r) r$pct_removed)
  valid_pcts <- as.numeric(valid_pcts[!is.na(valid_pcts)])
  it_removal_range_pct <- if (length(valid_pcts) >= 2) round(max(valid_pcts) - min(valid_pcts), 3) else NA_real_

  results$it_min_passing       <- if (is.na(min_passing)) NULL else min_passing
  results$it_removal_range_pct <- if (is.na(it_removal_range_pct)) NULL else it_removal_range_pct
  results$it_recoverable       <- isTRUE(!is.na(it_removal_range_pct) && it_removal_range_pct >= 10)
  results
}

# ------------------------------------------------------------------------------
# [E] Full gate cascade
# ------------------------------------------------------------------------------

gate_cascade <- function(file_path, qc_config, sample_id) {
  max_pct  <- as.numeric(qc_config$max_pct_removed  %||% 85)
  min_cells <- as.integer(qc_config$min_cells_final %||% 5000L)

  result <- tryCatch({
    qc_config_permissive <- modifyList(qc_config, list(max_pct_removed = 100))
    qc_res <- run_sample_qc(file_path, qc_config_permissive, sample_id)
    m      <- qc_res$qc_metrics
    n_fin  <- m$n_final
    pct_tot <- round(100 * (1 - n_fin / m$n_raw), 3)

    reason <- if (pct_tot > max_pct) {
      sprintf("max_pct_removed (%.1f%% > %.0f%%)", pct_tot, max_pct)
    } else if (n_fin < min_cells) {
      sprintf("min_cells_final (%d < %d)", n_fin, min_cells)
    } else {
      "passed"
    }

    gp <- m$gate_params %||% list()
    list(
      n_raw              = m$n_raw,
      n_after_peacoqc    = m$n_after_peacoqc,
      n_after_singlet    = m$n_after_singlet,
      n_after_debris     = m$n_after_debris,
      n_after_viability  = m$n_after_viability,
      n_final            = n_fin,
      pct_peacoqc_removed  = round(m$pct_peacoqc_removed   %||% NA_real_, 2),
      pct_singlet_removed  = round(m$pct_singlet_removed    %||% NA_real_, 2),
      pct_debris_removed   = round(m$pct_debris_removed     %||% NA_real_, 2),
      pct_viability_removed = round(m$pct_viability_removed %||% NA_real_, 2),
      pct_total_removed    = pct_tot,
      peacoqc_IT_limit_used = m$peacoqc_IT_limit_used %||% NA_real_,
      peacoqc_retry_triggered = isTRUE(m$peacoqc_retry_triggered),
      viability_channel    = m$viability_channel_raw %||% NA_character_,
      viability_threshold  = m$viability_threshold   %||% NA_real_,
      singlet_gate_center  = gp$singlet_center %||% NA_real_,
      singlet_gate_mad     = gp$singlet_mad    %||% NA_real_,
      singlet_gate_lower   = gp$singlet_lower  %||% NA_real_,
      singlet_gate_upper   = gp$singlet_upper  %||% NA_real_,
      singlet_efficiency_pct = if (!is.null(m$n_raw) && !is.na(m$n_raw) && m$n_raw > 0)
        round(100 * (m$n_after_singlet %||% m$n_raw) / m$n_raw, 2) else NA_real_,
      fsc_min_thresh       = gp$fsc_min_thresh %||% NA_real_,
      ssc_max_thresh       = gp$ssc_max_thresh %||% NA_real_,
      viability_suspicious = isTRUE(
        (m$pct_viability_removed %||% NA_real_) < 2 |
        (m$pct_viability_removed %||% NA_real_) > 30
      ),
      exclusion_reason     = reason
    )
  }, error = function(e) {
    list(
      n_raw = NA_integer_, n_final = NA_integer_,
      pct_total_removed = NA_real_,
      exclusion_reason  = paste("read_error:", conditionMessage(e))
    )
  })

  result
}

# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------

args       <- commandArgs(trailingOnly = TRUE)
config     <- yaml::read_yaml("config/global_params.yml")
qc_cfg     <- config$qc
max_pct    <- as.numeric(qc_cfg$max_pct_removed %||% 85)
peacoqc_base <- qc_cfg$peacoqc %||% list(
  MAD = 6, IT_limit = 0.6, remove_zeros = FALSE, min_cells = 150L
)

# Output path
diag_dir <- "results/diagnostics"
if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)
ts       <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_path <- if (length(args) > 0) args[1] else file.path(diag_dir, paste0("fcs_diagnostic_", ts, ".json"))

message("=== FCS DIAGNOSTIC SCRIPT ===")
message(sprintf("[Diag] Output: %s", out_path))

# Load previous QC outcomes if available
qc_filters_file <- file.path(config$directories$processed, "qc_cell_filters.rds")
prev_qc <- if (file.exists(qc_filters_file)) {
  message(sprintf("[Diag] Loading previous QC outcomes from %s", qc_filters_file))
  readRDS(qc_filters_file)
} else {
  message("[Diag] No qc_cell_filters.rds found — all files marked as 'unknown'")
  NULL
}

# Discover all FCS files
raw_dirs  <- config$directories$raw
all_files <- do.call(rbind, lapply(names(raw_dirs), function(batch) {
  fcs <- list.files(raw_dirs[[batch]], pattern = "\\.fcs$",
                    full.names = TRUE, recursive = FALSE)
  if (length(fcs) == 0L) return(NULL)
  data.frame(path = fcs, batch = batch, stringsAsFactors = FALSE)
}))

if (is.null(all_files) || nrow(all_files) == 0L) {
  message("[Diag] No FCS files found across any batch directory. Exiting.")
  quit(save = "no", status = 0)
}

message(sprintf("[Diag] Found %d FCS files across %d batches.",
                nrow(all_files), length(raw_dirs)))

# Determine previous QC outcome per file
get_prev_outcome <- function(fname) {
  if (is.null(prev_qc)) return("unknown")
  entry <- prev_qc[[fname]]
  if (is.null(entry)) return("not_in_last_run")
  if (!is.null(entry$status) && entry$status == "EXCLUDED") return("EXCLUDED")
  if (!is.null(entry$n_final) && !is.na(entry$n_final)) return("PASSED")
  return("unknown")
}

# Per-file analysis (parallelised)
n_cores <- min(detectCores(logical = FALSE), as.integer(config$n_cores %||% 4L))

# Progress log — workers write here; monitor with: tail -f <path>
prog_log <- file.path(diag_dir, paste0("progress_", ts, ".log"))
writeLines(sprintf("=== FCS Diagnostic progress log — %d files, %d cores ===",
                   nrow(all_files), n_cores), prog_log)
message(sprintf("[Diag] Parallelising over %d physical cores.", n_cores))
message(sprintf("[Diag] Follow progress in a second terminal:  tail -f %s", prog_log))

log_progress <- function(...) {
  line <- paste0(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))
  cat(line, "\n", file = prog_log, append = TRUE)
}

analyse_file <- function(i) {
  fpath  <- all_files$path[i]
  batch  <- all_files$batch[i]
  fname  <- basename(fpath)
  prev   <- get_prev_outcome(fname)
  n_tot  <- nrow(all_files)

  log_progress("[%d/%d] START  %s (%s)", i, n_tot, fname, batch)

  entry <- list(path = fpath, batch = batch, prev_qc_outcome = prev)

  ff <- tryCatch(
    flowCore::read.FCS(fpath, transformation = FALSE,
                       truncate_max_range = FALSE, emptyValue = FALSE),
    error = function(e) NULL
  )

  if (is.null(ff)) {
    entry$error      <- "Failed to read FCS file"
    entry$qc_outcome <- "read_error"
    log_progress("[%d/%d] ERROR  %s — could not read FCS", i, n_tot, fname)
    return(list(name = fname, entry = entry))
  }

  entry$fcs_keywords      <- extract_keywords(ff)
  entry$n_raw_events      <- nrow(flowCore::exprs(ff))
  entry$per_channel_stats <- channel_stats(ff)

  entry$temporal_analysis <- tryCatch(
    temporal_analysis(ff),
    error = function(e) list(error = conditionMessage(e))
  )

  # Suppress PeacoQC's per-file progress bars — they interleave badly across workers
  entry$peacoqc_sweep <- tryCatch(
    suppressMessages(peacoqc_sweep(ff, peacoqc_base, max_pct)),
    error = function(e) list(error = conditionMessage(e))
  )

  entry$gate_cascade <- suppressMessages(gate_cascade(fpath, qc_cfg, fname))

  cascade_reason   <- entry$gate_cascade$exclusion_reason %||% "unknown"
  entry$qc_outcome <- if (grepl("^passed", cascade_reason)) "PASSED" else "EXCLUDED"

  log_progress("[%d/%d] DONE   %s — %s (n_final=%s)",
               i, n_tot, fname, entry$qc_outcome,
               entry$gate_cascade$n_final %||% "NA")

  list(name = fname, entry = entry)
}

file_results_raw <- mclapply(
  seq_len(nrow(all_files)),
  analyse_file,
  mc.cores       = n_cores,
  mc.preschedule = FALSE  # dynamic scheduling — excluded files (faster) free cores early
)

# Collect into named list; NULL results indicate mclapply worker crash
file_results <- list()
for (res in file_results_raw) {
  if (!is.null(res) && !inherits(res, "try-error")) {
    file_results[[res$name]] <- res$entry
  }
}
log_progress("All %d files processed. Building summary...", nrow(all_files))

# ------------------------------------------------------------------------------
# Comparison summary
# ------------------------------------------------------------------------------

excluded_names <- names(Filter(function(x) identical(x$qc_outcome, "EXCLUDED"), file_results))
passing_names  <- names(Filter(function(x) identical(x$qc_outcome, "PASSED"),   file_results))

pull_metric <- function(fname, path_fn) {
  result <- tryCatch(path_fn(file_results[[fname]]), error = function(e) NA_real_)
  # NULL access (missing key in error-path gate_cascade) must become NA, not NULL,
  # so sapply always returns a vector rather than a list.
  if (is.null(result)) NA_real_ else result
}

mean_or_null <- function(vals) {
  if (is.list(vals)) vals <- unlist(vals)
  v <- as.numeric(vals[!is.na(vals)])
  if (length(v) == 0) NULL else round(mean(v), 3)
}

excl_pct_pqc <- sapply(excluded_names, pull_metric,
                       path_fn = function(x) x$gate_cascade$pct_peacoqc_removed)
pass_pct_pqc <- sapply(passing_names,  pull_metric,
                       path_fn = function(x) x$gate_cascade$pct_peacoqc_removed)

excl_cv  <- sapply(excluded_names, pull_metric,
                   path_fn = function(x) x$temporal_analysis$cv_events_per_bin)
pass_cv  <- sapply(passing_names,  pull_metric,
                   path_fn = function(x) x$temporal_analysis$cv_events_per_bin)

excl_drift <- sapply(excluded_names, pull_metric,
                     path_fn = function(x) x$temporal_analysis$highest_drift_score)
pass_drift <- sapply(passing_names,  pull_metric,
                     path_fn = function(x) x$temporal_analysis$highest_drift_score)

excl_dur <- sapply(excluded_names, pull_metric,
                   path_fn = function(x) x$fcs_keywords$acquisition_duration_sec)
pass_dur <- sapply(passing_names,  pull_metric,
                   path_fn = function(x) x$fcs_keywords$acquisition_duration_sec)

excl_sing <- sapply(excluded_names, pull_metric,
                    path_fn = function(x) x$gate_cascade$singlet_efficiency_pct)
pass_sing <- sapply(passing_names,  pull_metric,
                    path_fn = function(x) x$gate_cascade$singlet_efficiency_pct)

excl_via <- sapply(excluded_names, pull_metric,
                   path_fn = function(x) x$gate_cascade$pct_viability_removed)
pass_via <- sapply(passing_names,  pull_metric,
                   path_fn = function(x) x$gate_cascade$pct_viability_removed)

# Most common exclusion reason in excluded files
excl_reasons <- sapply(excluded_names,
                       function(n) file_results[[n]]$gate_cascade$exclusion_reason %||% "unknown")

# IT_limit min passing across excluded files (NULL means irrecoverable)
excl_it_min <- sapply(excluded_names, function(n) {
  sw <- file_results[[n]]$peacoqc_sweep
  v  <- sw$it_min_passing
  if (is.null(v)) NA_real_ else as.numeric(v)
})
n_irrecoverable <- sum(!sapply(excluded_names, function(n)
  isTRUE(file_results[[n]]$peacoqc_sweep$it_recoverable)))

# Per-batch summary — breakdown by batch for adequacy assessment
per_batch_summary <- lapply(names(raw_dirs), function(b) {
  batch_names <- basename(all_files$path[all_files$batch == b])
  n_total     <- length(batch_names)
  n_excl      <- sum(batch_names %in% excluded_names)
  n_pass      <- sum(batch_names %in% passing_names)
  pass_n      <- intersect(batch_names, passing_names)

  pass_cells <- sapply(pass_n, function(n)
    file_results[[n]]$gate_cascade$n_final %||% NA_real_)
  mean_cells <- mean_or_null(pass_cells)

  list(
    n_files              = n_total,
    n_excluded           = n_excl,
    n_passing            = n_pass,
    pass_rate            = if (n_total > 0) round(n_pass / n_total, 3) else NA_real_,
    mean_n_final_passing = mean_cells,
    adequate_for_baseline = isTRUE(n_pass >= 5 && !is.null(mean_cells) && mean_cells >= 5000),
    mean_pct_peacoqc_removed = mean_or_null(sapply(pass_n, function(n)
      file_results[[n]]$gate_cascade$pct_peacoqc_removed %||% NA_real_)),
    mean_singlet_efficiency_pct = mean_or_null(sapply(pass_n, function(n)
      file_results[[n]]$gate_cascade$singlet_efficiency_pct %||% NA_real_)),
    mean_pct_viability_removed = mean_or_null(sapply(pass_n, function(n)
      file_results[[n]]$gate_cascade$pct_viability_removed %||% NA_real_))
  )
})
names(per_batch_summary) <- names(raw_dirs)

# Statistical adequacy projection — answers "can we get defensible ±2SD ranges?"
baseline_adequacy_projection <- list(
  note = paste(
    "n_passing >= 5 per batch and mean_n_final > 5000 are the minimum criteria",
    "for statistically defensible ±2SD clinical flags.",
    "adequate_for_baseline = FALSE means the baseline will have wide CIs."
  ),
  per_batch = lapply(per_batch_summary, function(b) list(
    n_passing             = b$n_passing,
    mean_cells_per_sample = b$mean_n_final_passing,
    estimated_total_cells = if (!is.null(b$mean_n_final_passing))
      as.integer(round(b$n_passing * b$mean_n_final_passing)) else NULL,
    adequate_for_baseline = b$adequate_for_baseline
  ))
)

comparison_summary <- list(
  n_files_total   = nrow(all_files),
  n_excluded      = length(excluded_names),
  n_passing       = length(passing_names),
  excluded_files  = excluded_names,
  passing_files   = passing_names,
  exclusion_reasons = as.list(excl_reasons),
  n_irrecoverable_by_peacoqc = n_irrecoverable,
  per_batch_summary           = per_batch_summary,
  baseline_adequacy_projection = baseline_adequacy_projection,
  metrics_comparison = list(
    mean_pct_peacoqc_removed = list(
      excluded = mean_or_null(excl_pct_pqc),
      passing  = mean_or_null(pass_pct_pqc)
    ),
    mean_cv_events_per_bin = list(
      excluded = mean_or_null(excl_cv),
      passing  = mean_or_null(pass_cv)
    ),
    mean_highest_drift_score = list(
      excluded = mean_or_null(excl_drift),
      passing  = mean_or_null(pass_drift)
    ),
    mean_acquisition_duration_sec = list(
      excluded = mean_or_null(excl_dur),
      passing  = mean_or_null(pass_dur)
    ),
    mean_singlet_efficiency_pct = list(
      excluded = mean_or_null(excl_sing),
      passing  = mean_or_null(pass_sing)
    ),
    mean_pct_viability_removed = list(
      excluded = mean_or_null(excl_via),
      passing  = mean_or_null(pass_via)
    )
  )
)

# Build final report
report <- list(
  report_generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  config = list(
    IT_limit          = as.numeric(peacoqc_base$IT_limit %||% 0.6),
    MAD               = as.numeric(peacoqc_base$MAD %||% 6),
    max_pct_removed   = max_pct,
    min_cells_final   = as.integer(qc_cfg$min_cells_final %||% 5000L),
    batches           = names(raw_dirs)
  ),
  comparison_summary = comparison_summary,
  files = file_results
)

jsonlite::write_json(report, out_path, pretty = TRUE, auto_unbox = TRUE, na = "null")
message(sprintf("\n[Diag] Done. Report written to: %s", out_path))
message(sprintf("[Diag] %d excluded, %d passing, %d irrecoverable by PeacoQC.",
                length(excluded_names), length(passing_names), n_irrecoverable))
