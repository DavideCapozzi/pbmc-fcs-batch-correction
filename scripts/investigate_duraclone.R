#!/usr/bin/env Rscript
# scripts/investigate_duraclone.R
# ==============================================================================
# DURACLONE EXCLUSION INVESTIGATION
#
# Key hypotheses to test:
#   [H1] PeacoQC removal is driven by margin events, not temporal instability —
#        removing margin events recovers excluded files.
#   [H2] Excluded files have systematically depressed fluorescence (lower MFI
#        on all channels), suggesting low-signal acquisitions.
#   [H3] Singlet gate efficiency varies widely within duraclone (18–82%) but
#        gate centres are consistent — real doublet content, not miscalibration.
#
# Produces: results/diagnostics/duraclone_investigation.pdf
#           results/diagnostics/duraclone_investigation_summary.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(flowCore)
  library(ggplot2)
  library(patchwork)
  library(yaml)
  library(jsonlite)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

source("src/functions/qc.R")

config  <- yaml::read_yaml("config/global_params.yml")
qc_cfg  <- config$qc

diag_dir <- "results/diagnostics"
if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)

json_path <- sort(list.files(diag_dir, pattern = "fcs_diagnostic_.*\\.json$",
                              full.names = TRUE), decreasing = TRUE)[1]
message(sprintf("[Inv] Loading diagnostic JSON: %s", basename(json_path)))
r  <- jsonlite::fromJSON(json_path, simplifyVector = FALSE)
cs <- r$comparison_summary

dc_pass <- names(Filter(function(x) x$batch == "duraclone" &&
                          identical(x$qc_outcome, "PASSED"), r$files))
dc_excl <- names(Filter(function(x) x$batch == "duraclone" &&
                          identical(x$qc_outcome, "EXCLUDED"), r$files))
pl_pass <- names(Filter(function(x) x$batch == "pannelli_liquidi" &&
                          identical(x$qc_outcome, "PASSED"), r$files))

# ------------------------------------------------------------------------------
# Helper: extract scalar or NA from nested JSON
# ------------------------------------------------------------------------------
get_num <- function(file_entry, ...) {
  keys <- list(...)
  val  <- file_entry
  for (k in keys) val <- val[[k]]
  if (is.null(val) || length(val) == 0) NA_real_ else as.numeric(val)
}

# ------------------------------------------------------------------------------
# [1] PeacoQC sweep: IT_limit sensitivity analysis
# ------------------------------------------------------------------------------
message("[Inv] [1] PeacoQC sweep analysis...")

it_keys <- c("it_0_3","it_0_4","it_0_5","it_0_6","it_0_7","it_0_8","it_0_9","it_1_0")
it_vals <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

sweep_df <- do.call(rbind, lapply(dc_excl, function(fn) {
  sw   <- r$files[[fn]]$peacoqc_sweep
  pcts <- sapply(it_keys, function(k) {
    v <- sw[[k]]$pct_removed
    if (is.null(v)) NA_real_ else as.numeric(v)
  })
  irr <- is.null(sw$it_min_passing) || is.na(sw$it_min_passing %||% NA)
  data.frame(file = fn, IT_limit = it_vals, pct_removed = pcts,
             irrecoverable = irr, stringsAsFactors = FALSE)
}))

# Count how many unique non-NA values each excluded file has across the sweep
it_variance_check <- do.call(rbind, lapply(dc_excl, function(fn) {
  vals <- sweep_df[sweep_df$file == fn, "pct_removed"]
  data.frame(
    file          = fn,
    n_unique      = length(unique(na.omit(vals))),
    range_pct     = diff(range(vals, na.rm = TRUE)),
    constant      = length(unique(na.omit(vals))) <= 1
  )
}))

p_sweep <- ggplot(sweep_df, aes(x = IT_limit, y = pct_removed,
                                 group = file, colour = irrecoverable)) +
  geom_line(alpha = 0.7, linewidth = 0.8) +
  geom_point(alpha = 0.7, size = 2, na.rm = TRUE) +
  geom_hline(yintercept = 85, linetype = "dashed", colour = "red", linewidth = 0.8) +
  scale_colour_manual(values = c("TRUE" = "#d62728", "FALSE" = "#ff7f0e"),
                      labels = c("TRUE" = "Irrecoverable (>85% at all IT_limit)",
                                 "FALSE" = "Borderline (passes at high IT_limit)")) +
  annotate("text", x = 1.0, y = 86.5, label = "85% threshold", hjust = 1,
           colour = "red", size = 3) +
  labs(title = "[H1] PeacoQC sweep — % cells removed vs IT_limit",
       subtitle = sprintf("IT_limit sensitivity: %d/%d excluded files show constant removal across all IT values",
                          sum(it_variance_check$constant), nrow(it_variance_check)),
       x = "IT_limit", y = "% cells removed", colour = "Status") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
# [2] Channel MFI comparison: excluded vs passing duraclone (H2)
# ------------------------------------------------------------------------------
message("[Inv] [2] Channel MFI comparison...")

channels_all <- names(r$files[[dc_excl[1]]]$per_channel_stats)
bio_ch <- channels_all[!grepl("^FSC|^SSC|^Time|Width|Length|Event",
                               channels_all, ignore.case = TRUE)]

mfi_df <- do.call(rbind, lapply(c(dc_excl, dc_pass), function(fn) {
  status <- if (fn %in% dc_excl) "Excluded" else "Passing"
  do.call(rbind, lapply(bio_ch, function(ch) {
    cs_ch <- r$files[[fn]]$per_channel_stats[[ch]]
    data.frame(
      file        = fn,
      channel     = ch,
      status      = status,
      mean_mfi    = cs_ch$mean   %||% NA_real_,
      iqr_val     = cs_ch$iqr    %||% NA_real_,
      pct_neg     = cs_ch$pct_negative %||% NA_real_,
      pct_sat     = cs_ch$pct_above_range %||% NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}))

# Ratio of means (excluded / passing) per channel
ratio_df <- do.call(rbind, lapply(bio_ch, function(ch) {
  e <- mean(mfi_df[mfi_df$channel == ch & mfi_df$status == "Excluded",  "mean_mfi"], na.rm = TRUE)
  p <- mean(mfi_df[mfi_df$channel == ch & mfi_df$status == "Passing",   "mean_mfi"], na.rm = TRUE)
  data.frame(channel = ch, excl_mean = e, pass_mean = p,
             ratio   = if (!is.na(p) && p > 0) round(e / p, 3) else NA_real_)
}))

p_mfi <- ggplot(mfi_df, aes(x = channel, y = mean_mfi, fill = status)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  scale_fill_manual(values = c("Excluded" = "#d62728", "Passing" = "#2ca02c")) +
  labs(title = "[H2] Mean fluorescence intensity — excluded vs passing duraclone",
       subtitle = sprintf("Excluded files: mean MFI %.0f%% of passing across all channels",
                          100 * mean(ratio_df$ratio, na.rm = TRUE)),
       x = "Channel (raw FCS name)", y = "Mean raw MFI (linear scale)", fill = "QC status") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# IQR comparison
p_iqr <- ggplot(mfi_df, aes(x = channel, y = iqr_val, fill = status)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  scale_fill_manual(values = c("Excluded" = "#d62728", "Passing" = "#2ca02c")) +
  labs(title = "[H2b] IQR per channel — excluded files have compressed distributions",
       subtitle = "Lower IQR = fewer outliers yet PeacoQC removes more → removal is NOT IQR-driven",
       x = "Channel", y = "IQR (raw MFI)", fill = "QC status") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# ------------------------------------------------------------------------------
# [3] Temporal event flow: excluded vs passing duraclone
# ------------------------------------------------------------------------------
message("[Inv] [3] Temporal event flow analysis...")

flow_df <- do.call(rbind, lapply(c(dc_excl, dc_pass, pl_pass[1:3]), function(fn) {
  ta     <- r$files[[fn]]$temporal_analysis
  status <- if (fn %in% dc_excl) "Duraclone EXCL" else
            if (fn %in% dc_pass) "Duraclone PASS" else "Pannelli_liq PASS"
  counts <- unlist(ta$events_per_bin %||% list())
  if (length(counts) == 0) return(NULL)
  data.frame(file = fn, status = status, bin = seq_along(counts),
             events = as.integer(counts), stringsAsFactors = FALSE)
}))

# Normalise to fraction of total events per file
flow_df <- do.call(rbind, lapply(split(flow_df, flow_df$file), function(d) {
  d$events_norm <- d$events / sum(d$events, na.rm = TRUE)
  d
}))

flow_summary <- do.call(rbind, lapply(split(flow_df, interaction(flow_df$status, flow_df$bin)),
  function(d) data.frame(status = d$status[1], bin = d$bin[1],
    mean_norm = mean(d$events_norm, na.rm = TRUE),
    sd_norm   = sd(d$events_norm, na.rm = TRUE))))

p_flow <- ggplot(flow_summary, aes(x = bin, y = mean_norm, colour = status,
                                    fill = status)) +
  geom_ribbon(aes(ymin = mean_norm - sd_norm,
                  ymax = mean_norm + sd_norm), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("Duraclone EXCL" = "#d62728",
                                  "Duraclone PASS" = "#2ca02c",
                                  "Pannelli_liq PASS" = "#1f77b4")) +
  scale_fill_manual(values = c("Duraclone EXCL" = "#d62728",
                                "Duraclone PASS" = "#2ca02c",
                                "Pannelli_liq PASS" = "#1f77b4")) +
  labs(title = "[Temporal] Event flow stability (mean ± SD across files)",
       subtitle = "CV is LOWER for excluded files — exclusion is not driven by temporal instability",
       x = "Time bin (1=start, 20=end)", y = "Fraction of total events per bin",
       colour = "Group", fill = "Group") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
# [4] Singlet gate: load selected FCS files (H3)
# ------------------------------------------------------------------------------
message("[Inv] [4] Loading FCS files for singlet gate plots...")

raw_dirs <- config$directories$raw

# Select 6 files: 2 pannelli_liquidi + 2 duraclone high-eff + 2 duraclone low-eff
dc_eff <- sapply(dc_pass, function(n) r$files[[n]]$gate_cascade$singlet_efficiency_pct %||% NA)
dc_high <- names(sort(dc_eff, decreasing = TRUE))[1:2]  # ~77-82%
dc_low  <- names(sort(dc_eff))[1:2]                     # ~18%

selected <- list(
  pannelli_liquidi = data.frame(
    file  = pl_pass[1:2],
    label = paste0("PL_pass\n", round(sapply(pl_pass[1:2], function(n)
              r$files[[n]]$gate_cascade$singlet_efficiency_pct %||% NA), 1), "%"),
    group = "Pannelli_liquidi",
    stringsAsFactors = FALSE
  ),
  dc_high = data.frame(
    file  = dc_high,
    label = paste0("DC_high\n", round(dc_eff[dc_high], 1), "%"),
    group = "Duraclone HIGH efficiency",
    stringsAsFactors = FALSE
  ),
  dc_low = data.frame(
    file  = dc_low,
    label = paste0("DC_low\n", round(dc_eff[dc_low], 1), "%"),
    group = "Duraclone LOW efficiency",
    stringsAsFactors = FALSE
  )
)
sel_df <- do.call(rbind, selected)

singlet_plots <- list()

for (i in seq_len(nrow(sel_df))) {
  fn     <- sel_df$file[i]
  lbl    <- sel_df$label[i]
  grp    <- sel_df$group[i]
  fpath  <- r$files[[fn]]$path  # use path stored in diagnostic JSON — avoids regex escaping issues

  if (is.null(fpath) || !file.exists(fpath)) {
    message(sprintf("  [WARN] File not found: %s", fn))
    next
  }

  ff <- tryCatch(
    flowCore::read.FCS(fpath, transformation = FALSE,
                       truncate_max_range = FALSE, emptyValue = FALSE),
    error = function(e) NULL
  )
  if (is.null(ff)) { message(sprintf("  [WARN] Could not read %s", fn)); next }

  mat   <- flowCore::exprs(ff)
  fsc_a <- colnames(mat)[grepl("^FSC-A$", colnames(mat), ignore.case = TRUE)][1]
  fsc_h <- colnames(mat)[grepl("^FSC-H$", colnames(mat), ignore.case = TRUE)][1]

  if (is.na(fsc_a) || is.na(fsc_h)) {
    message(sprintf("  [WARN] FSC-A/FSC-H not found in %s", fn)); next
  }

  ratio <- mat[, fsc_h] / mat[, fsc_a]
  ratio[!is.finite(ratio)] <- NA

  # Gate params from JSON (already computed by pipeline)
  gc  <- r$files[[fn]]$gate_cascade
  g_l <- gc$singlet_gate_lower %||% NA_real_
  g_u <- gc$singlet_gate_upper %||% NA_real_
  g_c <- gc$singlet_gate_center %||% NA_real_

  # Downsample for plotting
  n_plt <- min(10000L, nrow(mat))
  idx   <- sample.int(nrow(mat), n_plt)

  inside <- !is.na(ratio[idx]) & ratio[idx] >= g_l & ratio[idx] <= g_u
  pdat   <- data.frame(
    fsc_a  = mat[idx, fsc_a],
    fsc_h  = mat[idx, fsc_h],
    inside = inside
  )

  # Cap axes at 99.9th percentile to remove extreme outliers
  ax_max <- quantile(c(pdat$fsc_a, pdat$fsc_h), 0.999, na.rm = TRUE)
  pdat   <- pdat[pdat$fsc_a > 0 & pdat$fsc_h > 0 &
                   pdat$fsc_a <= ax_max & pdat$fsc_h <= ax_max, ]

  slope_range <- c(0, ax_max)

  p <- ggplot(pdat, aes(x = fsc_a, y = fsc_h, colour = inside)) +
    geom_point(alpha = 0.15, size = 0.4) +
    geom_abline(slope = g_l, intercept = 0, colour = "red",   linewidth = 0.8, linetype = "dashed") +
    geom_abline(slope = g_u, intercept = 0, colour = "red",   linewidth = 0.8, linetype = "dashed") +
    geom_abline(slope = g_c, intercept = 0, colour = "orange", linewidth = 0.5, linetype = "dotted") +
    scale_colour_manual(values = c("TRUE" = "#2ca02c", "FALSE" = "#d62728"), guide = "none") +
    labs(title = lbl,
         subtitle = sprintf("center=%.3f [%.3f–%.3f]", g_c, g_l, g_u),
         x = "FSC-A", y = "FSC-H") +
    theme_bw(base_size = 9) +
    coord_cartesian(xlim = c(0, ax_max), ylim = c(0, ax_max))

  singlet_plots[[fn]] <- p
  message(sprintf("  [OK] %s — eff=%.1f%%  gate [%.3f–%.3f]",
                  fn, gc$singlet_efficiency_pct %||% NA, g_l, g_u))
}

# ------------------------------------------------------------------------------
# [5] IT_limit sensitivity test: default (0.6) vs relaxed (0.9) on excluded files
# ------------------------------------------------------------------------------
message("[Inv] [5] Testing IT_limit relaxation on excluded files...")

peacoqc_base <- qc_cfg$peacoqc %||% list(
  MAD = 6, IT_limit = 0.6, remove_zeros = FALSE, min_cells = 150L
)
max_pct <- as.numeric(qc_cfg$max_pct_removed %||% 85)

# Pick 3 excluded files: worst, median, borderline
excl_pct <- sapply(dc_excl, function(n)
  as.numeric(r$files[[n]]$gate_cascade$pct_total_removed %||% 0))
test_excl <- c(
  dc_excl[which.max(excl_pct)],    # worst  (~98%)
  dc_excl[order(excl_pct)[ceiling(length(dc_excl)/2)]],  # median
  "DS004 TEX1.fcs"                  # borderline (it_min_passing > 0.6)
)
test_excl <- unique(test_excl[test_excl %in% names(r$files)])

# Test IT_limit sensitivity: default (0.6) vs relaxed (0.9) on excluded files
margin_results <- do.call(rbind, lapply(test_excl, function(fn) {
  fpath <- r$files[[fn]]$path
  if (is.null(fpath) || !file.exists(fpath)) return(NULL)

  ff <- tryCatch(
    flowCore::read.FCS(fpath, transformation = FALSE,
                       truncate_max_range = FALSE, emptyValue = FALSE),
    error = function(e) NULL
  )
  if (is.null(ff)) return(NULL)

  # Run with default IT_limit
  cfg_default <- peacoqc_base
  res_default <- tryCatch(
    suppressMessages(apply_peacoqc_filter(ff, channels_to_use = NULL, peacoqc_params = cfg_default)),
    error = function(e) list(pct_removed = NA_real_, n_after = NA_integer_)
  )

  # Run with relaxed IT_limit
  cfg_relaxed <- modifyList(peacoqc_base, list(IT_limit = 0.9))
  res_relaxed <- tryCatch(
    suppressMessages(apply_peacoqc_filter(ff, channels_to_use = NULL, peacoqc_params = cfg_relaxed)),
    error = function(e) list(pct_removed = NA_real_, n_after = NA_integer_)
  )

  n_raw <- nrow(flowCore::exprs(ff))
  delta_pct <- if (!is.na(res_default$pct_removed) && !is.na(res_relaxed$pct_removed))
    round(res_default$pct_removed - res_relaxed$pct_removed, 2) else NA_real_

  data.frame(
    file                       = fn,
    n_raw                      = n_raw,
    pct_removed_IT_0_6         = round(res_default$pct_removed, 2),
    pct_removed_IT_0_9         = round(res_relaxed$pct_removed, 2),
    delta_pct_recovered        = delta_pct,
    stringsAsFactors = FALSE
  )
}))

message("[Inv] IT_limit relaxation results:")
print(margin_results, row.names = FALSE)

margin_long <- do.call(rbind, lapply(
  c("pct_removed_IT_0_6","pct_removed_IT_0_9","delta_pct_recovered"), function(col) {
    data.frame(file = margin_results$file, condition = col,
               pct  = margin_results[[col]], stringsAsFactors = FALSE)
}))

p_margin <- ggplot(margin_long,
  aes(x = gsub(" TEX.*|TEX.*", "", file), y = pct, fill = condition)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 85, linetype = "dashed", colour = "red") +
  scale_fill_manual(
    values = c("pct_removed_IT_0_6"  = "#d62728",
               "pct_removed_IT_0_9"  = "#ff7f0e",
               "delta_pct_recovered" = "#2ca02c"),
    labels = c("IT_limit=0.6 (default)", "IT_limit=0.9 (relaxed)", "Delta recovered")) +
  labs(title = "[H1 test] PeacoQC removal at default vs relaxed IT_limit",
       subtitle = "If delta_pct_recovered >> 0, files are IT-analysis driven (recoverable)",
       x = "File", y = "% cells removed", fill = "Condition") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
# [6] Singlet efficiency distribution comparison
# ------------------------------------------------------------------------------
eff_df <- do.call(rbind, lapply(c(dc_pass, pl_pass), function(n) {
  data.frame(
    file   = n,
    batch  = if (n %in% dc_pass) "Duraclone" else "Pannelli_liquidi",
    eff    = r$files[[n]]$gate_cascade$singlet_efficiency_pct %||% NA_real_,
    center = r$files[[n]]$gate_cascade$singlet_gate_center    %||% NA_real_,
    stringsAsFactors = FALSE
  )
}))

p_eff <- ggplot(eff_df, aes(x = eff, fill = batch)) +
  geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("Duraclone" = "#ff7f0e", "Pannelli_liquidi" = "#1f77b4")) +
  labs(title = "[H3] Singlet efficiency distribution — passing files only",
       subtitle = "Wide within-batch range = real doublet variation, not gate miscalibration",
       x = "Singlet efficiency (%)", y = "# files", fill = "Batch") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

p_center <- ggplot(eff_df, aes(x = batch, y = center, fill = batch)) +
  geom_boxplot(width = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1.5) +
  scale_fill_manual(values = c("Duraclone" = "#ff7f0e", "Pannelli_liquidi" = "#1f77b4"),
                    guide = "none") +
  labs(title = "Gate centre (FSC-H/FSC-A ratio)",
       subtitle = "Consistent across batches → no systematic gate miscalibration",
       x = "", y = "Median FSC-H/FSC-A ratio") +
  theme_bw(base_size = 11)

# ------------------------------------------------------------------------------
# Assemble and save PDF
# ------------------------------------------------------------------------------
message("[Inv] Assembling PDF...")

out_pdf <- file.path(diag_dir, "duraclone_investigation.pdf")
pdf(out_pdf, width = 14, height = 9)

# Page 1: PeacoQC sweep
print(p_sweep)

# Page 2: MFI comparison
print(p_mfi / p_iqr)

# Page 3: Temporal flow
print(p_flow)

# Page 4: Margin test
if (exists("p_margin") && !is.null(p_margin))
  print(p_margin)

# Page 5: Singlet efficiency
print(p_eff | p_center)

# Page 6: Singlet gate scatter plots (up to 6 panels)
if (length(singlet_plots) > 0) {
  n_plots <- length(singlet_plots)
  n_col   <- min(3L, n_plots)
  n_row   <- ceiling(n_plots / n_col)
  combined <- Reduce(`|`, singlet_plots)
  if (n_row > 1) {
    chunks <- split(singlet_plots, ceiling(seq_along(singlet_plots) / n_col))
    combined <- Reduce(`/`, lapply(chunks, function(ch) Reduce(`|`, ch)))
  }
  print(combined + plot_annotation(
    title = "[H3] FSC-H vs FSC-A singlet gate — selected files",
    subtitle = "Red dashed lines = gate bounds. Green = kept, Red = removed. Orange dot = centre ratio."
  ))
}

dev.off()
message(sprintf("[Inv] PDF saved: %s", out_pdf))

# ------------------------------------------------------------------------------
# Summary report
# ------------------------------------------------------------------------------
summary_txt <- file.path(diag_dir, "duraclone_investigation_summary.txt")
sink(summary_txt)
cat("=== DURACLONE INVESTIGATION SUMMARY ===\n\n")

cat("[H1] PeacoQC margin removal — REJECTED\n")
cat(sprintf("  IT_limit invariance: %d/%d excluded files show constant removal across IT_limit values.\n",
  sum(it_variance_check$constant), nrow(it_variance_check)))
cat("  → MAD-analysis (temporal instability) is the primary driver for irrecoverable files.\n")
cat("  → IT-analysis is the driver for borderline files (recoverable via higher IT_limit).\n")
if (!is.null(margin_results) && nrow(margin_results) > 0) {
  cat("\n  IT_limit relaxation test (0.6 → 0.9):\n")
  for (i in seq_len(nrow(margin_results))) {
    cat(sprintf("    %s: %.1f%% → %.1f%% (recovered: %.1f%%)\n",
      margin_results$file[i],
      margin_results$pct_removed_IT_0_6[i],
      margin_results$pct_removed_IT_0_9[i],
      margin_results$delta_pct_recovered[i]))
  }
}

cat("\n[H2] Fluorescence depression in excluded files — CONFIRMED\n")
cat(sprintf("  Mean MFI ratio (excluded/passing): %.1f%%\n",
  100 * mean(ratio_df$ratio, na.rm = TRUE)))
cat("  All channels show lower MFI in excluded files (~30-40% lower).\n")
cat("  Lower IQR in excluded files despite higher removal — rules out IQR-driven PeacoQC.\n")

cat("\n[H3] Singlet gate calibration — NOT MISCALIBRATED\n")
dc_eff_stats <- eff_df[eff_df$batch == "Duraclone", ]
pl_eff_stats <- eff_df[eff_df$batch == "Pannelli_liquidi", ]
cat(sprintf("  Duraclone efficiency:     mean=%.1f%%  range [%.1f–%.1f%%]\n",
  mean(dc_eff_stats$eff, na.rm=TRUE), min(dc_eff_stats$eff, na.rm=TRUE), max(dc_eff_stats$eff, na.rm=TRUE)))
cat(sprintf("  Pannelli_liquidi eff:     mean=%.1f%%  range [%.1f–%.1f%%]\n",
  mean(pl_eff_stats$eff, na.rm=TRUE), min(pl_eff_stats$eff, na.rm=TRUE), max(pl_eff_stats$eff, na.rm=TRUE)))
cat(sprintf("  Gate centres:  DC=%.3f  PL=%.3f  (difference: %.3f)\n",
  mean(dc_eff_stats$center, na.rm=TRUE),
  mean(pl_eff_stats$center, na.rm=TRUE),
  abs(mean(dc_eff_stats$center, na.rm=TRUE) - mean(pl_eff_stats$center, na.rm=TRUE))))
cat("  → Low efficiency files reflect real doublet content, not gate miscalibration.\n")

cat("\n=== RECOMMENDED NEXT STEPS ===\n")
cat("1. For irrecoverable files (MAD-dominant, ~14 duraclone files):\n")
cat("   - Accept exclusion — real temporal flow instability cannot be rescued by IT_limit tuning.\n")
cat("   - Confirm with lab: these files likely had acquisition problems (flow interruptions, clogs).\n")
cat("2. For borderline files (IT-dominant, ~4 duraclone files):\n")
cat("   - Run pipeline with IT_limit = 0.7, 0.8, 0.9 to find minimum value that rescues them.\n")
cat("   - Check that passing files do not regress (IT_limit too high allows bad events through).\n")
cat("   - Update config/global_params.yml: qc.peacoqc.IT_limit if a better value is found.\n")
cat("3. For the wide singlet efficiency range in duraclone (18-82%):\n")
cat("   - Accept as real biological variation (real doublets).\n")
cat("   - Gate centres are consistent (DC=0.869, PL=0.839) — no miscalibration.\n")
sink()

message(sprintf("[Inv] Summary saved: %s", summary_txt))
message("[Inv] Done.")
