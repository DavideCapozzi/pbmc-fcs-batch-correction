#!/usr/bin/env Rscript
# scripts/diagnose_integration.R
# ==============================================================================
# POST-PIPELINE INTEGRATION QUALITY DIAGNOSTIC
# Reads processed pipeline artifacts — no pipeline re-execution needed.
#
# Outputs (results/diagnostics/):
#   integration_diagnostic_{ts}.pdf   — 6-page report
#   integration_diagnostic_{ts}.json  — structured summary + recommendations
#
# Usage (from project root):
#   Rscript scripts/diagnose_integration.R
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(readxl)
  library(jsonlite)
  library(yaml)
  library(ggrepel)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b
safe_div <- function(num, den, default = NA_real_, eps = 1e-10)
  ifelse(abs(den) > eps, num / den, default)

# ── 0. Setup ──────────────────────────────────────────────────────────────────

config    <- yaml::read_yaml("config/global_params.yml")
ts        <- format(Sys.time(), "%Y%m%d_%H%M%S")
proc_dir  <- config$directories$processed   %||% "results/processed"
integ_dir <- config$directories$integration %||% "results/integration"
out_dir   <- "results/diagnostics"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_pdf  <- file.path(out_dir, paste0("integration_diagnostic_", ts, ".pdf"))
out_json <- file.path(out_dir, paste0("integration_diagnostic_", ts, ".json"))

message("=== INTEGRATION DIAGNOSTIC ===")

# ── 1. Load artifacts ─────────────────────────────────────────────────────────

mt     <- readRDS(file.path(proc_dir, "cluster_match_table.rds"))
bd     <- readRDS(file.path(proc_dir, "stratified_baseline_dictionary.rds"))
cent   <- readRDS(file.path(proc_dir, "cluster_centroids.rds"))
fm_obj <- readRDS(file.path(proc_dir, "frequency_matrix.rds"))
labs   <- readRDS(file.path(proc_dir, "auto_population_labels.rds"))

fm <- if (is.list(fm_obj) && !is.data.frame(fm_obj)) fm_obj$freq_matrix else fm_obj

val <- tryCatch(
  readxl::read_xlsx(file.path(integ_dir, "integration_validation_report.xlsx"),
                    sheet = "ClusterMatch"),
  error = function(e) { warning("[Diag] Validation report unreadable: ", conditionMessage(e)); NULL }
)

mt$pop_id <- paste0("pop_", mt$matched_population_id)
if (!is.null(val)) val$pop_id <- paste0("pop_", val$matched_population_id)

matched_pops <- sort(unique(mt$pop_id[mt$is_matched]))
batch_a      <- unique(mt$batch_a[mt$is_matched])[1]
batch_b      <- unique(mt$batch_b[mt$is_matched])[1]
markers      <- colnames(cent[[1]])

message(sprintf("[Diag] %d populations | batches: %s / %s | markers: %s",
  length(matched_pops), batch_a, batch_b, paste(markers, collapse = " ")))

# ── 2. Per-population summary ─────────────────────────────────────────────────

pop_summary <- do.call(rbind, lapply(matched_pops, function(p) {
  r  <- mt[mt$pop_id == p & mt$is_matched, ]
  v  <- if (!is.null(val)) val[val$pop_id == p, ] else data.frame()
  ba <- bd[bd$batch == batch_a & bd$population == p, ]
  bb <- bd[bd$batch == batch_b & bd$population == p, ]

  a_freq <- if (nrow(ba) > 0) ba$mean_freq else NA_real_
  a_sd   <- if (nrow(ba) > 0) ba$sd_freq   else NA_real_
  a_n    <- if (nrow(ba) > 0) ba$n_samples  else NA_integer_
  b_freq <- if (nrow(bb) > 0) bb$mean_freq  else NA_real_
  b_sd   <- if (nrow(bb) > 0) bb$sd_freq    else NA_real_
  b_n    <- if (nrow(bb) > 0) bb$n_samples  else NA_integer_

  a_cv  <- safe_div(a_sd, a_freq)
  b_cv  <- safe_div(b_sd, b_freq)
  ratio <- if (!is.na(a_freq) && !is.na(b_freq) && min(a_freq, b_freq) > 1e-4)
    max(a_freq, b_freq) / min(a_freq, b_freq) else NA_real_

  data.frame(
    pop_id   = p,
    label    = if (p %in% names(labs)) labs[[p]] else p,
    cos      = r$cosine_similarity[[1]],
    clust_a  = r$cluster_a[[1]],
    clust_b  = r$cluster_b[[1]],
    mean_d   = if (nrow(v) > 0) as.numeric(v$mean_marker_delta[[1]]) else NA_real_,
    max_d    = if (nrow(v) > 0) as.numeric(v$max_marker_delta[[1]])  else NA_real_,
    mk_max   = if (nrow(v) > 0) as.character(v$marker_with_max_delta[[1]]) else NA_character_,
    d_pass   = if (nrow(v) > 0) as.logical(v$delta_pass[[1]]) else NA,
    jsd      = if (nrow(v) > 0) as.numeric(v$jsd[[1]]) else NA_real_,
    j_pass   = if (nrow(v) > 0) as.logical(v$jsd_pass[[1]]) else NA,
    a_freq = a_freq, a_cv = a_cv, a_n = a_n,
    b_freq = b_freq, b_cv = b_cv, b_n = b_n,
    ratio  = ratio,
    stringsAsFactors = FALSE
  )
}))

pop_summary$status <- with(pop_summary, {
  red <- (
    (!is.na(cos)   & cos < 0.87)       |
    (!is.na(max_d) & max_d > 5.0)      |
    (!is.na(ratio) & ratio > 10)       |
    (!is.na(a_cv)  & a_cv > 3.0)       |
    (!is.na(b_cv)  & b_cv > 3.0)
  )
  green <- (
    (!is.na(cos)   & cos >= 0.95)                        &
    (is.na(j_pass) | (!is.na(j_pass) & j_pass == TRUE))  &
    (is.na(a_cv)   | a_cv < 1.5)                         &
    (is.na(b_cv)   | b_cv < 1.5)                         &
    (is.na(ratio)  | ratio < 5)
  )
  ifelse(red, "RED", ifelse(green, "GREEN", "YELLOW"))
})

STATUS_COL <- c(GREEN = "#27ae60", YELLOW = "#e67e22", RED = "#c0392b")
n_green <- sum(pop_summary$status == "GREEN")
n_yel   <- sum(pop_summary$status == "YELLOW")
n_red   <- sum(pop_summary$status == "RED")
message(sprintf("[Diag] Scorecard: %d GREEN | %d YELLOW | %d RED", n_green, n_yel, n_red))

# ── 3. Centroid extraction & per-marker delta matrix ──────────────────────────

.get_cent_row <- function(batch_name, cid) {
  m <- cent[[batch_name]]
  rid <- as.character(cid)
  if (rid %in% rownames(m)) m[rid, markers, drop = TRUE] else setNames(rep(NA_real_, length(markers)), markers)
}

cent_a <- do.call(rbind, lapply(pop_summary$clust_a, .get_cent_row, batch_name = batch_a))
cent_b <- do.call(rbind, lapply(pop_summary$clust_b, .get_cent_row, batch_name = batch_b))
rownames(cent_a) <- rownames(cent_b) <- pop_summary$pop_id

delta_mat <- abs(cent_a - cent_b)

# ── 4. Threshold saturation ───────────────────────────────────────────────────

pooled_cent <- rbind(cent_a, cent_b)

thresh_df <- do.call(rbind, lapply(markers, function(m) {
  vals <- pooled_cent[, m]
  vals <- vals[!is.na(vals)]
  if (length(vals) < 4)
    return(data.frame(marker = m, threshold = NA_real_, pct_above = NA_real_,
                      saturated = NA, stringsAsFactors = FALSE))
  km <- tryCatch(kmeans(vals, centers = 2, nstart = 20, iter.max = 100), error = function(e) NULL)
  if (is.null(km))
    return(data.frame(marker = m, threshold = NA_real_, pct_above = NA_real_,
                      saturated = NA, stringsAsFactors = FALSE))
  thr <- mean(sort(km$centers))
  pct <- mean(vals > thr)
  data.frame(marker = m, threshold = thr, pct_above = pct,
             saturated = pct > 0.75 | pct < 0.25, stringsAsFactors = FALSE)
}))

sat_markers <- thresh_df$marker[!is.na(thresh_df$saturated) & thresh_df$saturated]
message(sprintf("[Diag] Saturated markers (%d/%d): %s",
  length(sat_markers), length(markers),
  if (length(sat_markers) == 0) "none" else paste(sat_markers, collapse = ", ")))

# ── 5. Gate cascade from step 02 log ─────────────────────────────────────────

read_latest_step02_log <- function(logs_dir) {
  candidates <- c(
    list.files(file.path(logs_dir, "steps"), pattern = "^02_.*\\.json$", full.names = TRUE),
    unlist(lapply(list.dirs(logs_dir, recursive = FALSE), function(rd)
      list.files(file.path(rd, "steps"), pattern = "^02_.*\\.json$", full.names = TRUE)))
  )
  if (length(candidates) == 0L) return(NULL)
  f <- candidates[which.max(file.mtime(candidates))]
  tryCatch(jsonlite::fromJSON(f, simplifyDataFrame = FALSE), error = function(e) NULL)
}

gate_log  <- read_latest_step02_log(config$directories$logs %||% "results/logs")
gate_data <- NULL

if (!is.null(gate_log) && !is.null(gate_log$metrics$cd3_gate_per_sample)) {
  cd3_s <- gate_log$metrics$cd3_gate_per_sample
  cd8_s <- gate_log$metrics$cd8_gate_per_sample %||% list()
  min_cells <- as.integer(config$qc$min_cells_final %||% 5000L)

  gate_rows <- lapply(names(cd3_s), function(sid) {
    cd3 <- cd3_s[[sid]]
    cd8 <- cd8_s[[sid]] %||% list()
    sb  <- fm$batch[fm$sample_id == sid]
    data.frame(
      sample_id    = sid,
      batch        = if (length(sb) > 0) sb[1] else NA_character_,
      cd3_thr      = as.numeric(cd3$threshold  %||% NA_real_),
      cd3_pct      = as.numeric(cd3$pct_cd3    %||% NA_real_),
      cd3_n_after  = as.integer(cd3$n_after    %||% NA_integer_),
      cd8_thr      = as.numeric(cd8$threshold  %||% NA_real_),
      cd8_pct      = as.numeric(cd8$pct_cd8    %||% NA_real_),
      cd8_n_after  = as.integer(cd8$n_after    %||% NA_integer_),
      stringsAsFactors = FALSE
    )
  })
  gate_data <- do.call(rbind, gate_rows)
  gate_data$cd3_anomaly    <- !is.na(gate_data$cd3_thr) & gate_data$cd3_thr < 2.0
  gate_data$cd3_low_cells  <- !is.na(gate_data$cd3_n_after) & gate_data$cd3_n_after < min_cells
  gate_data$has_cd8        <- !is.na(gate_data$cd8_thr)
  message(sprintf("[Diag] Gate log loaded: %d samples, %d CD3 anomalies, %d below min_cells",
                  nrow(gate_data), sum(gate_data$cd3_anomaly), sum(gate_data$cd3_low_cells)))
} else {
  message("[Diag] No step 02 gate log found — skipping gate cascade page.")
}

# ── 6. Plot functions ─────────────────────────────────────────────────────────

# Page 1 — Population scorecard ------------------------------------------------
make_scorecard <- function(ps) {
  metric_defs <- list(
    list(name = "Cosine",      col = "cos",   fmt = function(x) sprintf("%.3f", x),
         ok = function(x) !is.na(x) & x >= 0.95),
    list(name = "MaxΔ",        col = "max_d", fmt = function(x) sprintf("%.1f", x),
         ok = function(x) !is.na(x) & x < 2.0),
    list(name = "Δpass",       col = "d_pass",fmt = function(x) ifelse(isTRUE(x), "✓", "✗"),
         ok = function(x) !is.na(x) & isTRUE(x)),
    list(name = "JSDpass",     col = "j_pass",fmt = function(x) ifelse(isTRUE(x), "✓", "✗"),
         ok = function(x) !is.na(x) & isTRUE(x)),
    list(name = sprintf("CV\n%s", batch_a), col = "a_cv", fmt = function(x) sprintf("%.2f", x),
         ok = function(x) !is.na(x) & x < 1.5),
    list(name = sprintf("CV\n%s", batch_b), col = "b_cv", fmt = function(x) sprintf("%.2f", x),
         ok = function(x) !is.na(x) & x < 1.5),
    list(name = "Freq\nratio", col = "ratio", fmt = function(x) if (is.na(x)) "—" else sprintf("%.1f×", x),
         ok = function(x) !is.na(x) & x < 5)
  )
  n_m <- length(metric_defs)

  rows <- do.call(rbind, lapply(seq_len(nrow(ps)), function(i) {
    p <- ps[i, ]
    do.call(rbind, lapply(seq_along(metric_defs), function(j) {
      md  <- metric_defs[[j]]
      val <- p[[md$col]]
      data.frame(row = i, col = j,
                 val_txt = md$fmt(val),
                 ok      = md$ok(val),
                 stringsAsFactors = FALSE)
    }))
  }))
  rows$fill <- ifelse(rows$ok, "#d5f5e3", "#fadbd8")
  rows$tcol <- ifelse(rows$ok, "#1e8449", "#c0392b")

  # Label column (col = 0)
  label_df <- data.frame(
    row = seq_len(nrow(ps)), col = 0,
    val_txt = substr(ps$label, 1, 42),
    ok = TRUE, fill = "#ecf0f1", tcol = "#2c3e50", stringsAsFactors = FALSE)

  # Status column (col = n_m + 1)
  stat_df <- data.frame(
    row = seq_len(nrow(ps)), col = n_m + 1,
    val_txt = ps$status,
    ok = ps$status == "GREEN",
    fill = STATUS_COL[ps$status],
    tcol = "white", stringsAsFactors = FALSE)

  all_df <- rbind(label_df, rows, stat_df)

  header <- data.frame(
    col   = c(0, seq_len(n_m), n_m + 1),
    label = c("Population", sapply(metric_defs, `[[`, "name"), "Status"),
    stringsAsFactors = FALSE)

  ggplot(all_df, aes(x = col, y = -row)) +
    geom_tile(aes(fill = fill), color = "white", width = 0.97, height = 0.9) +
    scale_fill_identity() +
    geom_text(aes(label = val_txt, color = tcol), size = 2.5, lineheight = 0.85) +
    scale_color_identity() +
    geom_text(data = header, aes(x = col, y = 0.45, label = label),
              fontface = "bold", size = 2.8, lineheight = 0.82, inherit.aes = FALSE) +
    scale_x_continuous(expand = expansion(add = 0.55)) +
    scale_y_continuous(expand = expansion(add = 0.6)) +
    labs(title = "Population Integration Scorecard",
         subtitle = sprintf("%d GREEN  ·  %d YELLOW  ·  %d RED", n_green, n_yel, n_red)) +
    theme_void() +
    theme(plot.title    = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 3)),
          plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#7f8c8d"),
          plot.margin   = margin(10, 10, 10, 10))
}

# Page 2 — Marker delta heatmap ------------------------------------------------
make_delta_heatmap <- function(delta_mat, ps) {
  df <- do.call(rbind, lapply(rownames(delta_mat), function(p) {
    data.frame(pop_id = p, marker = colnames(delta_mat),
               delta = as.numeric(delta_mat[p, ]), stringsAsFactors = FALSE)
  }))
  df$pop_label <- factor(
    paste0(df$pop_id, " — ", substr(ps$label[match(df$pop_id, ps$pop_id)], 1, 38)),
    levels = rev(paste0(ps$pop_id, " — ", substr(ps$label, 1, 38))))
  df$border <- ifelse(!is.na(df$delta) & df$delta > 2.0, "fail", "ok")

  ggplot(df, aes(x = marker, y = pop_label)) +
    geom_tile(aes(fill = delta, color = border), linewidth = 0.6) +
    scale_fill_gradientn(
      colours  = c("#f7fbff", "#6baed6", "#08519c", "#e6550d", "#a50f15"),
      values   = scales::rescale(c(0, 1.5, 3, 5, 7)),
      na.value = "grey85",
      name     = "|Δ arcsinh|") +
    scale_color_manual(values = c(fail = "#c0392b", ok = "white"), guide = "none") +
    geom_text(aes(label = ifelse(is.na(delta), "", sprintf("%.1f", delta))),
              size = 2.9, color = "white", fontface = "bold") +
    labs(title    = "Per-Marker Absolute Delta Between Batches",
         subtitle = "Red border: delta > 2.0 arcsinh (fails delta_pass threshold)",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y   = element_text(size = 8),
          plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "#7f8c8d"),
          legend.position = "right")
}

# Page 3 — Threshold saturation ------------------------------------------------
make_threshold_plots <- function(pooled_cent, thresh_df) {
  plots <- lapply(markers, function(m) {
    vals <- pooled_cent[, m]
    vals <- vals[!is.na(vals)]
    td   <- thresh_df[thresh_df$marker == m, ]
    thr  <- if (nrow(td) > 0) td$threshold[1] else NA_real_
    pct  <- if (nrow(td) > 0) td$pct_above[1] else NA_real_
    sat  <- if (nrow(td) > 0) isTRUE(td$saturated[1]) else FALSE

    df_hist <- data.frame(x = vals)
    title_col <- if (sat) "#c0392b" else "#2c3e50"
    subtitle  <- if (is.na(pct)) "" else
      sprintf("%.0f%% above threshold%s", pct * 100, if (sat) "  [SATURATED]" else "")

    p <- ggplot(df_hist, aes(x = x)) +
      geom_histogram(bins = 12, fill = "#4a90d9", color = "white", alpha = 0.8) +
      labs(title = m, subtitle = subtitle, x = "arcsinh expression", y = "count") +
      theme_minimal(base_size = 9) +
      theme(plot.title    = element_text(face = "bold", color = title_col, size = 10),
            plot.subtitle = element_text(size = 8, color = if (sat) "#c0392b" else "#7f8c8d"),
            axis.title    = element_text(size = 8))
    if (!is.na(thr))
      p <- p + geom_vline(xintercept = thr, color = "#e74c3c", linetype = "dashed", linewidth = 1)
    p
  })
  wrap_plots(plots, ncol = 3) +
    plot_annotation(
      title    = "Auto-Phenotyping Threshold Saturation",
      subtitle = "Red dashed line = k-means threshold on pooled matched centroids. Saturated (>75% or <25% above) = all populations get same label.",
      theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                       plot.subtitle = element_text(size = 9, color = "#7f8c8d")))
}

# Page 4 — Frequency correlation scatter ----------------------------------------
make_freq_scatter <- function(ps) {
  ps$a_pct <- ps$a_freq * 100
  ps$b_pct <- ps$b_freq * 100
  lim_max  <- max(c(ps$a_pct, ps$b_pct), na.rm = TRUE) * 1.08

  r2_val <- tryCatch({
    ok <- !is.na(ps$a_pct) & !is.na(ps$b_pct)
    if (sum(ok) > 2) round(cor(ps$a_pct[ok], ps$b_pct[ok])^2, 3) else NA_real_
  }, error = function(e) NA_real_)

  ggplot(ps, aes(x = a_pct, y = b_pct, color = status)) +
    geom_abline(slope = 1, intercept = 0, color = "grey70", linetype = "dashed") +
    geom_hline(yintercept = 5, color = "#e8e8e8") +
    geom_vline(xintercept = 5, color = "#e8e8e8") +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = STATUS_COL, name = "Status") +
    ggrepel::geom_label_repel(
      aes(label = pop_id), size = 2.8, max.overlaps = 20,
      box.padding = 0.3, label.padding = 0.15, show.legend = FALSE) +
    scale_x_continuous(limits = c(0, lim_max), labels = function(x) paste0(x, "%")) +
    scale_y_continuous(limits = c(0, lim_max), labels = function(x) paste0(x, "%")) +
    labs(title    = "Cross-Batch Frequency Concordance",
         subtitle = if (!is.na(r2_val)) sprintf("R² = %.3f  (diagonal = perfect concordance)", r2_val)
                    else "Diagonal = perfect concordance",
         x = sprintf("Mean frequency — %s", batch_a),
         y = sprintf("Mean frequency — %s", batch_b)) +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "#7f8c8d"),
          legend.position = "right")
}

# Page 5 — CV distribution -----------------------------------------------------
make_cv_plot <- function(ps) {
  df_cv <- rbind(
    data.frame(pop_id = ps$pop_id, cv = ps$a_cv, batch = batch_a,
               status = ps$status, stringsAsFactors = FALSE),
    data.frame(pop_id = ps$pop_id, cv = ps$b_cv, batch = batch_b,
               status = ps$status, stringsAsFactors = FALSE)
  )
  df_cv <- df_cv[!is.na(df_cv$cv), ]
  df_cv$pop_id <- factor(df_cv$pop_id, levels = rev(ps$pop_id[order(rowMeans(cbind(ps$a_cv, ps$b_cv), na.rm = TRUE))]))

  ggplot(df_cv, aes(x = cv, y = pop_id, color = batch, shape = batch)) +
    geom_vline(xintercept = 1.5, color = "#e74c3c", linetype = "dashed", linewidth = 0.8) +
    geom_point(size = 3, alpha = 0.85, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c("#2980b9", "#e67e22"), name = "Batch") +
    scale_shape_manual(values = c(16, 17), name = "Batch") +
    annotate("text", x = 1.52, y = nlevels(df_cv$pop_id) + 0.3,
             label = "CV = 1.5", hjust = 0, size = 3, color = "#c0392b") +
    labs(title    = "Intra-Batch Variability (CV = SD / Mean frequency)",
         subtitle = "CV > 1.5: baseline too noisy for reliable Z-scoring",
         x = "Coefficient of Variation (SD / mean freq)",
         y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "#7f8c8d"),
          legend.position = "bottom")
}

# Page 7 — Gate cascade summary (CD3 / CD8 thresholds per sample) --------------
make_gate_plot <- function(gd) {
  batch_cols <- c("#2980b9", "#e67e22")
  batches    <- sort(unique(gd$batch[!is.na(gd$batch)]))
  names(batch_cols) <- batches

  # CD3 threshold distribution
  p_cd3_thr <- ggplot(gd[!is.na(gd$cd3_thr), ],
                      aes(x = batch, y = cd3_thr, color = batch)) +
    geom_jitter(aes(shape = cd3_anomaly), width = 0.15, size = 2.5, alpha = 0.85) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 8),
                       labels = c("Normal", "ANOMALY (thr < 2.0)"), name = NULL) +
    scale_color_manual(values = batch_cols, guide = "none") +
    geom_hline(yintercept = 2.0, linetype = "dashed", color = "#c0392b", linewidth = 0.7) +
    annotate("text", x = 0.55, y = 2.1, label = "anomaly < 2.0", hjust = 0,
             size = 2.8, color = "#c0392b") +
    geom_hline(yintercept = as.numeric(config$cd3_gate$fallback_threshold %||% 4.0),
               linetype = "dotted", color = "#8e44ad", linewidth = 0.7) +
    annotate("text", x = 0.55,
             y = as.numeric(config$cd3_gate$fallback_threshold %||% 4.0) + 0.1,
             label = "fallback", hjust = 0, size = 2.8, color = "#8e44ad") +
    ggrepel::geom_text_repel(
      data = gd[!is.na(gd$cd3_thr) & gd$cd3_anomaly, ],
      aes(label = substr(sample_id, 1, 20)), size = 2.3, color = "#c0392b",
      max.overlaps = 10, show.legend = FALSE) +
    labs(title = "CD3 k-means threshold per sample", x = NULL,
         y = "arcsinh threshold") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          legend.position = "bottom", legend.text = element_text(size = 8))

  # CD3 %positive distribution
  p_cd3_pct <- ggplot(gd[!is.na(gd$cd3_pct), ],
                      aes(x = batch, y = cd3_pct, color = batch)) +
    geom_jitter(aes(shape = cd3_low_cells), width = 0.15, size = 2.5, alpha = 0.85) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 8),
                       labels = c("OK", "< min_cells"), name = NULL) +
    scale_color_manual(values = batch_cols, guide = "none") +
    geom_hline(yintercept = as.numeric(config$cd3_gate$min_cd3_pct %||% 10),
               linetype = "dashed", color = "#c0392b", linewidth = 0.7) +
    labs(title = "CD3+ % per sample (post-gate)", x = NULL,
         y = "% CD3+ cells") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          legend.position = "bottom", legend.text = element_text(size = 8))

  has_cd8 <- any(!is.na(gd$cd8_thr))

  if (has_cd8) {
    p_cd8_thr <- ggplot(gd[!is.na(gd$cd8_thr), ],
                        aes(x = batch, y = cd8_thr, color = batch)) +
      geom_jitter(width = 0.15, size = 2.5, alpha = 0.85) +
      scale_color_manual(values = batch_cols, guide = "none") +
      geom_hline(yintercept = as.numeric(config$cd8_gate$fallback_threshold %||% 3.5),
                 linetype = "dotted", color = "#8e44ad", linewidth = 0.7) +
      annotate("text", x = 0.55,
               y = as.numeric(config$cd8_gate$fallback_threshold %||% 3.5) + 0.1,
               label = "fallback", hjust = 0, size = 2.8, color = "#8e44ad") +
      labs(title = "CD8 k-means threshold per sample", x = NULL,
           y = "arcsinh threshold") +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold", size = 11))

    p_cd8_pct <- ggplot(gd[!is.na(gd$cd8_pct), ],
                        aes(x = batch, y = cd8_pct, color = batch)) +
      geom_jitter(width = 0.15, size = 2.5, alpha = 0.85) +
      scale_color_manual(values = batch_cols, guide = "none") +
      geom_hline(yintercept = as.numeric(config$cd8_gate$min_cd8_pct %||% 10),
                 linetype = "dashed", color = "#c0392b", linewidth = 0.7) +
      labs(title = "CD8+ % per sample (post-CD3 gate)", x = NULL,
           y = "% CD8+ of CD3+ cells") +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold", size = 11))

    (p_cd3_thr | p_cd3_pct) / (p_cd8_thr | p_cd8_pct) +
      plot_annotation(
        title    = "Gate Cascade — CD3 / CD8 k-means Thresholds",
        subtitle = "Dashed red = warning threshold  ·  Dotted purple = fallback threshold  ·  Star = anomaly",
        theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                         plot.subtitle = element_text(size = 9, color = "#7f8c8d")))
  } else {
    (p_cd3_thr | p_cd3_pct) +
      plot_annotation(
        title    = "Gate Cascade — CD3 k-means Thresholds per Sample",
        subtitle = "Dashed red = min_cd3_pct warning  ·  Dotted purple = fallback threshold  ·  Star = anomaly",
        theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                         plot.subtitle = element_text(size = 9, color = "#7f8c8d")))
  }
}

# Page 6 — Per-sample frequency distributions ----------------------------------
make_sample_dist <- function(fm, ps, n_show = 9L) {
  # Select most informative populations: greens first, then lowest-CV yellows
  ps_ord <- ps[order(
    match(ps$status, c("GREEN", "YELLOW", "RED")),
    rowMeans(cbind(ps$a_cv, ps$b_cv), na.rm = TRUE)
  ), ]
  sel_pops <- head(ps_ord$pop_id, n_show)

  pop_cols <- intersect(sel_pops, colnames(fm))
  if (length(pop_cols) == 0L) {
    return(ggplot() + labs(title = "No frequency data available") + theme_void())
  }

  fm_long <- do.call(rbind, lapply(pop_cols, function(p) {
    data.frame(sample_id  = fm$sample_id,
               batch      = fm$batch,
               population = p,
               freq_pct   = fm[[p]] * 100,
               stringsAsFactors = FALSE)
  }))
  fm_long$pop_label <- factor(
    fm_long$population,
    levels = sel_pops,
    labels = substr(ps$label[match(sel_pops, ps$pop_id)], 1, 30))
  fm_long$status <- ps$status[match(fm_long$population, ps$pop_id)]

  ggplot(fm_long, aes(x = batch, y = freq_pct, color = batch)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.6) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
    scale_color_manual(values = c("#2980b9", "#e67e22"), guide = "none") +
    facet_wrap(~pop_label, scales = "free_y", ncol = 3) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(title    = "Per-Sample Frequency Distributions (top populations)",
         subtitle = "Ordered by integration quality (GREEN first); y-axes free",
         x = NULL, y = "Frequency (% of total cells)") +
    theme_minimal(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "#7f8c8d"),
          axis.text.x   = element_text(angle = 30, hjust = 1, size = 8),
          strip.text    = element_text(size = 8, face = "bold"),
          panel.spacing = unit(0.8, "lines"))
}

# ── 6. Build and save PDF ─────────────────────────────────────────────────────

message("[Diag] Rendering plots...")

p1 <- make_scorecard(pop_summary)
p2 <- make_delta_heatmap(delta_mat, pop_summary)
p3 <- make_threshold_plots(pooled_cent, thresh_df)
p4 <- make_freq_scatter(pop_summary)
p5 <- make_cv_plot(pop_summary)
p6 <- make_sample_dist(fm, pop_summary)

grDevices::pdf(out_pdf, width = 14, height = 10)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
if (!is.null(gate_data) && nrow(gate_data) > 0L) {
  print(make_gate_plot(gate_data))
}
invisible(grDevices::dev.off())

message(sprintf("[Diag] PDF written: %s", out_pdf))

# ── 7. JSON summary ───────────────────────────────────────────────────────────

clinically_usable  <- pop_summary$pop_id[pop_summary$status == "GREEN"]
borderline         <- pop_summary$pop_id[pop_summary$status == "YELLOW"]
unreliable         <- pop_summary$pop_id[pop_summary$status == "RED"]

freq_paradox <- pop_summary[!is.na(pop_summary$ratio) & pop_summary$ratio > 5,
                             c("pop_id", "label", "ratio", "a_freq", "b_freq")]

# Rough n estimate: CV ~ 1/(sqrt(n)*something); target CV < 1.0 needs ~4x current n
high_cv_pops <- pop_summary[!is.na(pop_summary$a_cv) & pop_summary$a_cv > 1.5, ]
n_needed_estimate <- if (nrow(high_cv_pops) > 0) {
  mean_cv  <- mean(c(high_cv_pops$a_cv, high_cv_pops$b_cv), na.rm = TRUE)
  target   <- 1.0
  # CV ~ sigma/mu; SD scales as 1/sqrt(n); so n_needed ~ n_current * (cv/target)^2
  current_n <- mean(c(pop_summary$a_n, pop_summary$b_n), na.rm = TRUE)
  round(current_n * (mean_cv / target)^2)
} else NULL

recommendations <- c(
  if (length(clinically_usable) > 0)
    sprintf("Use for clinical Z-scoring: %s", paste(clinically_usable, collapse = ", ")),
  if (length(sat_markers) > 0)
    sprintf("Ignore functional tags from saturated markers: %s — threshold falls outside bimodal range, all populations receive same positive/negative label",
            paste(sat_markers, collapse = ", ")),
  if (nrow(freq_paradox) > 0)
    sprintf("Frequency paradox populations (ratio > 5x): %s — cosine similarity is high but panels resolve this compartment differently; exclude from cross-batch comparison",
            paste(freq_paradox$pop_id, collapse = ", ")),
  if (length(unreliable) > 0)
    sprintf("Unreliable populations (CV > 3 or max_delta > 5): %s — too noisy or too discordant for a defensible baseline",
            paste(unreliable, collapse = ", ")),
  if (!is.null(n_needed_estimate))
    sprintf("Estimated additional samples needed to bring CV below 1.0 in high-CV populations: ~%d per batch (current mean CV in flagged pops: %.2f)",
            max(0L, n_needed_estimate - round(mean(c(pop_summary$a_n[1], pop_summary$b_n[1]), na.rm = TRUE))),
            mean(c(high_cv_pops$a_cv, high_cv_pops$b_cv), na.rm = TRUE))
)

summary_json <- list(
  generated_at         = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  n_matched_populations = length(matched_pops),
  batches              = list(a = batch_a, b = batch_b),
  markers              = markers,
  scorecard = list(
    green  = clinically_usable,
    yellow = borderline,
    red    = unreliable
  ),
  threshold_saturation = list(
    saturated_markers     = sat_markers,
    per_marker = lapply(seq_len(nrow(thresh_df)), function(i) {
      as.list(thresh_df[i, ])
    })
  ),
  frequency_paradox = lapply(seq_len(nrow(freq_paradox)), function(i) {
    list(pop_id = freq_paradox$pop_id[i],
         label  = freq_paradox$label[i],
         ratio  = round(freq_paradox$ratio[i], 2),
         a_freq = round(freq_paradox$a_freq[i], 4),
         b_freq = round(freq_paradox$b_freq[i], 4))
  }),
  n_samples_needed_estimate = n_needed_estimate,
  gate_cascade_summary = if (!is.null(gate_data) && nrow(gate_data) > 0L) {
    list(
      n_samples       = nrow(gate_data),
      n_cd3_anomalies = sum(gate_data$cd3_anomaly),
      n_below_min_cells = sum(gate_data$cd3_low_cells),
      anomalous_samples = gate_data$sample_id[gate_data$cd3_anomaly],
      low_cell_samples  = gate_data$sample_id[gate_data$cd3_low_cells],
      per_batch = lapply(split(gate_data, gate_data$batch), function(b) list(
        n_samples       = nrow(b),
        mean_cd3_thr    = round(mean(b$cd3_thr, na.rm = TRUE), 3),
        sd_cd3_thr      = round(sd(b$cd3_thr,   na.rm = TRUE), 3),
        mean_cd3_pct    = round(mean(b$cd3_pct,  na.rm = TRUE), 1),
        n_cd3_anomalies = sum(b$cd3_anomaly),
        mean_cd8_thr    = if (any(!is.na(b$cd8_thr))) round(mean(b$cd8_thr, na.rm = TRUE), 3) else NULL,
        mean_cd8_pct    = if (any(!is.na(b$cd8_pct))) round(mean(b$cd8_pct, na.rm = TRUE), 1) else NULL
      ))
    )
  } else NULL,
  recommendations      = recommendations
)

jsonlite::write_json(summary_json, out_json, pretty = TRUE, auto_unbox = TRUE)
message(sprintf("[Diag] JSON written: %s", out_json))

message("\n=== RECOMMENDATIONS ===")
for (r in recommendations) message("• ", r)
message(sprintf("\nFull report: %s", out_pdf))
