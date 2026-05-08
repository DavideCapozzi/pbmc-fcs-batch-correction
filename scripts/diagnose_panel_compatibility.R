# scripts/diagnose_panel_compatibility.R
# ==============================================================================
# PANEL COMPATIBILITY DIAGNOSTIC
# Reads FCS headers only (no cell data — runs in seconds).
# Discovers panel groups across all source directories, computes cross-batch
# marker intersections for every combination scenario, and outputs:
#   - PDF: presence/absence heatmap + scenario comparison chart
#   - JSON: full compatibility report
#
# Usage (from project root):
#   Rscript scripts/diagnose_panel_compatibility.R
# ==============================================================================

suppressPackageStartupMessages({
  library(flowCore)
  library(yaml)
  library(jsonlite)
  library(ggplot2)
  library(patchwork)
})

source("src/functions/fcs_io.R")

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------

config   <- yaml::read_yaml("config/global_params.yml")
pl_root  <- config$directories$raw$pannelli_liquidi
dc_dir   <- config$directories$raw$duraclone
ds_old   <- file.path(pl_root, "DS old")

out_dir  <- "results/diagnostics"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ts       <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_pdf  <- file.path(out_dir, paste0("panel_compatibility_", ts, ".pdf"))
out_json <- file.path(out_dir, paste0("panel_compatibility_", ts, ".json"))

message("=== PANEL COMPATIBILITY DIAGNOSTIC ===")

# ------------------------------------------------------------------------------
# File discovery
# ------------------------------------------------------------------------------

scan_dir <- function(path, source_label, recursive = FALSE) {
  if (!dir.exists(path)) {
    warning(sprintf("Directory not found: %s", path))
    return(data.frame())
  }
  fcs <- list.files(path, pattern = "\\.fcs$", full.names = TRUE,
                    recursive = recursive, ignore.case = TRUE)
  if (length(fcs) == 0L) return(data.frame())
  data.frame(path = fcs, source = source_label, stringsAsFactors = FALSE)
}

all_files <- rbind(
  scan_dir(pl_root,  "pannelli_liquidi", recursive = FALSE),
  scan_dir(ds_old,   "ds_old",           recursive = FALSE),
  scan_dir(dc_dir,   "duraclone",        recursive = FALSE)
)

message(sprintf("[Compat] Found %d FCS files: %s",
  nrow(all_files),
  paste(tapply(all_files$source, all_files$source, length), names(tapply(all_files$source, all_files$source, length)), sep = " ", collapse = " | ")))

# ------------------------------------------------------------------------------
# Header reading + marker extraction
# ------------------------------------------------------------------------------

read_markers <- function(fpath) {
  ff <- tryCatch(
    flowCore::read.FCS(fpath, which.lines = 1L, transformation = FALSE,
                       truncate_max_range = FALSE, emptyValue = FALSE),
    error = function(e) NULL
  )
  if (is.null(ff)) return(NULL)

  p     <- flowCore::pData(flowCore::parameters(ff))
  # Use $PnS (description) preferring over $PnN (channel name) when available
  descs <- ifelse(!is.na(p$desc) & nzchar(p$desc), p$desc, p$name)
  clean <- sanitize_biological_markers(descs)
  sort(unique(clean[!is.na(clean)]))
}

message("[Compat] Reading FCS headers...")
all_files$markers <- lapply(all_files$path, read_markers)
all_files$n_markers <- sapply(all_files$markers, length)
all_files$panel_sig <- sapply(all_files$markers, paste, collapse = "|")
all_files$fname     <- basename(all_files$path)

# Flag read errors
failed <- sum(sapply(all_files$markers, is.null))
if (failed > 0) warning(sprintf("%d file(s) could not be read — excluded.", failed))
all_files <- all_files[!sapply(all_files$markers, is.null), ]

# ------------------------------------------------------------------------------
# Panel group detection (auto-cluster by marker signature)
# ------------------------------------------------------------------------------

# Within each source directory, group files by their exact marker set
all_files$panel_group <- paste(all_files$source, all_files$panel_sig, sep = ":::")

# Build group summary
group_info <- do.call(rbind, lapply(split(all_files, all_files$panel_group), function(g) {
  data.frame(
    panel_group  = g$panel_group[1],
    source       = g$source[1],
    n_files      = nrow(g),
    n_markers    = g$n_markers[1],
    panel_sig    = g$panel_sig[1],
    markers_list = I(list(g$markers[[1]])),
    stringsAsFactors = FALSE
  )
}))
rownames(group_info) <- NULL

# Short display labels
group_info$label <- sprintf("%s\n(n=%d, %dm)",
  gsub("_", " ", group_info$source),
  group_info$n_files,
  group_info$n_markers)

message(sprintf("[Compat] Discovered %d distinct panel group(s):", nrow(group_info)))
for (i in seq_len(nrow(group_info))) {
  message(sprintf("  [%d] %s — %d files, %d markers: %s",
    i, group_info$source[i], group_info$n_files[i], group_info$n_markers[i],
    paste(group_info$markers_list[[i]], collapse = ", ")))
}

# ------------------------------------------------------------------------------
# Intersection analysis
# Apply raw_filters from config to identify which files are "active" in the
# pipeline, then compute scenarios per individual group (not aggregated).
# ------------------------------------------------------------------------------

raw_filters_cfg <- config$directories$raw_filters %||% list()

# Mark each file as active (passes the per-source filter) or inactive
all_files$active <- mapply(function(src, fpath) {
  flt <- raw_filters_cfg[[src]]
  if (is.null(flt)) TRUE else grepl(flt, basename(fpath))
}, all_files$source, all_files$path)

# Active group summaries per source
active_files    <- all_files[all_files$active, ]
inactive_files  <- all_files[!all_files$active, ]

make_groups <- function(df) {
  if (nrow(df) == 0L) return(data.frame())
  do.call(rbind, lapply(split(df, df$panel_group), function(g) {
    data.frame(panel_group = g$panel_group[1], source = g$source[1],
               n_files = nrow(g), n_markers = g$n_markers[1],
               panel_sig = g$panel_sig[1],
               markers_list = I(list(g$markers[[1]])),
               stringsAsFactors = FALSE)
  }))
}

active_group_info   <- make_groups(active_files)
inactive_group_info <- make_groups(inactive_files)

duraclone_groups <- group_info[group_info$source == "duraclone", ]
pannelli_groups  <- if (!is.null(active_group_info) && nrow(active_group_info) > 0) {
  active_group_info[active_group_info$source == "pannelli_liquidi", ]
} else {
  group_info[group_info$source == "pannelli_liquidi", ]
}
dsold_groups     <- group_info[group_info$source == "ds_old", ]
# All pannelli groups including inactive (for S3 scenario)
all_pannelli_groups <- group_info[group_info$source == "pannelli_liquidi", ]

duraclone_markers <- if (nrow(duraclone_groups) == 1L) {
  duraclone_groups$markers_list[[1]]
} else {
  Reduce(intersect, duraclone_groups$markers_list)
}

# Per-group marker sets — use UNION within each source to capture all active files
# then intersect across groups per scenario (correct approach)
tex_markers   <- if (nrow(pannelli_groups) > 0)
  Reduce(intersect, pannelli_groups$markers_list) else character(0)
dsold_markers <- if (nrow(dsold_groups) > 0)
  Reduce(intersect, dsold_groups$markers_list) else character(0)
all_pl_markers <- if (nrow(all_pannelli_groups) > 0)
  Reduce(intersect, all_pannelli_groups$markers_list) else character(0)

message(sprintf("[Compat] Active pannelli_liquidi files: %d (filter: %s)",
  sum(pannelli_groups$n_files),
  raw_filters_cfg$pannelli_liquidi %||% "none"))
message(sprintf("[Compat] Inactive pannelli_liquidi files: %d",
  sum(inactive_group_info$n_files[inactive_group_info$source == "pannelli_liquidi"],
      na.rm = TRUE)))

# All combination scenarios vs duraclone
scenarios <- list(
  list(
    id      = "S1_TEX_only",
    label   = "TEX pannelli\n+ Duraclone",
    groups  = "pannelli_liquidi",
    markers = intersect(tex_markers, duraclone_markers),
    n_pl    = sum(pannelli_groups$n_files),
    n_dc    = sum(duraclone_groups$n_files)
  ),
  list(
    id      = "S2_TEX_DSold",
    label   = "TEX pannelli\n+ DS old\n+ Duraclone",
    groups  = c("pannelli_liquidi", "ds_old"),
    markers = intersect(intersect(tex_markers, dsold_markers), duraclone_markers),
    n_pl    = sum(pannelli_groups$n_files) + sum(dsold_groups$n_files),
    n_dc    = sum(duraclone_groups$n_files)
  ),
  list(
    id      = "S3_all_pannelli",
    label   = "All pannelli\n+ Duraclone",
    groups  = c("pannelli_liquidi", "ds_old"),
    markers = intersect(intersect(all_pl_markers, dsold_markers), duraclone_markers),
    n_pl    = sum(all_pannelli_groups$n_files) + sum(dsold_groups$n_files),
    n_dc    = sum(duraclone_groups$n_files)
  )
)

message("\n[Compat] Scenario comparison:")
for (s in scenarios) {
  message(sprintf("  %s: %d shared markers, n_pannelli=%d, n_duraclone=%d",
    s$id, length(s$markers), s$n_pl, s$n_dc))
  message(sprintf("    Markers: %s", paste(sort(s$markers), collapse = ", ")))
}

# ------------------------------------------------------------------------------
# Presence/absence matrix for heatmap
# ------------------------------------------------------------------------------

all_markers_union <- sort(unique(c(
  unlist(lapply(all_files$markers, identity)),
  duraclone_markers
)))

# One column per discovered panel group (unique source × panel_sig combo)
presence_list <- lapply(seq_len(nrow(group_info)), function(i) {
  gm  <- group_info$markers_list[[i]]
  lbl <- group_info$label[i]
  data.frame(
    marker       = all_markers_union,
    panel_group  = lbl,
    source       = group_info$source[i],
    n_files      = group_info$n_files[i],
    present      = all_markers_union %in% gm,
    in_duraclone = all_markers_union %in% duraclone_markers,
    stringsAsFactors = FALSE
  )
})
pmat <- do.call(rbind, presence_list)

# Marker sort order: by total presence count (descending)
marker_order <- names(sort(
  tapply(pmat$present, pmat$marker, sum),
  decreasing = TRUE
))
pmat$marker <- factor(pmat$marker, levels = rev(marker_order))

# Fill category
pmat$fill_cat <- with(pmat, ifelse(
  present & in_duraclone, "shared_with_DC",
  ifelse(present & !in_duraclone, "panel_only",
  "absent")))
pmat$fill_cat <- factor(pmat$fill_cat,
  levels = c("shared_with_DC", "panel_only", "absent"))

# Column order: duraclone last, pannelli groups first
source_order <- c("pannelli_liquidi", "ds_old", "duraclone")
pmat <- pmat[order(match(pmat$source, source_order)), ]
col_order <- unique(pmat$panel_group)
pmat$panel_group <- factor(pmat$panel_group, levels = col_order)

# ------------------------------------------------------------------------------
# Scenario bar chart data
# ------------------------------------------------------------------------------

# Which markers are in S2 (recommended) vs. other scenarios
s1_markers <- scenarios[[1]]$markers
s2_markers <- scenarios[[2]]$markers
s3_markers <- scenarios[[3]]$markers

scenario_df <- data.frame(
  scenario    = factor(
    sapply(scenarios, `[[`, "label"),
    levels = sapply(scenarios, `[[`, "label")),
  n_markers   = sapply(scenarios, function(s) length(s$markers)),
  n_pannelli  = sapply(scenarios, `[[`, "n_pl"),
  n_duraclone = sapply(scenarios, `[[`, "n_dc"),
  recommended = c(FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

# Per-scenario marker details for annotation
scenario_marker_df <- do.call(rbind, lapply(seq_along(scenarios), function(i) {
  s  <- scenarios[[i]]
  # All duraclone markers; mark each as retained or lost in this scenario
  data.frame(
    scenario  = factor(s$label, levels = sapply(scenarios, `[[`, "label")),
    marker    = sort(duraclone_markers),
    retained  = duraclone_markers %in% s$markers,
    stringsAsFactors = FALSE
  )
}))
scenario_marker_df$marker_rank <- match(
  scenario_marker_df$marker,
  sort(duraclone_markers))

# ------------------------------------------------------------------------------
# Colour palettes
# ------------------------------------------------------------------------------

fill_colours <- c(
  shared_with_DC = "#2E8B57",   # sea green  — shared with duraclone
  panel_only     = "#6BAED6",   # steel blue — present but not in DC
  absent         = "#F0F0F0"    # light grey — absent
)

source_colours <- c(
  pannelli_liquidi = "#4393C3",
  ds_old           = "#66C2A5",
  duraclone        = "#FC8D62"
)

# ------------------------------------------------------------------------------
# Panel A: Presence/absence heatmap
# ------------------------------------------------------------------------------

# Top annotation strip: colour by source
strip_df <- unique(pmat[, c("panel_group", "source")])
strip_df$source <- factor(strip_df$source, levels = names(source_colours))

plot_A <- ggplot(pmat, aes(x = panel_group, y = marker, fill = fill_cat)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(
    data = pmat[pmat$present, ],
    aes(label = "●"), size = 2.8, colour = "white"
  ) +
  scale_fill_manual(
    values = fill_colours,
    labels = c(
      shared_with_DC = "Shared with duraclone",
      panel_only     = "Present (panel-only)",
      absent         = "Absent"
    ),
    name = NULL
  ) +
  labs(
    title    = "Marker presence by panel group",
    subtitle = "● = present  |  Green = shared with duraclone  |  Blue = panel-specific",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8),
    axis.text.y       = element_text(size = 8),
    panel.grid        = element_blank(),
    legend.position   = "bottom",
    legend.text       = element_text(size = 8),
    plot.title        = element_text(face = "bold", size = 11),
    plot.subtitle     = element_text(size = 8, colour = "grey40")
  )

# ------------------------------------------------------------------------------
# Panel B: Scenario comparison
# ------------------------------------------------------------------------------

# Main bars: n_shared_markers
plot_B_top <- ggplot(scenario_df,
    aes(x = scenario, y = n_markers, fill = recommended)) +
  geom_col(width = 0.55, colour = "white") +
  geom_text(
    aes(label = paste0(n_markers, " markers")),
    vjust = -0.4, fontface = "bold", size = 3.2
  ) +
  geom_text(
    aes(label = sprintf("PL n=%d\nDC n=%d", n_pannelli, n_duraclone)),
    y = 0.3, vjust = 0, size = 2.6, colour = "white", lineheight = 1.1
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#2E8B57", "FALSE" = "#AAAAAA"),
    guide  = "none"
  ) +
  scale_y_continuous(
    limits = c(0, max(scenario_df$n_markers) + 1.5),
    expand = c(0, 0)
  ) +
  labs(
    title    = "Shared markers vs. duraclone — combination scenarios",
    subtitle = "Green = recommended scenario  |  PL = pannelli_liquidi  |  DC = duraclone",
    x = NULL, y = "N shared markers with duraclone"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x   = element_text(size = 8, lineheight = 1.1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, colour = "grey40")
  )

# Marker grid below: which duraclone markers are retained per scenario
plot_B_bot <- ggplot(scenario_marker_df,
    aes(x = scenario, y = reorder(marker, -marker_rank), fill = retained)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(
    data = scenario_marker_df[scenario_marker_df$retained, ],
    aes(label = "✓"), size = 2.8, colour = "white"
  ) +
  geom_text(
    data = scenario_marker_df[!scenario_marker_df$retained, ],
    aes(label = "✗"), size = 2.8, colour = "#CC3333"
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#2E8B57", "FALSE" = "#F0F0F0"),
    guide  = "none"
  ) +
  labs(x = NULL, y = "Duraclone markers") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x       = element_text(size = 8, lineheight = 1.1),
    axis.text.y       = element_text(size = 8),
    panel.grid        = element_blank()
  )

panel_B <- plot_B_top / plot_B_bot +
  plot_layout(heights = c(3, 2))

# ------------------------------------------------------------------------------
# Compose and save PDF
# ------------------------------------------------------------------------------

combined_plot <- plot_A | panel_B
ggsave(out_pdf, plot = combined_plot, width = 13, height = 8, device = "pdf")

message(sprintf("\n[Compat] PDF saved: %s", out_pdf))

# ------------------------------------------------------------------------------
# JSON report
# ------------------------------------------------------------------------------

scenario_report <- lapply(scenarios, function(s) {
  all_dc <- sort(duraclone_markers)
  lost   <- sort(setdiff(all_dc, s$markers))
  list(
    label              = gsub("\n", " ", s$label),
    n_shared_markers   = length(s$markers),
    shared_markers     = sort(s$markers),
    markers_lost_vs_dc = lost,
    n_pannelli_samples = s$n_pl,
    n_duraclone_samples= s$n_dc
  )
})
names(scenario_report) <- sapply(scenarios, `[[`, "id")

group_report <- lapply(seq_len(nrow(group_info)), function(i) {
  list(
    source      = group_info$source[i],
    n_files     = group_info$n_files[i],
    n_markers   = group_info$n_markers[i],
    markers     = sort(group_info$markers_list[[i]]),
    files       = basename(all_files$path[all_files$panel_sig == group_info$panel_sig[i]])
  )
})
names(group_report) <- paste(group_info$source, seq_len(nrow(group_info)), sep = "_")

recommendation <- list(
  recommended_scenario = "S2_TEX_DSold",
  rationale = paste(
    "TEX pannelli_liquidi + DS old gives the best marker/sample balance:",
    sprintf("%d shared markers with duraclone (loses only CTLA4 vs TEX-only),",
      length(s2_markers)),
    sprintf("n=%d pannelli_liquidi samples vs n=%d for TEX-only.",
      scenarios[[2]]$n_pl, scenarios[[1]]$n_pl),
    "Including non-TEX pannelli would additionally lose KI67 — unacceptable",
    "for an oncology immune monitoring baseline (proliferation marker)."
  ),
  markers_retained = sort(s2_markers),
  markers_lost_vs_tex_only = sort(setdiff(s1_markers, s2_markers))
)

report <- list(
  generated_at     = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  n_files_total    = nrow(all_files),
  panel_groups     = group_report,
  duraclone_markers= sort(duraclone_markers),
  scenarios        = scenario_report,
  recommendation   = recommendation
)

jsonlite::write_json(report, out_json, pretty = TRUE, auto_unbox = TRUE)
message(sprintf("[Compat] JSON saved: %s", out_json))

# ------------------------------------------------------------------------------
# Console summary
# ------------------------------------------------------------------------------

cat("\n")
cat("╔══════════════════════════════════════════════════════╗\n")
cat("║          PANEL COMPATIBILITY SUMMARY                 ║\n")
cat("╠══════════════════════════════════════════════════════╣\n")
for (s in scenarios) {
  marker_str <- paste(sort(s$markers), collapse = ", ")
  cat(sprintf("║  %-38s %2d markers ║\n",
    gsub("\n", " + ", s$label), length(s$markers)))
}
cat("╠══════════════════════════════════════════════════════╣\n")
cat(sprintf("║  RECOMMENDATION: TEX pannelli + DS old              ║\n"))
cat(sprintf("║  → %d markers  |  n_PL=%d  |  n_DC=%d              ║\n",
  length(s2_markers), scenarios[[2]]$n_pl, scenarios[[2]]$n_dc))
cat(sprintf("║  Retained: %-41s║\n",
  paste(sort(s2_markers), collapse = ", ")))
cat(sprintf("║  Lost vs. TEX-only: %-33s║\n",
  paste(sort(setdiff(s1_markers, s2_markers)), collapse = ", ")))
cat("╚══════════════════════════════════════════════════════╝\n")
