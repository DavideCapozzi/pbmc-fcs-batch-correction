# scripts/diagnose_fcs_panels.R
# ==============================================================================
# FLUOROCHROME-AWARE PANEL COMPATIBILITY DIAGNOSTIC
#
# Reads FCS headers only (no cell data — runs in ~30 seconds for 100 files).
# Groups files by filename bracket-prefix (e.g. "[TEX]") or source directory,
# then sub-splits each group by unique marker signature so every distinct panel
# configuration gets its own column in the output.
#
# Output (results/diagnostics/panels_{timestamp}/):
#   fcs_panels_{timestamp}.pdf   — 3 pages
#       Page 1: Fluorochrome Heatmap (one column per panel sub-group)
#       Page 2: Panel Group & Sample Inventory (files listed per group)
#       Page 3: Cross-Batch Marker Integration Readiness (per-group columns)
#   fcs_panels_{timestamp}.xlsx  — 5 sheets
#   fcs_panels_{timestamp}.json  — machine-readable full report
#
# Usage (from project root):
#   Rscript scripts/diagnose_fcs_panels.R
# ==============================================================================

suppressPackageStartupMessages({
  library(flowCore)
  library(yaml)
  library(ggplot2)
  library(writexl)
  library(jsonlite)
})

source("src/functions/fcs_io.R")   # provides sanitize_biological_markers()

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==============================================================================
# Style constants
# ==============================================================================

COL_HDR      <- "#2E4057"   # table header background
COL_ACTIVE   <- "#1A7A3A"   # ACTIVE / SHARED_SAME_DYE
COL_EXCLUDED <- "#CC2222"   # EXCLUDED / INTRA_CONFLICT
COL_INACTIVE <- "#707070"   # INACTIVE / BATCH_SPECIFIC
COL_WARN     <- "#CC8800"   # SHARED_DIFF_DYE
COL_CONFLICT <- "#D73027"   # within-batch fluorochrome conflict (heatmap)
COL_MISTO    <- "#FC8D62"   # intra-group inconsistency (heatmap)
COL_OK       <- "#2E8B57"   # present and consistent

FILL_PAL <- c(
  CONFLICT = COL_CONFLICT, MISTO = COL_MISTO, OK = COL_OK,
  ABSENT = "#F5F5F5", INACTIVE = "#DADADA", EXCLUDED = "#B8B8B8")
TEXT_PAL <- c(
  CONFLICT = "white", MISTO = "white", OK = "white",
  ABSENT = "#BBBBBB", INACTIVE = "#2A2A2A", EXCLUDED = "#111111")

CAT_COLOURS <- c(
  SHARED_SAME_DYE = COL_ACTIVE,
  SHARED_DIFF_DYE = COL_WARN,
  BATCH_SPECIFIC  = COL_INACTIVE,
  INTRA_CONFLICT  = COL_EXCLUDED)

HDR_BG <- c(ACTIVE = COL_HDR, INACTIVE = "#606060", EXCLUDED = "#883333")

# ==============================================================================
# Config
# ==============================================================================

config            <- yaml::read_yaml("config/global_params.yml")
raw_dirs          <- config$directories$raw
raw_filters_cfg   <- config$directories$raw_filters   %||% list()
exclude_files_cfg <- config$directories$exclude_files %||% list()
batch_labels_cfg  <- config$batch_labels              %||% list()

ts          <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir_run <- file.path("results/diagnostics", paste0("panels_", ts))
if (!dir.exists(out_dir_run)) dir.create(out_dir_run, recursive = TRUE)
out_pdf  <- file.path(out_dir_run, paste0("fcs_panels_", ts, ".pdf"))
out_xlsx <- file.path(out_dir_run, paste0("fcs_panels_", ts, ".xlsx"))
out_json <- file.path(out_dir_run, paste0("fcs_panels_", ts, ".json"))

message("=== FCS PANEL FLUOROCHROME DIAGNOSTIC ===")

# ==============================================================================
# Helpers
# ==============================================================================

safe_kw <- function(kw, key) {
  val <- kw[[key]]
  if (is.null(val) || length(val) == 0L) NA_character_ else as.character(val[[1L]])
}

detect_file_group <- function(fname, dir_key) {
  m <- regmatches(fname, regexpr("^\\[[^\\]]+\\]", fname))
  if (length(m) > 0L && nchar(m) > 0L) m else dir_key
}

short_label <- function(s, max_len = 14L) {
  ifelse(nchar(s) > max_len, paste0(substr(s, 1L, max_len - 2L), ".."), s)
}

# ==============================================================================
# Step 1 — File inventory
# ==============================================================================

scan_dir_full <- function(dir_key, dir_path, raw_filter = NULL,
                           exclude_files_vec = NULL, batch_label = NULL) {
  if (!dir.exists(dir_path)) {
    warning(sprintf("[Panels] Directory not found: %s", dir_path))
    return(NULL)
  }
  fcs_paths <- list.files(dir_path, pattern = "\\.fcs$",
                           full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  if (length(fcs_paths) == 0L) return(NULL)
  fnames        <- basename(fcs_paths)
  in_exclude    <- if (!is.null(exclude_files_vec)) fnames %in% exclude_files_vec else rep(FALSE, length(fnames))
  passes_filter <- if (!is.null(raw_filter)) grepl(raw_filter, fnames) else rep(TRUE, length(fnames))
  status <- ifelse(in_exclude, "EXCLUDED", ifelse(passes_filter, "ACTIVE", "INACTIVE"))
  data.frame(file = fnames, file_path = fcs_paths, dir_key = dir_key,
             batch_label = batch_label %||% dir_key, status = status,
             stringsAsFactors = FALSE)
}

message("[Panels] Scanning directories...")
all_files_list <- Filter(Negate(is.null), lapply(names(raw_dirs), function(dk) {
  scan_dir_full(dk, raw_dirs[[dk]], raw_filters_cfg[[dk]],
                unlist(exclude_files_cfg[[dk]]), batch_labels_cfg[[dk]])
}))
if (length(all_files_list) == 0L) {
  message("[Panels] No FCS files found. Configured paths:")
  for (dk in names(raw_dirs)) message(sprintf("    %s: %s", dk, raw_dirs[[dk]]))
  quit(save = "no", status = 0)
}
all_files <- do.call(rbind, all_files_list)
rownames(all_files) <- NULL
all_files$group_label <- mapply(detect_file_group, all_files$file, all_files$dir_key,
                                 USE.NAMES = FALSE)
message(sprintf("[Panels] Found %d FCS files: %d ACTIVE | %d INACTIVE | %d EXCLUDED",
  nrow(all_files), sum(all_files$status == "ACTIVE"),
  sum(all_files$status == "INACTIVE"), sum(all_files$status == "EXCLUDED")))

# ==============================================================================
# Step 2 — FCS header extraction
# ==============================================================================

read_channel_info <- function(fpath) {
  ff <- tryCatch(
    flowCore::read.FCS(fpath, which.lines = 1L, transformation = FALSE,
                       truncate_max_range = FALSE, emptyValue = FALSE),
    error = function(e) NULL)
  if (is.null(ff)) return(NULL)
  p  <- flowCore::pData(flowCore::parameters(ff))
  kw <- flowCore::keyword(ff)
  descs   <- ifelse(!is.na(p$desc) & nzchar(p$desc), p$desc, p$name)
  markers <- sanitize_biological_markers(descs)
  extract_fluoro <- function(d) {
    if (is.na(d) || !grepl(" ", d, fixed = TRUE)) return(NA_character_)
    gsub("-[AHW]$", "", trimws(sub("^[^ ]+ ", "", d)))
  }
  fluorochromes <- vapply(descs, extract_fluoro, character(1L), USE.NAMES = FALSE)
  list(
    channels = data.frame(channel = p$name, desc_raw = descs, marker = markers,
                          fluorochrome = fluorochromes, stringsAsFactors = FALSE),
    instrument = safe_kw(kw, "$CYT"),
    acq_date   = safe_kw(kw, "$DATE"),
    n_events   = suppressWarnings(as.integer(safe_kw(kw, "$TOT"))),
    n_channels = nrow(p))
}

message(sprintf("[Panels] Reading %d FCS headers...", nrow(all_files)))
header_info <- lapply(seq_len(nrow(all_files)), function(i) {
  h <- read_channel_info(all_files$file_path[i])
  message(sprintf("    [%d/%d] %s: %s", i, nrow(all_files), all_files$file[i],
    if (is.null(h)) "FAILED"
    else sprintf("%d ch | %s events | %s", h$n_channels,
      format(h$n_events %||% NA_integer_, big.mark = ","), h$instrument %||% "?")))
  h
})
all_files$instrument    <- vapply(header_info, function(h) if (is.null(h)) NA_character_ else h$instrument %||% NA_character_, character(1L))
all_files$acq_date      <- vapply(header_info, function(h) if (is.null(h)) NA_character_ else h$acq_date   %||% NA_character_, character(1L))
all_files$n_events      <- vapply(header_info, function(h) if (is.null(h)) NA_integer_  else as.integer(h$n_events   %||% NA_integer_), integer(1L))
all_files$n_channels    <- vapply(header_info, function(h) if (is.null(h)) NA_integer_  else as.integer(h$n_channels %||% NA_integer_), integer(1L))
all_files$n_bio_markers <- vapply(header_info, function(h) {
  if (is.null(h)) NA_integer_ else sum(!is.na(h$channels$marker))
}, integer(1L))

# ==============================================================================
# Step 2.5 — Panel sub-grouping by marker signature
# ==============================================================================
# Within each group_label, detect distinct marker configurations. If > 1 exists,
# rename to {group_label}_1, {group_label}_2, etc. (most files = _1). If only one
# configuration exists, the name is left unchanged (no suffix).

all_files$marker_sig <- vapply(seq_len(nrow(all_files)), function(i) {
  h <- header_info[[i]]
  if (is.null(h)) return(NA_character_)
  mks <- sort(unique(h$channels$marker[!is.na(h$channels$marker)]))
  paste(mks, collapse = "|")
}, character(1L))

for (gl in unique(all_files$group_label)) {
  idx        <- which(all_files$group_label == gl)
  sigs       <- all_files$marker_sig[idx]
  valid_sigs <- sigs[!is.na(sigs)]
  unique_sigs <- names(sort(table(valid_sigs), decreasing = TRUE))
  if (length(unique_sigs) > 1L) {
    sig_map <- setNames(paste0(gl, "_", seq_along(unique_sigs)), unique_sigs)
    new_labels <- sig_map[sigs]
    new_labels[is.na(new_labels)] <- paste0(gl, "_unknown")
    all_files$group_label[idx] <- new_labels
    message(sprintf("[Panels] Sub-split '%s' → %s", gl,
      paste(paste0(unique_sigs, " → ", sig_map), collapse = " | ")))
  }
}
message(sprintf("[Panels] Panel groups after sub-split: %s",
  paste(sort(unique(all_files$group_label)), collapse = ", ")))

# ==============================================================================
# Step 3 — Long-format channel data (biological channels only)
# ==============================================================================

channel_long <- do.call(rbind, Filter(Negate(is.null), lapply(seq_len(nrow(all_files)), function(i) {
  h <- header_info[[i]]
  if (is.null(h)) return(NULL)
  df             <- h$channels
  df$file        <- all_files$file[i]
  df$dir_key     <- all_files$dir_key[i]
  df$batch_label <- all_files$batch_label[i]
  df$group_label <- all_files$group_label[i]
  df$status      <- all_files$status[i]
  df
})))
channel_long <- channel_long[!is.na(channel_long$marker), ]

# ==============================================================================
# Step 4 — Fluorochrome matrix: one row per (marker, group_label)
# ==============================================================================

all_markers <- sort(unique(channel_long$marker))
all_groups  <- sort(unique(all_files$group_label))

fluoro_long <- do.call(rbind, lapply(all_markers, function(mk) {
  do.call(rbind, lapply(all_groups, function(gl) {
    rows    <- channel_long[channel_long$marker == mk & channel_long$group_label == gl, ]
    fluoros <- unique(rows$fluorochrome[!is.na(rows$fluorochrome)])
    fluoro_val <- if (length(fluoros) == 0L)      NA_character_
                  else if (length(fluoros) == 1L)  fluoros
                  else                             "MISTO"
    data.frame(marker = mk, group_label = gl, fluorochrome = fluoro_val,
               stringsAsFactors = FALSE)
  }))
}))

group_meta <- do.call(rbind, lapply(all_groups, function(gl) {
  rows       <- all_files[all_files$group_label == gl, ]
  grp_status <- if ("ACTIVE"   %in% rows$status) "ACTIVE"
                else if ("INACTIVE" %in% rows$status) "INACTIVE"
                else "EXCLUDED"
  data.frame(group_label = gl, batch_label = rows$batch_label[1L],
             dir_key = rows$dir_key[1L], grp_status = grp_status,
             n_files = nrow(rows), n_active = sum(rows$status == "ACTIVE"),
             stringsAsFactors = FALSE)
}))

# Marker count per group (for group inventory header)
group_meta$n_bio_markers <- vapply(group_meta$group_label, function(gl) {
  length(unique(channel_long$marker[channel_long$group_label == gl]))
}, integer(1L))

# Representative marker_sig per group
group_meta$marker_sig <- vapply(group_meta$group_label, function(gl) {
  s <- all_files$marker_sig[all_files$group_label == gl & !is.na(all_files$marker_sig)]
  if (length(s) == 0L) NA_character_ else s[1L]
}, character(1L))

# Dominant cytometer model per group (most frequent non-NA value)
group_meta$instrument <- vapply(group_meta$group_label, function(gl) {
  instr <- all_files$instrument[all_files$group_label == gl & !is.na(all_files$instrument)]
  if (length(instr) == 0L) NA_character_
  else names(sort(table(instr), decreasing = TRUE))[1L]
}, character(1L))

fluoro_long <- merge(fluoro_long, group_meta, by = "group_label", all.x = TRUE)

# ==============================================================================
# Step 5 — Conflict detection (ACTIVE groups, within same batch_label)
# ==============================================================================

detect_conflict_cells <- function(fl) {
  fl_a   <- fl[fl$grp_status == "ACTIVE" & !is.na(fl$fluorochrome) & fl$fluorochrome != "MISTO", ]
  result <- list()
  for (bl in unique(fl_a$batch_label)) {
    fl_b <- fl_a[fl_a$batch_label == bl, ]
    for (mk in unique(fl_b$marker)) {
      fl_mk <- fl_b[fl_b$marker == mk, ]
      if (length(unique(fl_mk$fluorochrome)) > 1L) {
        for (gl in fl_mk$group_label) {
          result[[length(result) + 1L]] <- data.frame(
            batch_label = bl, marker = mk, group_label = gl, stringsAsFactors = FALSE)
        }
      }
    }
  }
  if (length(result) > 0L) unique(do.call(rbind, result))
  else data.frame(batch_label = character(), marker = character(), group_label = character())
}

build_conflict_summary <- function(fl, conflict_cells) {
  empty <- data.frame(batch_label = character(), marker = character(),
                      group_A = character(), fluorochrome_A = character(),
                      group_B = character(), fluorochrome_B = character(),
                      stringsAsFactors = FALSE)
  if (nrow(conflict_cells) == 0L) return(empty)
  fl_a <- fl[fl$grp_status == "ACTIVE" & !is.na(fl$fluorochrome) & fl$fluorochrome != "MISTO", ]
  keys <- unique(conflict_cells[, c("batch_label", "marker")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    bl  <- keys$batch_label[i]; mk <- keys$marker[i]
    sub <- fl_a[fl_a$batch_label == bl & fl_a$marker == mk, ]
    sub <- sub[!duplicated(sub$group_label), ]
    if (nrow(sub) < 2L) return(NULL)
    data.frame(batch_label = bl, marker = mk,
               group_A = sub$group_label[1L], fluorochrome_A = sub$fluorochrome[1L],
               group_B = sub$group_label[2L], fluorochrome_B = sub$fluorochrome[2L],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

conflict_cells   <- detect_conflict_cells(fluoro_long)
conflict_summary <- build_conflict_summary(fluoro_long, conflict_cells)
message(sprintf("[Panels] Within-batch conflicts: %d across %d marker(s)",
  nrow(conflict_summary), length(unique(conflict_summary$marker))))

# ==============================================================================
# Step 5b — Cross-batch intersection analysis
# ==============================================================================

active_batches <- unique(group_meta$batch_label[group_meta$grp_status == "ACTIVE"])

batch_marker_fluoros <- do.call(rbind, lapply(active_batches, function(bl) {
  active_groups_in_batch <- group_meta$group_label[
    group_meta$batch_label == bl & group_meta$grp_status == "ACTIVE"]
  do.call(rbind, lapply(all_markers, function(mk) {
    rows    <- fluoro_long[as.character(fluoro_long$marker) == mk &
                           as.character(fluoro_long$group_label) %in% active_groups_in_batch, ]
    fluoros <- unique(rows$fluorochrome[!is.na(rows$fluorochrome) & rows$fluorochrome != "MISTO"])
    misto   <- any(!is.na(rows$fluorochrome) & rows$fluorochrome == "MISTO")
    present    <- length(fluoros) > 0L || misto
    fluoro_val <- if (misto)              "MISTO"
                  else if (length(fluoros) == 1L) fluoros
                  else if (length(fluoros) == 0L) NA_character_
                  else paste(sort(fluoros), collapse = "/")
    data.frame(batch_label = bl, marker = mk, present = present,
               fluorochrome = fluoro_val, stringsAsFactors = FALSE)
  }))
}))

cross_batch_summary <- do.call(rbind, lapply(all_markers, function(mk) {
  rows           <- batch_marker_fluoros[batch_marker_fluoros$marker == mk, ]
  in_all         <- all(rows$present)
  has_intra      <- any(!is.na(rows$fluorochrome) & rows$fluorochrome == "MISTO")
  fluoros_clean  <- unique(rows$fluorochrome[!is.na(rows$fluorochrome) & rows$fluorochrome != "MISTO"])
  same_fluoro    <- in_all && !has_intra && length(fluoros_clean) == 1L
  category <- if (has_intra)        "INTRA_CONFLICT"
              else if (!in_all)      "BATCH_SPECIFIC"
              else if (same_fluoro)  "SHARED_SAME_DYE"
              else                   "SHARED_DIFF_DYE"
  data.frame(marker = mk, in_all_batches = in_all, same_fluorochrome = same_fluoro,
             different_fluorochrome = in_all && !has_intra && !same_fluoro,
             category = category, stringsAsFactors = FALSE)
}))

cat_order_lvls <- c("SHARED_SAME_DYE", "SHARED_DIFF_DYE", "BATCH_SPECIFIC", "INTRA_CONFLICT")
cross_batch_summary <- cross_batch_summary[
  order(match(cross_batch_summary$category, cat_order_lvls), cross_batch_summary$marker), ]
rownames(cross_batch_summary) <- NULL

message(sprintf(
  "[Panels] Cross-batch (%d batches): %d SHARED_SAME_DYE | %d SHARED_DIFF_DYE | %d BATCH_SPECIFIC | %d INTRA_CONFLICT",
  length(active_batches),
  sum(cross_batch_summary$category == "SHARED_SAME_DYE"),
  sum(cross_batch_summary$category == "SHARED_DIFF_DYE"),
  sum(cross_batch_summary$category == "BATCH_SPECIFIC"),
  sum(cross_batch_summary$category == "INTRA_CONFLICT")))

# ==============================================================================
# Step 6 — Assign cell_type and display label
# ==============================================================================

conflict_key <- paste(conflict_cells$marker, conflict_cells$group_label)
fluoro_long$cell_type <- with(fluoro_long, {
  mk_gl <- paste(marker, group_label)
  ifelse(grp_status == "INACTIVE", "INACTIVE",
  ifelse(grp_status == "EXCLUDED", "EXCLUDED",
  ifelse(is.na(fluorochrome),      "ABSENT",
  ifelse(fluorochrome == "MISTO",  "MISTO",
  ifelse(mk_gl %in% conflict_key,  "CONFLICT", "OK")))))
})
fluoro_long$cell_label <- with(fluoro_long,
  ifelse(is.na(fluorochrome), "—",
  ifelse(fluorochrome == "MISTO", "MISTO", fluorochrome)))

# ==============================================================================
# Step 7 — Factor ordering for plotting
# ==============================================================================

marker_counts <- tapply(!is.na(fluoro_long$fluorochrome), fluoro_long$marker, sum)
marker_order  <- names(sort(marker_counts, decreasing = TRUE))
fluoro_long$marker <- factor(fluoro_long$marker, levels = rev(marker_order))

status_rank <- c(ACTIVE = 1L, INACTIVE = 2L, EXCLUDED = 3L)
go_df <- unique(fluoro_long[, c("group_label", "batch_label", "grp_status")])
go_df <- go_df[order(status_rank[go_df$grp_status], go_df$batch_label, go_df$group_label), ]
group_order <- go_df$group_label

fluoro_long$group_label <- factor(fluoro_long$group_label, levels = group_order)

n_per_group    <- setNames(group_meta$n_files,      group_meta$group_label)
instr_per_group <- setNames(
  vapply(group_order, function(gl) {
    v <- group_meta$instrument[group_meta$group_label == gl]
    if (length(v) == 0L || is.na(v)) "?" else v
  }, character(1L)),
  group_order)

col_labels <- setNames(
  paste0(group_order, "\n",
         instr_per_group[group_order], "\n",
         "[", go_df$grp_status, "]  n=", n_per_group[group_order]),
  group_order)
fluoro_long$col_label <- factor(col_labels[as.character(fluoro_long$group_label)],
                                 levels = col_labels[group_order])

batch_for_group  <- setNames(go_df$batch_label, go_df$group_label)[group_order]
batch_boundaries <- which(diff(match(batch_for_group, unique(batch_for_group))) != 0L) + 0.5

# ==============================================================================
# Page 1 — Fluorochrome Heatmap
# ==============================================================================

vlines <- if (length(batch_boundaries) > 0L) {
  geom_vline(xintercept = batch_boundaries, colour = "#222222", linewidth = 1.0, linetype = "dashed")
} else NULL

p1 <- ggplot(fluoro_long, aes(x = col_label, y = marker, fill = cell_type)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = cell_label, colour = cell_type), size = 2.8, fontface = "bold") +
  vlines +
  scale_fill_manual(values = FILL_PAL, name = "Status",
    labels = c(CONFLICT = "Within-batch conflict", MISTO = "Intra-group inconsistency",
               OK = "Present (consistent)", ABSENT = "Absent",
               INACTIVE = "Inactive (filtered)", EXCLUDED = "Excluded")) +
  scale_colour_manual(values = TEXT_PAL, guide = "none") +
  labs(title = "FCS Panel Fluorochrome Compatibility",
       subtitle = paste0(
         "Column headers: group name / cytometer / [status]  n  |  ",
         "GREEN = consistent  |  RED = within-batch conflict  |  ORANGE = intra-group mix  |  ",
         "LIGHT GREY = inactive  |  DARK GREY = excluded  |  fluorochrome shown when present in panel  |  ",
         "— = absent from panel  |  dashed line = batch boundary"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x     = element_text(size = 8, lineheight = 1.2),
        axis.text.y     = element_text(size = 9, face = "bold"),
        panel.grid      = element_blank(),
        legend.position = "bottom",
        legend.title    = element_text(size = 8.5),
        legend.text     = element_text(size = 8),
        legend.key.size = unit(0.45, "cm"),
        plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(size = 8.5, colour = "grey40"),
        plot.margin     = margin(12, 16, 12, 12))

# ==============================================================================
# Page 2 — Panel Group & Sample Inventory
# ==============================================================================

make_group_inventory_page <- function(af, gm, go, page = 1L, n_pages = 1L) {
  GROUP_H  <- 1.4    # header bar height
  FILE_H   <- 0.85   # height per file row
  GAP_H    <- 0.6    # gap between groups
  N_COLS   <- 3L     # files per row

  rects_df <- list()
  texts_df <- list()
  cur_y    <- 0

  for (gl in go) {
    meta     <- gm[gm$group_label == gl, ]
    all_f    <- sort(af$file[af$group_label == gl])
    n_f      <- length(all_f)
    n_rows   <- ceiling(n_f / N_COLS)
    hdr_col  <- HDR_BG[meta$grp_status]
    n_bio    <- meta$n_bio_markers

    instr_str <- if (is.na(meta$instrument)) "?" else meta$instrument
    hdr_txt <- sprintf("  %s   ·   %s   ·   batch: %s   ·   [%s]   ·   %d / %d files   ·   %d bio markers",
      gl, instr_str, meta$batch_label, meta$grp_status, meta$n_active, meta$n_files, n_bio)

    rects_df[[length(rects_df) + 1L]] <- data.frame(
      xmin = -0.05, xmax = N_COLS + 0.05,
      ymin = -(cur_y + GROUP_H), ymax = -cur_y,
      fill = hdr_col, stringsAsFactors = FALSE)

    texts_df[[length(texts_df) + 1L]] <- data.frame(
      x = 0.05, y = -(cur_y + GROUP_H / 2L),
      label = hdr_txt, size = 3.4, colour = "white",
      fontface = "bold", hjust = 0L, stringsAsFactors = FALSE)

    cur_y <- cur_y + GROUP_H

    for (fi in seq_along(all_f)) {
      row_i <- (fi - 1L) %/% N_COLS
      col_i <- (fi - 1L) %%  N_COLS
      fname <- all_f[fi]
      disp  <- if (nchar(fname) > 34L) paste0(substr(fname, 1L, 31L), "...") else fname
      txt_col <- switch(af$status[af$file == fname & af$group_label == gl][1L],
        ACTIVE   = "#1A1A1A",
        EXCLUDED = "#CC2222",
        "#777777")  # INACTIVE

      texts_df[[length(texts_df) + 1L]] <- data.frame(
        x = col_i + 0.08, y = -(cur_y + row_i * FILE_H + FILE_H / 2L),
        label = disp, size = 2.7, colour = txt_col,
        fontface = "plain", hjust = 0L, stringsAsFactors = FALSE)
    }
    cur_y <- cur_y + n_rows * FILE_H + GAP_H
  }

  rect_df <- do.call(rbind, rects_df)
  text_df <- do.call(rbind, texts_df)

  pg_subtitle <- "Files listed by panel sub-group, alphabetical order"
  if (n_pages > 1L)
    pg_subtitle <- sprintf("%s  (page %d / %d)", pg_subtitle, page, n_pages)

  ggplot() +
    geom_rect(data = rect_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              colour = NA) +
    scale_fill_identity() +
    geom_text(data = text_df,
              aes(x = x, y = y, label = label, size = I(size),
                  colour = colour, fontface = fontface, hjust = hjust)) +
    scale_colour_identity() +
    expand_limits(x = c(0, N_COLS), y = c(-cur_y + GAP_H / 2, 0.4)) +
    labs(title = sprintf("Panel Group & Sample Inventory  (%d groups, %d files total)",
                         length(go), sum(af$group_label %in% go)),
         subtitle = pg_subtitle, x = NULL, y = NULL) +
    theme_void(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 13, margin = margin(b = 6)),
          plot.subtitle = element_text(size = 8.5, colour = "grey40", margin = margin(b = 8)),
          plot.margin   = margin(12, 16, 12, 16))
}

# Split group_order across pages if total row count is large
GROUP_H <- 1.4; FILE_H <- 0.85; GAP_H <- 0.6; N_COLS <- 3L
group_y_heights <- vapply(group_order, function(gl) {
  n_f <- sum(all_files$group_label == gl)
  GROUP_H + ceiling(n_f / N_COLS) * FILE_H + GAP_H
}, numeric(1L))

max_y_per_page <- 46
page_assign    <- integer(length(group_order))
cur_pg <- 1L; cur_height <- 0
for (i in seq_along(group_order)) {
  if (cur_height + group_y_heights[i] > max_y_per_page && cur_height > 0) {
    cur_pg    <- cur_pg + 1L
    cur_height <- 0
  }
  page_assign[i] <- cur_pg
  cur_height <- cur_height + group_y_heights[i]
}
n_inv_pages <- max(page_assign)

inv_plots <- lapply(seq_len(n_inv_pages), function(pg) {
  grps <- group_order[page_assign == pg]
  make_group_inventory_page(all_files, group_meta, grps, page = pg, n_pages = n_inv_pages)
})

# ==============================================================================
# Shared table helper
# ==============================================================================

make_table_gg <- function(df, title, subtitle = NULL,
                           status_col    = "Status",
                           status_colours = c(ACTIVE = COL_ACTIVE, EXCLUDED = COL_EXCLUDED),
                           status_default = COL_INACTIVE,
                           text_size     = 3.0) {
  cols         <- colnames(df)
  n_rows       <- nrow(df)
  n_cols       <- length(cols)
  stat_col_idx <- match(status_col, cols)
  hdr_text_sz  <- text_size + 0.5

  hdr <- data.frame(row = 0L, col = seq_len(n_cols), label = cols,
                    is_hdr = TRUE, status_val = NA_character_, stringsAsFactors = FALSE)
  body <- do.call(rbind, lapply(seq_len(n_rows), function(r) {
    data.frame(row = r, col = seq_len(n_cols),
               label = vapply(df[r, ], as.character, character(1L)),
               is_hdr = FALSE,
               status_val = if (!is.na(stat_col_idx)) df[[status_col]][r] else NA_character_,
               stringsAsFactors = FALSE)
  }))
  tbl <- rbind(hdr, body)
  tbl$bg <- ifelse(tbl$is_hdr, COL_HDR,
              ifelse(tbl$row %% 2L == 0L, "#F5F5F5", "white"))
  tbl$fc <- ifelse(tbl$is_hdr, "white",
              ifelse(!is.na(tbl$status_val) & !is.na(stat_col_idx) & tbl$col == stat_col_idx,
                ifelse(tbl$label %in% names(status_colours),
                  status_colours[tbl$label], status_default),
              "#1A1A1A"))
  tbl$fw  <- ifelse(tbl$is_hdr, "bold", "plain")
  tbl$tsz <- ifelse(tbl$is_hdr, hdr_text_sz, text_size)

  ttl_parts <- c(title, subtitle)
  ttl_parts <- ttl_parts[!is.null(ttl_parts) & !is.na(ttl_parts)]

  p <- ggplot(tbl, aes(x = col, y = -row)) +
    geom_tile(aes(fill = bg), colour = "white", linewidth = 0.25) +
    geom_text(aes(label = label, colour = fc, fontface = fw, size = I(tsz))) +
    scale_fill_identity() +
    scale_colour_identity() +
    scale_x_continuous(expand = c(0.02, 0)) +
    labs(title = ttl_parts[1L],
         subtitle = if (length(ttl_parts) > 1L) ttl_parts[2L] else NULL,
         x = NULL, y = NULL) +
    theme_void(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 13, margin = margin(b = 4)),
          plot.subtitle = element_text(size = 8.5, colour = "grey40", margin = margin(b = 6)),
          plot.margin   = margin(12, 16, 12, 12))
  p
}

# ==============================================================================
# Page 3 — Cross-Batch Marker Integration Readiness
# ==============================================================================

make_batch_intersection_page <- function(fl, cbs, grp_order, active_bl) {
  ttl     <- "Cross-Batch Marker Integration Readiness"
  sub_ttl <- sprintf("Active batches: %s   |   GREEN = safe  |  YELLOW = assay effect  |  GREY = batch-specific  |  RED = intra-conflict",
                     paste(active_bl, collapse = " vs "))

  if (length(active_bl) < 2L) {
    return(ggplot(data.frame(x=0.5,y=0.5,l="Only one active batch — no cross-batch comparison possible."),
                  aes(x,y,label=l)) +
             geom_text(size=5, colour="grey50") + expand_limits(x=c(0,1),y=c(0,1)) +
             labs(title=ttl) + theme_void() +
             theme(plot.title=element_text(face="bold", size=13, margin=margin(b=12))))
  }

  # Build display: Marker | <group_1> | <group_2> | ... | Category
  df_disp <- data.frame(Marker = cbs$marker, stringsAsFactors = FALSE)
  for (pg in grp_order) {
    col_vals <- vapply(cbs$marker, function(mk) {
      v <- fl$cell_label[as.character(fl$marker) == mk & as.character(fl$group_label) == pg]
      if (length(v) > 0L) v[1L] else "—"
    }, character(1L))
    df_disp[[short_label(pg, 14L)]] <- col_vals
  }
  df_disp$Category <- cbs$category

  make_table_gg(df_disp,
    title      = ttl,
    subtitle   = sub_ttl,
    status_col = "Category",
    status_colours = CAT_COLOURS,
    status_default = COL_INACTIVE,
    text_size  = 3.5)
}

p3 <- make_batch_intersection_page(fluoro_long, cross_batch_summary, group_order, active_batches)

# ==============================================================================
# Save PDF (2 + n_inv_pages pages)
# ==============================================================================

n_pdf_pages <- 1L + n_inv_pages + 1L
message(sprintf("[Panels] Writing PDF: %s (%d pages)", out_pdf, n_pdf_pages))
grDevices::pdf(out_pdf, width = 14, height = 11)
print(p1)
for (pp in inv_plots) print(pp)
print(p3)
invisible(grDevices::dev.off())
message(sprintf("[Panels] PDF written (%d pages).", n_pdf_pages))

# ==============================================================================
# JSON output
# ==============================================================================

build_json_report <- function() {
  # Per-group fluorochrome lookup for per_marker structured field
  grp_fluoro_lookup <- function(mk) {
    vals <- vapply(group_order, function(pg) {
      v <- fluoro_long$fluorochrome[as.character(fluoro_long$marker) == mk &
                                    as.character(fluoro_long$group_label) == pg]
      if (length(v) > 0L && !is.na(v[1L])) v[1L] else NA_character_
    }, character(1L))
    as.list(vals)
  }

  grp_summary <- lapply(group_order, function(gl) {
    meta  <- group_meta[group_meta$group_label == gl, ]
    fl    <- fluoro_long[as.character(fluoro_long$group_label) == gl, ]
    mk_ch <- as.character(fl$marker)
    files_all    <- sort(all_files$file[all_files$group_label == gl])
    files_active <- sort(all_files$file[all_files$group_label == gl & all_files$status == "ACTIVE"])
    list(
      group_label            = gl,
      batch_label            = meta$batch_label,
      dir_key                = meta$dir_key,
      status                 = meta$grp_status,
      n_files                = meta$n_files,
      n_active               = meta$n_active,
      instrument             = meta$instrument %||% NA_character_,
      panel_signature        = meta$marker_sig %||% NA_character_,
      n_bio_markers          = meta$n_bio_markers,
      files                  = files_all,
      files_active           = files_active,
      markers_present        = sort(mk_ch[!is.na(fl$fluorochrome)]),
      markers_absent         = sort(mk_ch[is.na(fl$fluorochrome)]),
      markers_intra_conflict = sort(mk_ch[!is.na(fl$fluorochrome) & fl$fluorochrome == "MISTO"])
    )
  })

  batch_specific_list <- lapply(active_batches, function(bl) {
    bs_markers <- cross_batch_summary$marker[cross_batch_summary$category == "BATCH_SPECIFIC"]
    present_here <- vapply(bs_markers, function(mk) {
      idx <- batch_marker_fluoros$batch_label == bl & batch_marker_fluoros$marker == mk
      any(batch_marker_fluoros$present[idx])
    }, logical(1L))
    bs_markers[present_here]
  })
  names(batch_specific_list) <- active_batches

  cbi <- list(
    active_batches   = active_batches,
    n_active_batches = length(active_batches),
    shared_same_dye  = cross_batch_summary$marker[cross_batch_summary$category == "SHARED_SAME_DYE"],
    shared_diff_dye  = cross_batch_summary$marker[cross_batch_summary$category == "SHARED_DIFF_DYE"],
    intra_conflict   = cross_batch_summary$marker[cross_batch_summary$category == "INTRA_CONFLICT"],
    batch_specific   = batch_specific_list,
    per_marker = lapply(seq_len(nrow(cross_batch_summary)), function(i) {
      mk <- cross_batch_summary$marker[i]
      c(as.list(cross_batch_summary[i, c("marker","in_all_batches","same_fluorochrome",
                                          "different_fluorochrome","category")]),
        list(per_group = grp_fluoro_lookup(mk)))
    })
  )

  conflicts_list <- if (nrow(conflict_summary) > 0L) {
    lapply(seq_len(nrow(conflict_summary)), function(i) as.list(conflict_summary[i, ]))
  } else list()

  files_list <- lapply(seq_len(nrow(all_files)), function(i) {
    h <- header_info[[i]]
    channels_list <- if (!is.null(h)) {
      ch <- h$channels[!is.na(h$channels$marker), ]
      lapply(seq_len(nrow(ch)), function(j)
        list(channel = ch$channel[j], desc_raw = ch$desc_raw[j], marker = ch$marker[j],
             fluorochrome = if (is.na(ch$fluorochrome[j])) NULL else ch$fluorochrome[j]))
    } else list()
    list(file = all_files$file[i], group_label = all_files$group_label[i],
         batch_label = all_files$batch_label[i], status = all_files$status[i],
         instrument    = if (is.na(all_files$instrument[i]))    NULL else all_files$instrument[i],
         acq_date      = if (is.na(all_files$acq_date[i]))      NULL else all_files$acq_date[i],
         n_events      = if (is.na(all_files$n_events[i]))      NULL else all_files$n_events[i],
         n_channels    = if (is.na(all_files$n_channels[i]))    NULL else all_files$n_channels[i],
         n_bio_markers = if (is.na(all_files$n_bio_markers[i])) NULL else all_files$n_bio_markers[i],
         panel_signature = if (is.na(all_files$marker_sig[i]))  NULL else all_files$marker_sig[i],
         channels = channels_list)
  })
  names(files_list) <- all_files$file

  list(generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
       n_files = nrow(all_files), n_active = sum(all_files$status == "ACTIVE"),
       n_inactive = sum(all_files$status == "INACTIVE"), n_excluded = sum(all_files$status == "EXCLUDED"),
       n_markers = length(all_markers), n_groups = length(all_groups),
       n_batches = length(unique(group_meta$batch_label)), n_active_batches = length(active_batches),
       within_batch_conflicts = nrow(conflict_summary),
       group_summary = grp_summary, cross_batch_intersection = cbi,
       conflicts = conflicts_list, files = files_list)
}

message(sprintf("[Panels] Writing JSON: %s", out_json))
json_report <- build_json_report()
writeLines(jsonlite::toJSON(json_report, auto_unbox = TRUE, pretty = TRUE, na = "null"), out_json)
message("[Panels] JSON written.")

# ==============================================================================
# Excel output (5 sheets)
# ==============================================================================

# Sheet 1: fluorochrome matrix (wide, one column per group)
sheet1_rows <- lapply(marker_order, function(mk) {
  vals <- setNames(lapply(group_order, function(gl) {
    v <- fluoro_long$cell_label[as.character(fluoro_long$marker) == mk &
                                as.character(fluoro_long$group_label) == gl]
    if (length(v) == 0L) "—" else v[1L]
  }), group_order)
  as.data.frame(c(list(Marker = mk), vals), stringsAsFactors = FALSE, check.names = FALSE)
})
sheet1 <- do.call(rbind, sheet1_rows)

# Sheet 2: full file inventory
sheet2 <- all_files[order(match(all_files$status, c("ACTIVE","EXCLUDED","INACTIVE")),
                           all_files$batch_label, all_files$group_label, all_files$file),
                     c("file","group_label","batch_label","status",
                       "instrument","acq_date","n_events","n_channels","n_bio_markers")]
colnames(sheet2) <- c("File","Group","Batch","Status","Instrument","Date","N_Events","N_Channels","N_Bio_Markers")
rownames(sheet2) <- NULL

# Sheet 3: conflicts
if (nrow(conflict_summary) > 0L) {
  sheet3 <- conflict_summary
  colnames(sheet3) <- c("Batch","Marker","Group_A","Fluorochrome_A","Group_B","Fluorochrome_B")
  rownames(sheet3) <- NULL
} else {
  sheet3 <- data.frame(Message = "No within-batch fluorochrome conflicts detected.",
                       stringsAsFactors = FALSE)
}

# Sheet 4: per-file × biological channel detail
sheet4 <- channel_long[order(channel_long$batch_label, channel_long$group_label,
                              channel_long$file, channel_long$marker),
                        c("file","group_label","batch_label","status","channel","desc_raw","marker","fluorochrome")]
colnames(sheet4) <- c("File","Group","Batch","Status","Channel","Description_Raw","Marker","Fluorochrome")
rownames(sheet4) <- NULL

# Sheet 5: cross-batch integration readiness
sheet5 <- cross_batch_summary[, c("marker","in_all_batches","same_fluorochrome","category")]
colnames(sheet5) <- c("Marker","In_All_Batches","Same_Fluorochrome","Category")
# Add per-group fluorochrome columns
for (pg in group_order) {
  col_vals <- vapply(cross_batch_summary$marker, function(mk) {
    v <- fluoro_long$cell_label[as.character(fluoro_long$marker) == mk &
                                as.character(fluoro_long$group_label) == pg]
    if (length(v) > 0L) v[1L] else "—"
  }, character(1L))
  sheet5[[pg]] <- col_vals
}
rownames(sheet5) <- NULL

message(sprintf("[Panels] Writing Excel: %s", out_xlsx))
writexl::write_xlsx(
  list(fluoro_matrix = sheet1, file_inventory = sheet2, conflicts = sheet3,
       channel_map = sheet4, batch_intersection = sheet5),
  path = out_xlsx)
message("[Panels] Excel written (5 sheets).")

# ==============================================================================
# Console summary
# ==============================================================================

n_shared_same <- sum(cross_batch_summary$category == "SHARED_SAME_DYE")
n_shared_diff <- sum(cross_batch_summary$category == "SHARED_DIFF_DYE")
n_batch_spec  <- sum(cross_batch_summary$category == "BATCH_SPECIFIC")

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║          FCS PANEL FLUOROCHROME DIAGNOSTIC SUMMARY          ║\n")
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Files: %-3d total   ACTIVE %-3d | INACTIVE %-3d | EXCL %-3d  ║\n",
  nrow(all_files), sum(all_files$status == "ACTIVE"),
  sum(all_files$status == "INACTIVE"), sum(all_files$status == "EXCLUDED")))
cat(sprintf("║  Panel sub-groups: %-42s║\n",
  paste(sort(unique(all_files$group_label)), collapse = ", ")))
cat(sprintf("║  Biological markers: %-40d║\n", length(all_markers)))
cat(sprintf("║  Within-batch conflicts: %-36d║\n", nrow(conflict_summary)))
cat(sprintf("║  Cross-batch: %d SHARED_SAME_DYE | %d SHARED_DIFF_DYE | %d BATCH_SPECIFIC\n",
  n_shared_same, n_shared_diff, n_batch_spec))
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Dir:  %s\n", basename(out_dir_run)))
cat(sprintf("║  PDF:  %s  (%d pages)\n", basename(out_pdf), n_pdf_pages))
cat(sprintf("║  XLSX: %s  (5 sheets)\n", basename(out_xlsx)))
cat(sprintf("║  JSON: %s\n", basename(out_json)))
cat("╚══════════════════════════════════════════════════════════════╝\n")
