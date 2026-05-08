# src/functions/auto_phenotyping.R
# ==============================================================================
# AUTO-PHENOTYPING MODULE
# Programmatic immunological labeling of matched cell populations from centroid
# expression data. Uses per-batch 1D k-means adaptive thresholds combined with
# a hierarchical T cell classification scheme.
#
# Pipeline:
#   1. compute_adaptive_thresholds()  — 1D k-means per marker per batch
#   2. classify_population()          — hierarchical immunological rules
#   3. compute_reliability_score()    — margin-based confidence in [0,1]
#   4. run_auto_phenotyping()         — orchestrator; returns report + labels
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. ADAPTIVE MARKER THRESHOLDS
# ------------------------------------------------------------------------------

#' Compute a positive/negative threshold per marker using 1D k-means (k=2).
#' For each marker column the threshold is set at the midpoint between the two
#' k-means cluster centres. A marker is declared uninformative (threshold = NA)
#' when the gap between centres is smaller than `uninformative_gap` arcsinh
#' units, preventing artificial splits on homogeneous populations.
#'
#' @param centroid_mat Numeric matrix [n_populations x n_markers].
#'   Row names = metacluster IDs; col names = marker names.
#' @param uninformative_gap numeric(1) Minimum separation between k-means
#'   centres to consider a marker informative. Default 0.5 arcsinh units.
#' @return Named numeric vector length n_markers. NA where uninformative.
compute_adaptive_thresholds <- function(centroid_mat, uninformative_gap = 0.5,
                                        kmeans_nstart = 10L, kmeans_iter_max = 50L) {
  stopifnot(is.matrix(centroid_mat), nrow(centroid_mat) >= 2L)

  apply(centroid_mat, 2L, function(vals) {
    valid <- vals[is.finite(vals)]
    if (length(unique(valid)) < 2L) return(NA_real_)
    km <- tryCatch(
      kmeans(valid, centers = 2L,
             nstart   = as.integer(kmeans_nstart),
             iter.max = as.integer(kmeans_iter_max)),
      error = function(e) NULL
    )
    if (is.null(km)) return(NA_real_)
    ctrs <- sort(km$centers[, 1L])
    if ((ctrs[2L] - ctrs[1L]) < uninformative_gap) return(NA_real_)
    mean(ctrs)
  })
}

# ------------------------------------------------------------------------------
# 2. SINGLE-POPULATION CLASSIFIER
# ------------------------------------------------------------------------------

#' Classify one population's immunophenotype using a hierarchical T cell scheme.
#'
#' Decision tree:
#'   Gate 1 — CD3:  T cell vs Non-T cell
#'   Gate 2 — CD8:  CD8+ cytotoxic vs CD4+ helper (CD3+CD8-)
#'   Gate 3 — CD45RA × CD28 (± CCR7): Naive / Central Memory / Effector Memory / TEMRA
#'   Gate 4 — Functional overlays (all lineages):
#'             PD1+CTLA4+ → Exhausted; CD137+ → Activated; KI67+ → Proliferating
#'
#' When a gate marker has threshold = NA (uninformative for this batch) the gate
#' is skipped and the uncertainty is reflected in the label text.
#'
#' @param expr Named numeric vector: centroid expression for one population.
#'   Names must match sanitized marker names (e.g. CD3, CD8, CD45RA ...).
#' @param thresholds Named numeric vector from compute_adaptive_thresholds().
#' @return Named list: $phenotype_label (character), $defining_markers
#'   (character vector of "+/-" tokens), $functional_tags (character vector).
classify_population <- function(expr, thresholds) {
  is_pos <- function(marker) {
    if (!marker %in% names(thresholds)) return(NA)
    thr <- thresholds[[marker]]
    val <- expr[[marker]]
    if (is.na(thr) || is.na(val)) return(NA)
    val > thr
  }

  defining <- character(0L)
  tags     <- character(0L)

  # Gate 1: T cell lineage
  cd3_pos <- is_pos("CD3")
  if (is.na(cd3_pos)) {
    base_label <- "Unknown (CD3 uninformative)"
  } else if (!cd3_pos) {
    base_label <- "Non-T cell"
    defining   <- c(defining, "CD3-")
  } else {
    defining <- c(defining, "CD3+")

    # Gate 2: CD8 vs CD4
    cd8_pos <- is_pos("CD8")
    if (is.na(cd8_pos)) {
      lineage_label <- "T cell (CD8 uninformative)"
    } else if (cd8_pos) {
      lineage_label <- "CD8+ T cell"
      defining      <- c(defining, "CD8+")
    } else {
      lineage_label <- "CD4+ T cell"
      defining      <- c(defining, "CD8-")
    }

    # Gate 3: differentiation state (CD45RA x CD28)
    cd45ra_pos <- is_pos("CD45RA")
    cd28_pos   <- is_pos("CD28")

    if (is.na(cd45ra_pos) || is.na(cd28_pos)) {
      diff_state <- "(differentiation undetermined)"
    } else if (cd45ra_pos && cd28_pos) {
      diff_state <- "Naive"
    } else if (!cd45ra_pos && cd28_pos) {
      diff_state <- "Central Memory"
    } else if (!cd45ra_pos && !cd28_pos) {
      diff_state <- "Effector Memory"
    } else {
      diff_state <- "TEMRA"
    }
    if (!is.na(cd45ra_pos)) defining <- c(defining, if (cd45ra_pos) "CD45RA+" else "CD45RA-")
    if (!is.na(cd28_pos))   defining <- c(defining, if (cd28_pos)   "CD28+"   else "CD28-")

    base_label <- if (grepl("undetermined", diff_state, fixed = TRUE)) {
      lineage_label
    } else {
      paste(diff_state, lineage_label)
    }
  }

  # Gate 4: functional overlays
  pd1_pos   <- is_pos("PD1")
  ctla4_pos <- is_pos("CTLA4")
  cd137_pos <- is_pos("CD137")
  ki67_pos  <- is_pos("KI67")

  if (!is.na(pd1_pos) && !is.na(ctla4_pos) && pd1_pos && ctla4_pos) {
    tags     <- c(tags, "Exhausted")
    defining <- c(defining, "PD1+", "CTLA4+")
  } else if (!is.na(pd1_pos) && pd1_pos) {
    tags     <- c(tags, "PD1+")
    defining <- c(defining, "PD1+")
  }
  if (!is.na(cd137_pos) && cd137_pos) {
    tags     <- c(tags, "Activated")
    defining <- c(defining, "CD137+")
  }
  if (!is.na(ki67_pos) && ki67_pos) {
    tags     <- c(tags, "Proliferating")
    defining <- c(defining, "KI67+")
  }

  full_label <- if (length(tags) > 0L) {
    paste0(base_label, " (", paste(tags, collapse = "/"), ")")
  } else {
    base_label
  }

  list(
    phenotype_label  = full_label,
    defining_markers = defining,
    functional_tags  = tags
  )
}

# ------------------------------------------------------------------------------
# 3. RELIABILITY SCORE
# ------------------------------------------------------------------------------

#' Compute a confidence score for a phenotype assignment in [0, 1].
#'
#' For each defining marker (e.g. "CD3+", "CD8-"), the margin is:
#'   margin_i = |centroid_value - threshold_i| / sd_across_all_populations
#' The reliability score is tanh(median(margins)), which maps 0 (on threshold)
#' to 0 and a large, consistent margin to ~1.
#'
#' @param expr Named numeric vector: centroid for one population.
#' @param thresholds Named numeric vector of adaptive thresholds (NA = uninformative).
#' @param defining_markers Character vector from classify_population()$defining_markers.
#'   Trailing +/- sign is stripped to recover the marker name.
#' @param all_centroids Numeric matrix [n_pop x n_markers]: full batch centroid
#'   matrix used to compute per-marker population-level SD.
#' @return numeric(1) in [0, 1]; NA if no defining markers have valid thresholds.
compute_reliability_score <- function(expr, thresholds, defining_markers, all_centroids) {
  marker_names <- unique(sub("[+-]$", "", defining_markers))
  marker_names <- intersect(marker_names, colnames(all_centroids))
  if (length(marker_names) == 0L) return(NA_real_)

  margins <- vapply(marker_names, function(m) {
    thr <- thresholds[[m]]
    val <- expr[[m]]
    if (is.na(thr) || is.na(val)) return(NA_real_)
    col_sd <- sd(all_centroids[, m], na.rm = TRUE)
    if (is.na(col_sd) || col_sd < 1e-9) return(NA_real_)
    abs(val - thr) / col_sd
  }, numeric(1L))

  valid <- margins[is.finite(margins)]
  if (length(valid) == 0L) return(NA_real_)
  tanh(median(valid))
}

# ------------------------------------------------------------------------------
# 4. ORCHESTRATOR
# ------------------------------------------------------------------------------

#' Run auto-phenotyping on all matched populations across all batches.
#'
#' For each batch:
#'   1. Compute adaptive thresholds via compute_adaptive_thresholds()
#'   2. For each matched population assigned to this batch, call
#'      classify_population() and compute_reliability_score()
#'   3. Assemble a per-row report data.frame
#' Then build a consensus population_labels vector (majority vote across batches;
#' discordant labels are flagged with "(batch-discordant)").
#'
#' @param centroids_list Named list [batch_name -> matrix(n_meta x n_markers)].
#'   Row names = metacluster IDs; col names = sanitized marker names.
#' @param match_table data.frame from match_metaclusters() with columns:
#'   batch_a, cluster_a, batch_b, cluster_b, matched_population_id, is_matched.
#' @param markers Character vector of marker names to use. If NULL, uses the
#'   intersection of all batch centroid column names.
#' @param uninformative_gap numeric(1) passed to compute_adaptive_thresholds().
#' @return list with two elements:
#'   $report — data.frame: population_id, batch, cluster_id, <marker cols>,
#'             phenotype_label, reliability_score, defining_markers
#'   $population_labels — named character vector (names = "pop_1" … "pop_N")
run_auto_phenotyping <- function(centroids_list,
                                 match_table,
                                 markers                   = NULL,
                                 uninformative_gap         = 0.5,
                                 low_reliability_threshold = 0.30,
                                 kmeans_nstart             = 10L,
                                 kmeans_iter_max           = 50L) {
  stopifnot(
    is.list(centroids_list), length(centroids_list) >= 1L,
    is.data.frame(match_table),
    "matched_population_id" %in% colnames(match_table)
  )

  if (is.null(markers)) {
    markers <- Reduce(intersect, lapply(centroids_list, colnames))
  }
  if (length(markers) == 0L) stop("[AutoPheno] No shared markers across batches.")

  matched <- match_table[
    !is.na(match_table$matched_population_id) &
    !is.na(match_table$is_matched) &
    match_table$is_matched == TRUE, ]

  if (nrow(matched) == 0L) stop("[AutoPheno] No matched populations found in match_table.")

  # Build global thresholds by pooling centroids from all matched populations across all batches.
  # This ensures the same decision boundary is used regardless of per-batch MFI scale differences.
  avail_m_global <- Reduce(intersect, lapply(centroids_list, colnames))
  avail_m_global <- intersect(markers, avail_m_global)

  global_thr_rows <- list()
  for (row_idx in seq_len(nrow(matched))) {
    row_m <- matched[row_idx, ]
    for (bat_col in c("batch_a", "batch_b")) {
      clust_col <- if (bat_col == "batch_a") "cluster_a" else "cluster_b"
      bat       <- row_m[[bat_col]]
      cid       <- as.character(row_m[[clust_col]])
      if (!is.na(bat) && bat %in% names(centroids_list) && !is.na(cid)) {
        cmat_tmp <- centroids_list[[bat]]
        if (cid %in% rownames(cmat_tmp)) {
          global_thr_rows[[length(global_thr_rows) + 1L]] <-
            cmat_tmp[cid, avail_m_global, drop = FALSE]
        }
      }
    }
  }
  global_thresholds <- if (length(global_thr_rows) >= 2L) {
    combined_mat <- do.call(rbind, global_thr_rows)
    compute_adaptive_thresholds(combined_mat, uninformative_gap,
                                kmeans_nstart   = kmeans_nstart,
                                kmeans_iter_max = kmeans_iter_max)
  } else {
    setNames(rep(NA_real_, length(avail_m_global)), avail_m_global)
  }

  all_rows            <- list()
  labels_by_pop       <- list()

  for (b in names(centroids_list)) {
    cmat          <- centroids_list[[b]]
    avail_markers <- intersect(markers, colnames(cmat))

    if (length(avail_markers) == 0L) {
      warning(sprintf("[AutoPheno] Batch '%s': no shared markers — skipped.", b))
      next
    }

    for (i in seq_len(nrow(matched))) {
      row    <- matched[i, ]
      pop_id <- paste0("pop_", row$matched_population_id)

      if (!is.na(row$batch_a) && row$batch_a == b) {
        cluster_id <- as.character(row$cluster_a)
      } else if (!is.na(row$batch_b) && row$batch_b == b) {
        cluster_id <- as.character(row$cluster_b)
      } else {
        next
      }

      if (!cluster_id %in% rownames(cmat)) {
        warning(sprintf("[AutoPheno] Cluster '%s' absent in batch '%s' — skipped.",
                        cluster_id, b))
        next
      }

      expr_vec <- cmat[cluster_id, avail_markers]
      pheno    <- classify_population(expr_vec, global_thresholds[avail_markers])
      rel      <- compute_reliability_score(
        expr_vec, global_thresholds[avail_markers],
        pheno$defining_markers,
        cmat[, avail_markers, drop = FALSE]
      )

      label_out <- if (!is.na(rel) && rel < low_reliability_threshold) {
        paste0(pheno$phenotype_label, " [low confidence]")
      } else if (is.na(rel)) {
        paste0(pheno$phenotype_label, " [unclassified]")
      } else {
        pheno$phenotype_label
      }

      centroid_vals <- setNames(rep(NA_real_, length(markers)), markers)
      centroid_vals[avail_markers] <- round(expr_vec, 4L)

      all_rows[[length(all_rows) + 1L]] <- c(
        list(population_id = pop_id, batch = b, cluster_id = cluster_id),
        as.list(centroid_vals),
        list(
          phenotype_label   = label_out,
          reliability_score = round(rel %||% NA_real_, 4L),
          defining_markers  = paste(pheno$defining_markers, collapse = "; ")
        )
      )

      labels_by_pop[[pop_id]] <- c(labels_by_pop[[pop_id]], label_out)
    }
  }

  if (length(all_rows) == 0L) {
    stop("[AutoPheno] No populations could be classified. Check centroid data.")
  }

  report <- do.call(rbind, lapply(all_rows, as.data.frame, stringsAsFactors = FALSE))
  rownames(report) <- NULL

  population_labels <- vapply(names(labels_by_pop), function(pid) {
    labs <- labels_by_pop[[pid]]
    uniq <- unique(labs)
    if (length(uniq) == 1L) return(uniq)
    best <- names(sort(table(labs), decreasing = TRUE))[1L]
    paste0(best, " (batch-discordant)")
  }, character(1L))

  list(report = report, population_labels = population_labels)
}

# %||% is defined in src/functions/utils.R — sourced by main.R before all steps.
