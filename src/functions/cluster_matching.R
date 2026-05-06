# src/functions/cluster_matching.R
# ==============================================================================
# CLUSTER MATCHING MODULE
# Description: Cosine-similarity-based metacluster matching across batches.
#              Operates entirely on marker-centroid profiles, never on raw MFI.
# Dependencies: SingleCellExperiment
# ==============================================================================

library(SingleCellExperiment)

#' @title Compute Metacluster Centroids
#' @description Computes mean arcsinh-transformed expression per metacluster.
#' @param sce SingleCellExperiment with metacluster_id in colData and exprs assay.
#' @param markers Character vector of markers to include.
#' @return Numeric matrix [n_metaclusters x n_markers]; rownames = metacluster IDs.
compute_metacluster_centroids <- function(sce, markers) {

  stopifnot(
    "metacluster_id" %in% names(colData(sce)),
    "exprs" %in% assayNames(sce),
    all(markers %in% rownames(sce))
  )

  exprs_mat      <- t(assay(sce, "exprs")[markers, , drop = FALSE])
  metacluster_ids <- as.character(sce$metacluster_id)
  unique_ids      <- sort(unique(metacluster_ids))

  centroid_mat <- do.call(rbind, lapply(unique_ids, function(k) {
    cells_k <- which(metacluster_ids == k)
    if (length(cells_k) == 0L) return(rep(NA_real_, length(markers)))
    colMeans(exprs_mat[cells_k, , drop = FALSE])
  }))

  rownames(centroid_mat) <- unique_ids
  colnames(centroid_mat) <- markers

  centroid_mat
}

#' @title Compute Cosine Similarity Matrix
#' @description Pairwise cosine similarity between two L2-normalized centroid matrices.
#' @param centroids_a Numeric matrix [k_A x m].
#' @param centroids_b Numeric matrix [k_B x m].
#' @return Numeric matrix [k_A x k_B].
compute_cosine_similarity <- function(centroids_a, centroids_b) {
  row_norm <- function(mat) {
    norms <- sqrt(rowSums(mat^2))
    norms[norms == 0] <- 1e-10
    mat / norms
  }
  row_norm(centroids_a) %*% t(row_norm(centroids_b))
}

#' @title Match Metaclusters Across Batches
#' @description Greedy bipartite matching on cosine similarity. Only two batches
#'   are supported (the two cohorts). Pairs below similarity_threshold are kept
#'   in the table with is_matched = FALSE.
#' @param centroids_list Named list of centroid matrices, one per batch.
#' @param markers Character vector of markers used for centroid computation.
#' @param similarity_threshold Numeric in [0,1]. Default 0.70.
#' @return data.frame: batch_a, cluster_a, batch_b, cluster_b,
#'   cosine_similarity, matched_population_id, is_matched.
match_metaclusters <- function(centroids_list, markers, similarity_threshold = 0.70) {

  stopifnot(length(centroids_list) == 2L)

  batch_names <- names(centroids_list)
  batch_a     <- batch_names[1L]
  batch_b     <- batch_names[2L]
  C_a         <- centroids_list[[batch_a]]
  C_b         <- centroids_list[[batch_b]]

  valid_a <- !apply(C_a, 1L, anyNA)
  valid_b <- !apply(C_b, 1L, anyNA)
  C_a_clean <- C_a[valid_a, , drop = FALSE]
  C_b_clean <- C_b[valid_b, , drop = FALSE]

  sim_mat <- compute_cosine_similarity(C_a_clean, C_b_clean)

  assigned_a <- integer(0L)
  assigned_b <- integer(0L)
  records    <- list()
  pop_id     <- 1L

  flat_order <- order(sim_mat, decreasing = TRUE)
  for (idx in flat_order) {
    ij <- arrayInd(idx, dim(sim_mat))
    i  <- ij[1L]
    j  <- ij[2L]
    if (i %in% assigned_a || j %in% assigned_b) next

    sim_val    <- sim_mat[i, j]
    is_matched <- sim_val >= similarity_threshold
    records[[length(records) + 1L]] <- data.frame(
      batch_a              = batch_a,
      cluster_a            = rownames(C_a_clean)[i],
      batch_b              = batch_b,
      cluster_b            = rownames(C_b_clean)[j],
      cosine_similarity    = sim_val,
      matched_population_id = if (is_matched) pop_id else NA_integer_,
      is_matched           = is_matched,
      stringsAsFactors     = FALSE
    )
    if (is_matched) pop_id <- pop_id + 1L
    assigned_a <- c(assigned_a, i)
    assigned_b <- c(assigned_b, j)
  }

  # Clusters with no partner
  unmatched_a <- setdiff(seq_len(nrow(C_a_clean)), assigned_a)
  for (i in unmatched_a) {
    records[[length(records) + 1L]] <- data.frame(
      batch_a = batch_a, cluster_a = rownames(C_a_clean)[i],
      batch_b = NA_character_, cluster_b = NA_character_,
      cosine_similarity = NA_real_, matched_population_id = NA_integer_,
      is_matched = FALSE, stringsAsFactors = FALSE
    )
  }
  unmatched_b <- setdiff(seq_len(nrow(C_b_clean)), assigned_b)
  for (j in unmatched_b) {
    records[[length(records) + 1L]] <- data.frame(
      batch_a = NA_character_, cluster_a = NA_character_,
      batch_b = batch_b, cluster_b = rownames(C_b_clean)[j],
      cosine_similarity = NA_real_, matched_population_id = NA_integer_,
      is_matched = FALSE, stringsAsFactors = FALSE
    )
  }

  do.call(rbind, records)
}

#' @title Validate Cluster Match Quality
#' @description Computes per-marker absolute delta (primary metric) and
#'   Jensen-Shannon Divergence (secondary/legacy metric) for each matched pair.
#'   Per-marker delta is interpretable in arcsinh units and robust to cross-batch
#'   MFI scale shifts; JSD on softmax centroids is sensitive to absolute values
#'   and retained for backward compatibility only.
#' @param match_table data.frame from match_metaclusters().
#' @param centroids_list Named list of centroid matrices.
#' @param markers Character vector of markers.
#' @param delta_threshold numeric(1) max_marker_delta threshold for delta_pass
#'   (default 2.0 arcsinh units). Values < 1.0 = negligible, 1.0-2.0 = moderate.
#' @return List: $summary (data.frame with delta + JSD columns), $jsd_values
#'   (named numeric), $unmatched_clusters (data.frame).
validate_cluster_matches <- function(match_table, centroids_list, markers,
                                     delta_threshold = 2.0) {

  softmax <- function(x) {
    x <- x - max(x, na.rm = TRUE)
    exp_x <- exp(x)
    exp_x / sum(exp_x, na.rm = TRUE)
  }

  kl_div <- function(p, q) {
    eps <- 1e-10
    p   <- pmax(p, eps)
    q   <- pmax(q, eps)
    sum(p * log(p / q))
  }

  jsd <- function(p, q) {
    m       <- (p + q) / 2
    jsd_val <- 0.5 * kl_div(p, m) + 0.5 * kl_div(q, m)
    jsd_val / log(2)
  }

  batch_names <- names(centroids_list)
  batch_a     <- batch_names[1L]
  batch_b     <- batch_names[2L]

  matched_rows <- match_table[!is.na(match_table$matched_population_id), ]

  jsd_values   <- setNames(numeric(nrow(matched_rows)), matched_rows$matched_population_id)
  summary_rows <- vector("list", nrow(matched_rows))

  for (k in seq_len(nrow(matched_rows))) {
    row <- matched_rows[k, ]
    c_a <- centroids_list[[batch_a]][row$cluster_a, markers]
    c_b <- centroids_list[[batch_b]][row$cluster_b, markers]

    if (anyNA(c_a) || anyNA(c_b)) {
      jsd_val    <- NA_real_
      mean_delta <- NA_real_
      max_delta  <- NA_real_
      marker_max <- NA_character_
    } else {
      jsd_val    <- jsd(softmax(c_a), softmax(c_b))
      delta_vec  <- abs(c_a - c_b)
      mean_delta <- mean(delta_vec, na.rm = TRUE)
      max_delta  <- max(delta_vec,  na.rm = TRUE)
      marker_max <- markers[which.max(delta_vec)]
    }

    jsd_values[k]   <- jsd_val
    summary_rows[[k]] <- data.frame(
      matched_population_id = row$matched_population_id,
      cluster_a             = row$cluster_a,
      cluster_b             = row$cluster_b,
      cosine_similarity     = row$cosine_similarity,
      mean_marker_delta     = round(mean_delta, 4L),
      max_marker_delta      = round(max_delta,  4L),
      marker_with_max_delta = marker_max,
      delta_pass            = !is.na(max_delta) && max_delta < delta_threshold,
      jsd                   = jsd_val,
      jsd_pass              = !is.na(jsd_val) && jsd_val < 0.1,
      stringsAsFactors      = FALSE
    )
  }

  summary_df <- do.call(rbind, summary_rows)
  unmatched  <- match_table[!isTRUE(match_table$is_matched), ]

  list(
    summary            = summary_df,
    jsd_values         = jsd_values,
    unmatched_clusters = unmatched
  )
}
