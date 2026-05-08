library(testthat)

.src("src", "functions", "utils.R")
.src("src", "functions", "cluster_matching.R")

# Helper: make a simple centroid matrix
make_centroids <- function(n_clusters, n_markers, seed = 1L) {
  set.seed(seed)
  mat <- matrix(rnorm(n_clusters * n_markers), nrow = n_clusters, ncol = n_markers)
  rownames(mat) <- as.character(seq_len(n_clusters))
  colnames(mat) <- paste0("M", seq_len(n_markers))
  mat
}

test_that("identical centroids yield cosine similarity = 1", {
  C <- make_centroids(3L, 5L)
  sim <- compute_cosine_similarity(C, C)
  expect_equal(unname(diag(sim)), rep(1, 3), tolerance = 1e-10)
})

test_that("match_metaclusters returns expected columns", {
  C_a <- make_centroids(4L, 5L, seed = 1L)
  C_b <- make_centroids(4L, 5L, seed = 2L)
  result <- match_metaclusters(list(batch_a = C_a, batch_b = C_b),
                               markers = paste0("M", 1:5),
                               similarity_threshold = 0.0)
  expected_cols <- c("batch_a", "cluster_a", "batch_b", "cluster_b",
                     "cosine_similarity", "matched_population_id", "is_matched")
  expect_true(all(expected_cols %in% colnames(result)))
})

test_that("no pairs matched when all similarity < threshold", {
  set.seed(42)
  # Force orthogonal-like centroids
  C_a <- diag(4)
  C_b <- -diag(4)
  rownames(C_a) <- rownames(C_b) <- as.character(1:4)
  colnames(C_a) <- colnames(C_b) <- paste0("M", 1:4)
  result <- match_metaclusters(list(batch_a = C_a, batch_b = C_b),
                               markers = paste0("M", 1:4),
                               similarity_threshold = 0.9)
  expect_true(all(!result$is_matched, na.rm = TRUE))
})

test_that("validate_cluster_matches respects configurable jsd_threshold", {
  # Identical centroids → JSD = 0, should pass any positive threshold
  C <- make_centroids(3L, 5L, seed = 10L)
  centroids_list <- list(batch_a = C, batch_b = C)
  match_table <- match_metaclusters(centroids_list,
                                    markers = paste0("M", 1:5),
                                    similarity_threshold = 0.0)
  val <- validate_cluster_matches(match_table, centroids_list,
                                  markers = paste0("M", 1:5),
                                  jsd_threshold = 0.001)
  if (nrow(val$summary) > 0L) {
    expect_true(all(val$summary$jsd_pass, na.rm = TRUE))
  }
})
