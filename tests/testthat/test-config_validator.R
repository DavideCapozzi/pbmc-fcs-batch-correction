library(testthat)

.src("src", "functions", "utils.R")
.src("src", "functions", "config_validator.R")

make_valid_config <- function() {
  list(
    n_cores     = 4L,
    random_seed = 42L,
    directories = list(
      raw         = list(batch_a = "/tmp"),
      intermediate = "/tmp",
      processed   = "/tmp",
      figures     = "/tmp",
      logs        = "/tmp",
      integration = "/tmp"
    ),
    markers  = list(transform_cofactor = 150),
    qc       = list(enabled = FALSE),
    clustering = list(xdim = 10L, ydim = 10L, n_metaclusters = 20L, seed = 1234L,
                      run_per_batch = TRUE),
    integration   = list(similarity_threshold = 0.70),
    auto_phenotyping = list(uninformative_gap = 0.5, low_reliability_threshold = 0.30),
    scoring   = list(z_score_threshold = 2.0),
    dimred    = list(n_neighbors = 15L, min_dist = 0.1, seed = 1234L),
    parallel  = list(n_workers = NULL)
  )
}

test_that("valid config passes without error", {
  cfg <- make_valid_config()
  expect_no_error(validate_config(cfg))
})

test_that("missing required key is caught", {
  cfg <- make_valid_config()
  cfg$directories <- NULL
  expect_error(validate_config(cfg), regexp = "directories")
})

test_that("negative cofactor is rejected", {
  cfg <- make_valid_config()
  cfg$markers$transform_cofactor <- -5
  expect_error(validate_config(cfg), regexp = "transform_cofactor")
})

test_that("similarity_threshold out of [0,1] is rejected", {
  cfg <- make_valid_config()
  cfg$integration$similarity_threshold <- 1.5
  expect_error(validate_config(cfg), regexp = "similarity_threshold")
})

test_that("z_score_threshold <= 0 is rejected", {
  cfg <- make_valid_config()
  cfg$scoring$z_score_threshold <- 0
  expect_error(validate_config(cfg), regexp = "z_score_threshold")
})

test_that("n_cores is cast to integer", {
  cfg <- make_valid_config()
  cfg$n_cores <- "4"
  result <- validate_config(cfg)
  expect_true(is.integer(result$n_cores))
})

test_that("clustering integers are cast correctly", {
  cfg <- make_valid_config()
  cfg$clustering$xdim <- list(10L)  # YAML can produce lists
  result <- validate_config(cfg)
  expect_true(is.integer(result$clustering$xdim))
})
