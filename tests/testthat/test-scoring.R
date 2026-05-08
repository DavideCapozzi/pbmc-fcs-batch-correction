library(testthat)

.src("src", "functions", "utils.R")
.src("src", "functions", "scoring.R")

# baseline_dict: batch + population (pop_N) + mean_freq + sd_freq
make_baseline <- function() {
  data.frame(
    batch      = c("A",     "A",     "B",     "B"),
    population = c("pop_1", "pop_2", "pop_1", "pop_2"),
    mean_freq  = c(0.20,    0.10,    0.25,    0.12),
    sd_freq    = c(0.05,    0.03,    0.06,    0.02),
    stringsAsFactors = FALSE
  )
}

# patient_freq_matrix: wide format — sample_id, batch, pop_N columns
make_patient <- function(pop1, pop2, batch = "A") {
  df <- data.frame(
    sample_id   = "S1",
    batch       = batch,
    total_cells = 1000L,
    stringsAsFactors = FALSE
  )
  df$pop_1 <- pop1
  df$pop_2 <- pop2
  df
}

test_that("EXPANDED flag fires when z > threshold", {
  # S1 batch A, pop_1: z = (0.40 - 0.20) / 0.05 = 4.0 → EXPANDED
  result <- compute_patient_zscores(make_patient(0.40, 0.10), make_baseline(), z_threshold = 2.0)
  row <- result[result$population == "pop_1", ]
  expect_equal(row$clinical_flag, "EXPANDED")
})

test_that("DEPLETED flag fires when z < -threshold", {
  # z = (0.05 - 0.20) / 0.05 = -3.0 → DEPLETED
  result <- compute_patient_zscores(make_patient(0.05, 0.10), make_baseline(), z_threshold = 2.0)
  row <- result[result$population == "pop_1", ]
  expect_equal(row$clinical_flag, "DEPLETED")
})

test_that("NORMAL flag fires when |z| <= threshold", {
  # z = (0.20 - 0.20) / 0.05 = 0 → NORMAL
  result <- compute_patient_zscores(make_patient(0.20, 0.10), make_baseline(), z_threshold = 2.0)
  row <- result[result$population == "pop_1", ]
  expect_equal(row$clinical_flag, "NORMAL")
})

test_that("UNDETERMINED flag when sd_freq is near zero", {
  baseline_zero_sd <- make_baseline()
  baseline_zero_sd$sd_freq[baseline_zero_sd$batch == "A" & baseline_zero_sd$population == "pop_1"] <- 0
  result <- compute_patient_zscores(make_patient(0.30, 0.10), baseline_zero_sd, z_threshold = 2.0)
  row <- result[result$population == "pop_1", ]
  expect_equal(row$clinical_flag, "UNDETERMINED")
})

test_that("safe_div returns NA for near-zero denominator", {
  expect_true(is.na(safe_div(1.0,  0.0)))
  expect_true(is.na(safe_div(1.0,  1e-12)))
  expect_equal(safe_div(4.0, 2.0), 2.0)
  expect_equal(safe_div(c(4, 0), c(2, 0), default = -99), c(2, -99))
})
