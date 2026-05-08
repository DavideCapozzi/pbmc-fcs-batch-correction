# src/functions/utils.R
# ==============================================================================
# SHARED UTILITIES
# Centralizes small helpers used across multiple pipeline modules.
# Source this file once in main.R before any step scripts.
# ==============================================================================

#' Null-coalescing operator.
#' Returns x if not NULL, otherwise y.
`%||%` <- function(x, y) if (!is.null(x)) x else y

#' Safe division with a guard for near-zero denominators.
#' @param num numeric vector — numerator
#' @param den numeric vector — denominator
#' @param default value returned when |den| <= eps (default NA_real_)
#' @param eps numeric tolerance for "zero" (default 1e-10)
#' @return numeric vector same length as num/den
safe_div <- function(num, den, default = NA_real_, eps = 1e-10) {
  ifelse(abs(den) > eps, num / den, default)
}
