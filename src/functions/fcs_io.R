# src/functions/fcs_io.R
# ==============================================================================
# I/O MODULE FOR FLOW CYTOMETRY DATA
# Description: Handles reading, parsing, and preprocessing of FCS files.
# Dependencies: flowCore, dplyr, stringr
# ==============================================================================

library(flowCore)
library(dplyr)
library(stringr)

# ------------------------------------------------------------------------------
# Helper: Determine Batch ID from Filename
# ------------------------------------------------------------------------------
get_batch_id <- function(filename, patterns) {
  fname <- basename(filename)
  for (lbl in names(patterns)) {
    if (grepl(patterns[[lbl]], fname, ignore.case = TRUE)) {
      return(patterns[[lbl]])
    }
  }
  return("Unknown")
}

# ------------------------------------------------------------------------------
# Helper: Detect Shared Markers from Headers
# ------------------------------------------------------------------------------
#' @title Detect Shared Markers
#' @description Scans headers of all FCS files to find common biological markers.
#' @param files Character vector of file paths.
#' @param exclude_list Character vector of channel names/descriptions to ignore.
#' @return Character vector of shared markers.
detect_shared_markers <- function(files, exclude_list = NULL) {
  message("[IO] Scanning file headers to detect shared markers...")
  
  marker_lists <- lapply(files, function(f) {
    # Read ONLY the header for speed
    ff <- tryCatch({
      flowCore::read.FCS(f, which.lines = 1, transformation = FALSE, dataset = 1)
    }, error = function(e) {
      warning(sprintf("[IO Warning] Could not read header of %s", basename(f)))
      return(NULL)
    })
    
    if (is.null(ff)) return(NULL)
    
    pdata <- flowCore::pData(flowCore::parameters(ff))
    
    # Prioritize Description ($desc), fall back to Name ($name)
    mk_names <- pdata$desc
    mk_names[is.na(mk_names)] <- pdata$name[is.na(mk_names)]
    
    # Basic cleaning: remove NAs and empty strings
    mk_names <- mk_names[!is.na(mk_names) & mk_names != ""]
    return(mk_names)
  })
  
  # Remove NULLs
  marker_lists <- marker_lists[!sapply(marker_lists, is.null)]
  
  if (length(marker_lists) == 0) stop("[IO Fatal] No valid FCS headers read.")
  
  # Compute Intersection
  shared <- Reduce(intersect, marker_lists)
  
  # Filter exclusions
  if (!is.null(exclude_list)) {
    shared <- setdiff(shared, exclude_list)
  }
  
  message(sprintf("    -> Found %d shared markers across %d files.", length(shared), length(files)))
  
  if (length(shared) < 3) {
    warning("[IO Warning] Very few shared markers found (<3). Check panel compatibility or exclusions.")
  }
  
  return(shared)
}

# ------------------------------------------------------------------------------
# Core: Read Single FCS File
# ------------------------------------------------------------------------------
read_single_fcs <- function(file_path, markers_to_keep, cofactor = 5) {
  
  # Read full file
  ff <- tryCatch({
    flowCore::read.FCS(file_path, transformation = FALSE, truncate_max_range = FALSE)
  }, error = function(e) {
    warning(sprintf("[IO Error] Failed to read %s: %s", basename(file_path), e$message))
    return(NULL)
  })
  
  if (is.null(ff)) return(NULL)
  
  # Extract Data
  exprs_mat <- flowCore::exprs(ff)
  pdata     <- flowCore::pData(flowCore::parameters(ff))
  
  # Create mapping: Channel Name -> Marker Description
  channel_map <- setNames(pdata$desc, pdata$name)
  channel_map[is.na(channel_map)] <- names(channel_map)[is.na(channel_map)]
  
  # Rename columns to Marker Descriptions
  colnames(exprs_mat) <- channel_map[colnames(exprs_mat)]
  
  # Subset to shared markers
  valid_cols <- intersect(colnames(exprs_mat), markers_to_keep)
  
  if (length(valid_cols) == 0) {
    warning(sprintf("[IO Error] No shared markers found in %s. Skipping.", basename(file_path)))
    return(NULL)
  }
  
  exprs_mat <- exprs_mat[, valid_cols, drop = FALSE]
  
  # Arcsinh Transformation
  if (!is.null(cofactor) && cofactor > 0) {
    exprs_mat <- asinh(exprs_mat / cofactor)
  }
  
  return(as.data.frame(exprs_mat))
}

# ------------------------------------------------------------------------------
# Main: Batch Load Directory
# ------------------------------------------------------------------------------
#' @title Read and Prepare All FCS Data
#' @description Orchestrates the loading of a directory of FCS files.
#' @param data_dir Path to raw data.
#' @param patterns Named list for batch detection.
#' @param cofactor Transformation cofactor.
#' @param exclude List of markers to exclude.
#' @return A tibble with concatenated data.
read_and_prep_data <- function(data_dir, patterns, cofactor = 5, exclude = NULL) {
  
  # 1. File Discovery
  fcs_files <- list.files(data_dir, pattern = "\\.fcs$", full.names = TRUE, recursive = TRUE)
  if (length(fcs_files) == 0) stop("[IO Fatal] No FCS files found.")
  
  # 2. Dynamic Marker Detection
  shared_markers <- detect_shared_markers(fcs_files, exclude)
  
  message(sprintf("[IO] Importing data using %d shared markers...", length(shared_markers)))
  message(sprintf("    -> Markers: %s...", paste(head(shared_markers, 5), collapse=", ")))
  
  # 3. Iterative Loading
  data_list <- lapply(seq_along(fcs_files), function(i) {
    f <- fcs_files[i]
    if (i %% 10 == 0) message(sprintf("    -> Processing file %d/%d...", i, length(fcs_files)))
    
    df <- read_single_fcs(f, shared_markers, cofactor)
    
    if (is.null(df)) return(NULL)
    
    # Attach Metadata
    df$sample_id <- basename(f)
    df$batch     <- get_batch_id(f, patterns)
    
    return(df)
  })
  
  # 4. Aggregation
  data_list <- data_list[!sapply(data_list, is.null)]
  
  if (length(data_list) == 0) stop("[IO Fatal] All files failed to import.")
  
  full_data <- dplyr::bind_rows(data_list)
  
  # Final Integrity Check
  if (any(full_data$batch == "Unknown")) {
    n_unk <- sum(full_data$batch == "Unknown")
    warning(sprintf("[IO Warning] %d cells detected with 'Unknown' batch ID.", n_unk))
  }
  
  return(as_tibble(full_data))
}
