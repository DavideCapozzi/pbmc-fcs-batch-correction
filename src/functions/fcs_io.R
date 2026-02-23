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
#' @description Orchestrates the loading of specific FCS directories. The directory name becomes the batch ID.
#' @param data_dirs Named list of paths to raw data (Name = Folder/Batch ID).
#' @param cofactor Transformation cofactor.
#' @param exclude List of markers to exclude.
#' @param test_mode Boolean to limit file ingestion.
#' @param test_limit Integer, max files per directory if test_mode is TRUE.
#' @return A tibble with concatenated data.
read_and_prep_data <- function(data_dirs, cofactor = 5, exclude = NULL, test_mode = FALSE, test_limit = 5) {
  
  # 1. File Discovery per Directory (Directory Name = Batch)
  all_files_df <- data.frame(file_path = character(), batch_id = character(), stringsAsFactors = FALSE)
  
  for (dir_name in names(data_dirs)) {
    dir_path <- data_dirs[[dir_name]]
    
    # Force recursive = FALSE to process only top-level files
    files_in_dir <- list.files(dir_path, pattern = "\\.fcs$", full.names = TRUE, recursive = FALSE)
    
    if (length(files_in_dir) == 0) {
      warning(sprintf("[IO Warning] No FCS files found in %s", dir_path))
      next
    }
    
    # Subsampling for test mode
    if (test_mode) {
      files_in_dir <- head(files_in_dir, test_limit)
      message(sprintf("[IO Test Mode] Taking first %d files from batch/folder: %s", length(files_in_dir), dir_name))
    }
    
    temp_df <- data.frame(file_path = files_in_dir, batch_id = dir_name, stringsAsFactors = FALSE)
    all_files_df <- rbind(all_files_df, temp_df)
  }
  
  fcs_files <- all_files_df$file_path
  if (length(fcs_files) == 0) stop("[IO Fatal] No FCS files found across all provided directories.")
  
  # 2. Dynamic Marker Detection
  shared_markers <- detect_shared_markers(fcs_files, exclude)
  
  message(sprintf("[IO] Importing data using %d shared markers...", length(shared_markers)))
  message(sprintf("    -> Markers: %s...", paste(head(shared_markers, 5), collapse=", ")))
  
  # 3. Iterative Loading
  data_list <- lapply(seq_along(fcs_files), function(i) {
    f <- fcs_files[i]
    b_id <- all_files_df$batch_id[i]
    if (i %% 5 == 0) message(sprintf("    -> Processing file %d/%d...", i, length(fcs_files)))
    
    df <- read_single_fcs(f, shared_markers, cofactor)
    
    if (is.null(df)) return(NULL)
    
    # Attach Metadata
    df$sample_id <- basename(f)
    df$batch     <- b_id # Folder name is now the batch directly
    
    return(df)
  })
  
  # 4. Aggregation
  data_list <- data_list[!sapply(data_list, is.null)]
  
  if (length(data_list) == 0) stop("[IO Fatal] All files failed to import.")
  
  full_data <- dplyr::bind_rows(data_list)
  
  return(as_tibble(full_data))
}
