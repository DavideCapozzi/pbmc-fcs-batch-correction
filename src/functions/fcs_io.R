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
# Helper: Semantic Normalizer for Cytometry Markers
# ------------------------------------------------------------------------------
#' @title Sanitize Biological Markers
#' @description Normalizes raw cytometer channel names to pure biological antigens.
#' @param raw_names Character vector of raw channel descriptions.
#' @return Character vector of cleaned antigen names (or NA for technical channels).
sanitize_biological_markers <- function(raw_names) {
  
  # 1. Flag physical, time, and redundant height/width channels
  regex_garbage <- "(?i)(^FSC|^SSC|^TIME|WIDTH|HEIGHT|LENGTH|\\-H$|\\-W$|\\-H\\b|\\-W\\b)"
  is_garbage <- grepl(regex_garbage, raw_names, perl = TRUE)
  
  # 2. Isolate the core biological antigen (remove fluorochrome after space)
  # Example: "Ki-67 FITC-A" -> "Ki-67", "CD3 APC-A700-A" -> "CD3"
  core_names <- sub(" .*$", "", raw_names)
  
  # 3. Standardize spelling variations (remove hyphens, force uppercase)
  # Example: "Ki-67" -> "KI67", "PD-1" -> "PD1", "Viakrome" -> "VIAKROME"
  clean_names <- toupper(gsub("-", "", core_names))
  
  # 4. Erase garbage channels completely
  clean_names[is_garbage] <- NA
  
  return(clean_names)
}

# ------------------------------------------------------------------------------
# Helper: Detect Shared Markers from Headers
# ------------------------------------------------------------------------------
detect_shared_markers <- function(files, exclude_list = NULL) {
  message("[IO] Scanning file headers to detect shared markers...")
  
  marker_lists <- lapply(files, function(f) {
    ff <- tryCatch({
      flowCore::read.FCS(f, which.lines = 1, transformation = FALSE, dataset = 1)
    }, error = function(e) {
      warning(sprintf("[IO Warning] Could not read header of %s", basename(f)))
      return(NULL)
    })
    
    if (is.null(ff)) return(NULL)
    
    pdata <- flowCore::pData(flowCore::parameters(ff))
    mk_names <- pdata$desc
    mk_names[is.na(mk_names)] <- pdata$name[is.na(mk_names)]
    
    # Semantic normalization to find pure biological intersections
    clean_names <- sanitize_biological_markers(mk_names)
    
    # Drop NAs (which are the physical and -H/-W channels)
    clean_names <- clean_names[!is.na(clean_names)]
    
    return(clean_names)
  })
  
  marker_lists <- marker_lists[!sapply(marker_lists, is.null)]
  if (length(marker_lists) == 0) stop("[IO Fatal] No valid FCS headers read.")
  
  shared <- Reduce(intersect, marker_lists)
  
  # YAML Explicit Exclusions
  if (!is.null(exclude_list)) {
    clean_exclude <- sanitize_biological_markers(unlist(exclude_list))
    shared <- shared[!shared %in% clean_exclude]
  }
  
  message(sprintf("    -> Found %d shared biological markers across %d files.", length(shared), length(files)))
  
  if (length(shared) < 3) {
    warning("[IO Warning] Very few shared markers found (<3). Check panel compatibility.")
  }
  
  return(shared)
}

# ------------------------------------------------------------------------------
# Core: Read Single FCS File
# ------------------------------------------------------------------------------
read_single_fcs <- function(file_path, markers_to_keep, cofactor = 5) {
  
  ff <- tryCatch({
    flowCore::read.FCS(file_path, transformation = FALSE, truncate_max_range = FALSE)
  }, error = function(e) {
    warning(sprintf("[IO Error] Failed to read %s: %s", basename(file_path), e$message))
    return(NULL)
  })
  
  if (is.null(ff)) return(NULL)
  
  exprs_mat <- flowCore::exprs(ff)
  pdata     <- flowCore::pData(flowCore::parameters(ff))
  
  channel_map <- setNames(pdata$desc, pdata$name)
  channel_map[is.na(channel_map)] <- names(channel_map)[is.na(channel_map)]
  
  # Extract raw descriptions and apply the exact same semantic normalization
  raw_desc <- channel_map[colnames(exprs_mat)]
  clean_desc <- sanitize_biological_markers(raw_desc)
  
  # Rename matrix columns to pure antigen names (e.g., "CD3", "KI67")
  colnames(exprs_mat) <- clean_desc
  
  # Remove the NA columns (the physical/garbage channels stripped by the normalizer)
  exprs_mat <- exprs_mat[, !is.na(colnames(exprs_mat)), drop = FALSE]
  
  # Subset to the shared biological markers
  valid_cols <- intersect(colnames(exprs_mat), markers_to_keep)
  
  if (length(valid_cols) == 0) {
    warning(sprintf("[IO Error] No shared markers found in %s. Skipping.", basename(file_path)))
    return(NULL)
  }
  
  exprs_mat <- exprs_mat[, valid_cols, drop = FALSE]
  
  if (!is.null(cofactor) && cofactor > 0) {
    exprs_mat <- asinh(exprs_mat / cofactor)
  }
  
  return(as.data.frame(exprs_mat))
}

# ------------------------------------------------------------------------------
# Main: Batch Load Directory
# ------------------------------------------------------------------------------
read_and_prep_data <- function(data_dirs, cofactor = 5, exclude = NULL, test_mode = FALSE, test_limit = 5) {
  
  all_files_df <- data.frame(file_path = character(), batch_id = character(), stringsAsFactors = FALSE)
  
  for (dir_name in names(data_dirs)) {
    dir_path <- data_dirs[[dir_name]]
    files_in_dir <- list.files(dir_path, pattern = "\\.fcs$", full.names = TRUE, recursive = FALSE)
    
    if (length(files_in_dir) == 0) next
    
    if (test_mode) {
      files_in_dir <- head(files_in_dir, test_limit)
    }
    
    temp_df <- data.frame(file_path = files_in_dir, batch_id = dir_name, stringsAsFactors = FALSE)
    all_files_df <- rbind(all_files_df, temp_df)
  }
  
  fcs_files <- all_files_df$file_path
  if (length(fcs_files) == 0) stop("[IO Fatal] No FCS files found.")
  
  shared_markers <- detect_shared_markers(fcs_files, exclude)
  
  message(sprintf("[IO] Importing data using %d shared markers...", length(shared_markers)))
  message(sprintf("    -> Markers: %s...", paste(head(shared_markers, 5), collapse=", ")))
  
  data_list <- lapply(seq_along(fcs_files), function(i) {
    f <- fcs_files[i]
    b_id <- all_files_df$batch_id[i]
    
    df <- read_single_fcs(f, shared_markers, cofactor)
    
    if (is.null(df)) return(NULL)
    
    message(sprintf("    -> [%d/%d] SUCCESS: %s | Batch: %s | Events: %d", 
                    i, length(fcs_files), basename(f), b_id, nrow(df)))
    
    df$sample_id <- basename(f)
    df$batch     <- b_id 
    
    return(df)
  })
  
  data_list <- data_list[!sapply(data_list, is.null)]
  if (length(data_list) == 0) stop("[IO Fatal] All files failed to import.")
  
  message("[IO] Aggregating final dataset...")
  full_data <- dplyr::bind_rows(data_list)
  
  return(as_tibble(full_data))
}