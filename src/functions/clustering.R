# src/functions/clustering.R
# ==============================================================================
# CLUSTERING MODULE
# Description: Wrapper for FlowSOM clustering and metaclustering.
# Dependencies: FlowSOM, ConsensusClusterPlus, SingleCellExperiment
# ==============================================================================

library(FlowSOM)
library(SingleCellExperiment)
library(ConsensusClusterPlus)

#' @title Run FlowSOM Clustering
#' @description Executes FlowSOM workflow: SOM grid training -> Metaclustering.
#' @param sce SingleCellExperiment object (corrected data).
#' @param markers Character vector of markers to use for clustering.
#' @param config List containing clustering params (xdim, ydim, n_metaclusters, seed).
#' @return Modified SCE object with 'cluster_id' and 'metacluster_id' in colData.
run_flowsom_clustering <- function(sce, markers, config) {
  
  message("[Clustering] Initializing FlowSOM...")
  
  # FIXED: Force integer types for critical configuration parameters.
  # This prevents "list" vs "integer" errors from YAML parsing inconsistencies.
  xdim_val <- as.integer(config$xdim)
  ydim_val <- as.integer(config$ydim)
  seed_val <- as.integer(config$seed)
  k_val    <- as.integer(config$n_metaclusters)
  
  # 1. Create FlowFrame for FlowSOM
  exprs_mat <- t(assay(sce, "exprs"))
  
  # 2. Build SOM (Self-Organizing Map)
  set.seed(seed_val)
  
  fsom <- tryCatch({
    FlowSOM::ReadInput(
      input = exprs_mat,
      transform = FALSE, 
      scale = FALSE      
    )
  }, error = function(e) stop("FlowSOM Input Prep failed: ", e$message))
  
  message(sprintf("[Clustering] Building SOM Grid %dx%d...", xdim_val, ydim_val))
  
  fsom <- tryCatch({
    FlowSOM::BuildSOM(
      fsom,
      colsToUse = markers,
      xdim = xdim_val, # Use sanitized integer
      ydim = ydim_val  # Use sanitized integer
    )
  }, error = function(e) stop("FlowSOM BuildSOM failed: ", e$message))
  
  # 3. Build MST
  fsom <- FlowSOM::BuildMST(fsom)
  
  # 4. Metaclustering
  message(sprintf("[Clustering] Consensus Metaclustering into %d populations...", k_val))
  
  meta_clustering <- tryCatch({
    suppressMessages(
      ConsensusClusterPlus::ConsensusClusterPlus(
        t(fsom$map$codes),
        maxK = k_val,
        reps = 100, 
        pItem = 0.9, 
        pFeature = 1, 
        title = tempdir(), 
        plot = "png", 
        verbose = FALSE,
        seed = seed_val
      )
    )
  }, error = function(e) stop("Metaclustering failed: ", e$message))
  
  metacluster_codes <- meta_clustering[[k_val]]$consensusClass
  
  # 5. Map back to individual cells
  cluster_assignments <- fsom$map$mapping[, 1]
  metacluster_assignments <- metacluster_codes[cluster_assignments]
  
  # 6. Append to SCE
  sce$cluster_id <- as.factor(cluster_assignments)
  sce$metacluster_id <- as.factor(metacluster_assignments)
  
  metadata(sce)$FlowSOM <- fsom
  
  return(sce)
}