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
  
  # 1. Create FlowFrame for FlowSOM
  # FlowSOM works natively with flowFrames or matrices.
  exprs_mat <- t(assay(sce, "exprs"))
  
  # 2. Build SOM (Self-Organizing Map)
  # Use set.seed for reproducibility of the starting weights
  set.seed(config$seed)
  
  fsom <- tryCatch({
    FlowSOM::ReadInput(
      input = exprs_mat,
      transform = FALSE, # Data is already transformed/corrected
      scale = FALSE      # cyCombine output is already on comparable scale
    )
  }, error = function(e) stop("FlowSOM Input Prep failed: ", e$message))
  
  message(sprintf("[Clustering] Building SOM Grid %dx%d...", config$xdim, config$ydim))
  
  fsom <- tryCatch({
    FlowSOM::BuildSOM(
      fsom,
      colsToUse = markers,
      xdim = config$xdim,
      ydim = config$ydim
    )
  }, error = function(e) stop("FlowSOM BuildSOM failed: ", e$message))
  
  # 3. Build MST (Minimal Spanning Tree) - Optional but good for viz
  fsom <- FlowSOM::BuildMST(fsom)
  
  # 4. Metaclustering (ConsensusClusterPlus)
  message(sprintf("[Clustering] Consensus Metaclustering into %d populations...", config$n_metaclusters))
  
  meta_clustering <- tryCatch({
    suppressMessages(
      ConsensusClusterPlus::ConsensusClusterPlus(
        t(fsom$map$codes),
        maxK = config$n_metaclusters,
        reps = 100, 
        pItem = 0.9, 
        pFeature = 1, 
        title = tempdir(), 
        plot = "png", 
        verbose = FALSE,
        seed = config$seed
      )
    )
  }, error = function(e) stop("Metaclustering failed: ", e$message))
  
  # Extract specific K
  k <- config$n_metaclusters
  metacluster_codes <- meta_clustering[[k]]$consensusClass
  
  # 5. Map back to individual cells
  # fsom$map$mapping contains the grid assignment (1..100) for each cell
  cluster_assignments <- fsom$map$mapping[, 1]
  metacluster_assignments <- metacluster_codes[cluster_assignments]
  
  # 6. Append to SCE
  sce$cluster_id <- as.factor(cluster_assignments)
  sce$metacluster_id <- as.factor(metacluster_assignments)
  
  # Store the FlowSOM object in metadata for plotting later
  metadata(sce)$FlowSOM <- fsom
  
  return(sce)
}