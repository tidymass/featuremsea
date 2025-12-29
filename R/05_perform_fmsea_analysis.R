#' Perform Iterative Feature-based Metabolite Set Enrichment Analysis (fMSEA)
#'
#' This function performs an iterative fMSEA analysis. It iteratively updates the
#' weighting of annotations based on the significance of metabolic feature modules (MFMs)
#' until convergence or a maximum number of iterations is reached.
#'
#' @param pathway_database A data frame containing pathway information. Must include columns
#'   `MFM_id`, `MFM_name`, `MFM_description`, and `pathway_class_all`.
#' @param annotation_table A data frame containing the initial annotation weights.
#' @param ranking_table A data frame containing the ranking of features (e.g., based on p-values or fold changes).
#' @param threads Integer. Number of threads to use for parallel computing. Default is `NULL` (auto-detect or serial).
#' @param min.compounds.num Integer. Minimum number of compounds required in a pathway. Default is 15.
#' @param max.compounds.num Integer. Maximum number of compounds allowed in a pathway. Default is 300.
#' @param id.col Character. The column name for compound IDs (e.g., "KEGG_ID"). Default is "KEGG_ID".
#' @param perm.num Integer. Number of permutations for significance testing. Default is 1000.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#' @param fdr.thr Numeric. FDR threshold for selecting significant modules. Default is 0.05.
#' @param max.iter.num Integer. Maximum number of iterations for the algorithm. Default is 10.
#' @param verbose Logical. Whether to print progress messages. Default is `TRUE`.
#'
#' @return A list containing:
#' \describe{
#'   \item{feature_metabolite_count}{The final count of feature-metabolite pairs.}
#'   \item{annotation_table_weighting}{The final weighted annotation table.}
#'   \item{significant_modules}{A data frame of significant MFMs with enrichment statistics (ES, NES, p-value, FDR).}
#'   \item{res_list}{The raw results list from the last iteration.}
#'   \item{converged}{Logical. Whether the algorithm converged.}
#'   \item{iterations_used}{Integer. The number of iterations performed.}
#' }
#'
#' @author Xiaotao Shen \email{xiaotao.shen@@outlook.com}
#' @author Yijiang Liu \email{ejoliu@@outlook.com}
#'
#' @importFrom dplyr left_join select %>%
#' @export
perform_fmsea_analysis <- function(
    pathway_database,            
    annotation_table,
    ranking_table,
    threads = NULL,
    min.compounds.num = 15,   
    max.compounds.num = 300,  
    id.col = "KEGG_ID",       
    perm.num = 1000,          
    seed = 123,
    fdr.thr = 0.05,           
    max.iter.num = 10,        
    verbose = TRUE
) {
  if (verbose) message("Starting FMSEA iterative analysis...")
  
  score_annotation_table <- annotation_table
  prev_feature_metabolite_count <- NULL
  
  converged <- FALSE
  iter_used <- 0L
  last_res_list <- NULL
  last_annotation_table_weighting <- NULL
  last_feature_metabolite_count <- NULL
  last_significant_mfm <- NULL    
  
  # Use max.iter.num here
  for (iter in seq_len(max.iter.num)) {
    iter_used <- iter
    if (verbose) message(sprintf("Iteration %d: computing enrichment with current score annotation...", iter))
    
    # 1) normalize score annotation
    annotation_table_norm <- normalize_score_annotation_new(score_annotation = score_annotation_table)
    
    # 2) compute results
    last_res_list <- parallel_computing_pathways_indexed_fast(
      pathway_dataset = pathway_database, 
      annotation_table = annotation_table_norm,
      ranking_table = ranking_table,
      min_compounds = min.compounds.num, # Map new param to internal param
      max_compounds = max.compounds.num, # Map new param to internal param
      id_col = id.col,                   # Map new param to internal param
      n_perm = perm.num,                 # Map new param to internal param
      seed = seed,
      return_perm = FALSE,
      n_cores = threads,
      verbose = verbose
    )
    
    # 3) select significant modules
    last_significant_mfm <- get_significant_mfm(
      last_res_list,
      fdr_threshold = fdr.thr # Map new param to internal param
    )
    
    if (nrow(last_significant_mfm) == 0) {
      if (verbose) message(sprintf("No significant modules found under the current FDR threshold. Stopping iteration."))
      break
    }
    
    # 4) count feature–metabolite pairs
    last_feature_metabolite_count <- get_fm_long_table(
      last_significant_mfm,
      last_res_list,
      id_col = id.col # Map new param to internal param
    )
    
    # 5) new weighting table
    last_annotation_table_weighting <- get_weighting_annotation_table_fast(
      annotation_table,
      last_feature_metabolite_count
    )
    
    # 6) convergence check
    if (!is.null(prev_feature_metabolite_count) &&
        counts_equal(prev_feature_metabolite_count, last_feature_metabolite_count)) {
      converged <- TRUE
      if (verbose) message(sprintf("Converged at iteration %d: feature_metabolite_count unchanged.", iter))
      break
    }
    
    prev_feature_metabolite_count <- last_feature_metabolite_count
    score_annotation_table <- last_annotation_table_weighting
    
    if (verbose) message(sprintf("Iteration %d complete. Proceeding to the next iteration.", iter))
  }
  
  if (!converged && verbose) {
    message(sprintf("Reached the maximum iteration limit (%d) without convergence.", max.iter.num))
  }
  
  significant_modules <- last_significant_mfm %>%
    dplyr::left_join(
      pathway_database[, c("MFM_id", "MFM_name", "MFM_description", "pathway_class_all")], 
      by = c("MFM_id" = "MFM_id")
    ) %>%
    dplyr::select(
      MFM_id, MFM_name, MFM_description, pathway_class_all,
      ES, NES, p_value, FDR
    )
  
  list(
    feature_metabolite_count   = last_feature_metabolite_count,
    annotation_table_weighting = last_annotation_table_weighting,
    significant_modules        = significant_modules,   
    res_list                   = last_res_list,
    converged                  = converged,
    iterations_used            = iter_used
  )
}