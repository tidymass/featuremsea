#' Perform Iterative Feature-based Metabolite Set Enrichment Analysis (fMSEA)
#'
#' This function performs an iterative fMSEA analysis. It iteratively updates the
#' weighting of annotations based on the significance of metabolic feature modules (MFMs)
#' until convergence or a maximum number of iterations is reached.
#'
#' @param pathway_dataset A data frame containing pathway information. Must include columns
#'   `MFM_id`, `MFM_name`, `MFM_description`, and `pathway_class_all`.
#' @param annotation_table A data frame containing the initial annotation weights.
#' @param ranking_table A data frame containing the ranking of features (e.g., based on p-values or fold changes).
#' @param n_cores Integer. Number of cores to use for parallel computing. Default is `NULL` (auto-detect or serial).
#' @param min_compounds Integer. Minimum number of compounds required in a pathway. Default is 15.
#' @param max_compounds Integer. Maximum number of compounds allowed in a pathway. Default is 300.
#' @param id_col Character. The column name for compound IDs (e.g., "KEGG_ID"). Default is "KEGG_ID".
#' @param n_perm Integer. Number of permutations for significance testing. Default is 1000.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#' @param fdr_threshold Numeric. FDR threshold for selecting significant modules. Default is 0.05.
#' @param max_iter Integer. Maximum number of iterations for the algorithm. Default is 10.
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
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#' @author Yijiang Liu \email{liuyj8066@gmail.com}
#'
#' @importFrom dplyr left_join select %>%
#' @export
perform_fmsea_analysis <- function(
    pathway_dataset,
    annotation_table,
    ranking_table,
    n_cores = NULL,
    min_compounds = 15,
    max_compounds = 300,
    id_col = "KEGG_ID",
    n_perm = 1000,
    seed = 123,
    fdr_threshold = 0.05,
    max_iter = 10,              # maximum iterations
    verbose = TRUE              # print progress
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
  
  for (iter in seq_len(max_iter)) {
    iter_used <- iter
    if (verbose) message(sprintf("Iteration %d: computing enrichment with current score annotation...", iter))
    
    # 1) normalize score annotation (use the NEW function)
    # Ensure normalize_score_annotation_new is defined in your package
    annotation_table_norm <- normalize_score_annotation_new(score_annotation = score_annotation_table)
    
    # 2) compute results (use the NEW fast indexed parallel function)
    # Ensure parallel_computing_pathways_indexed_fast is defined in your package
    last_res_list <- parallel_computing_pathways_indexed_fast(
      pathway_dataset = pathway_dataset,
      annotation_table = annotation_table_norm,
      ranking_table = ranking_table,
      min_compounds = min_compounds,
      max_compounds = max_compounds,
      id_col = id_col,
      n_perm = n_perm,
      seed = seed,
      n_cores = n_cores,
      verbose = verbose
    )
    
    # 3) select significant modules
    # Ensure get_significant_mfm is defined in your package
    last_significant_mfm <- get_significant_mfm(
      last_res_list,
      fdr_threshold = fdr_threshold
    )
    
    # 4) count feature–metabolite pairs
    # Ensure get_fm_long_table is defined in your package
    last_feature_metabolite_count <- get_fm_long_table(
      last_significant_mfm,
      last_res_list
    )
    
    # 5) new weighting table
    # Ensure get_weighting_annotation_table_fast is defined in your package
    last_annotation_table_weighting <- get_weighting_annotation_table_fast(
      annotation_table,
      last_feature_metabolite_count
    )
    
    # 6) convergence check
    # Ensure counts_equal is defined in your package
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
    message(sprintf("Reached the maximum iteration limit (%d) without convergence.", max_iter))
  }
  
  significant_modules <- last_significant_mfm %>%
    dplyr::left_join(
      pathway_dataset[, c("MFM_id", "MFM_name", "MFM_description", "pathway_class_all")],
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