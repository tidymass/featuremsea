#' Perform Iterative Feature-based Metabolite Set Enrichment Analysis (fMSEA)
#'
#' This function performs an iterative fMSEA analysis. It iteratively updates the
#' weighting of annotations based on the significance of metabolic feature modules (MFMs)
#' until convergence or a maximum number of iterations is reached.
#' 
#' It automatically detects the source of the provided pathway database object.
#' Valid sources are strictly: "KEGG", "SMPDB", "IMETPD", "Reactome", "Wikipathway".
#'
#' @param pathway_database An S4 object containing pathway information. Must contain a 
#'   `database_info` slot with a `source` field indicating the database type.
#' @param annotation_table A data frame containing the initial annotation weights.
#' @param ranking_table A data frame containing the ranking of features (e.g., based on p-values or fold changes).
#' @param threads Integer. Number of threads to use for parallel computing. Default is `NULL` (auto-detect or serial).
#' @param min.compounds.num Integer. Minimum number of compounds required in a pathway. Default is 15.
#' @param max.compounds.num Integer. Maximum number of compounds allowed in a pathway. Default is 300.
#' @param id.col Character. The column name for compound IDs (e.g., "KEGG_ID"). Default is "KEGG_ID".
#' @param perm.num Integer. Number of permutations for significance testing. Default is 1000.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#' @param fdr.thr Numeric. FDR threshold for selecting significant modules. Default is 0.05.
#' @param min.matched.features.num Integer. Post hoc filter: minimum number of features
#'   matched to a significant module for it to be reported. Default is 15. Set to 0 to disable.
#' @param min.leading.edge.num Integer. Post hoc filter: minimum number of features in the
#'   leading edge of a significant module for it to be reported. Default is 2. Set to 0 to disable.
#' @param max.iter.num Integer. Maximum number of iterations for the algorithm. Default is 3.
#' @param verbose Logical. Whether to print progress messages. Default is `TRUE`.
#'
#' @return A \code{featuremsea_object} containing the analysis results. 
#'   The \code{significant_modules} slot will contain columns: 
#'   \code{pathway_id}, \code{pathway_name}, \code{pathway_description}, etc.,
#'   together with \code{n_matched_features} and \code{n_leading_edge_features}.
#'
#' @details The two post hoc filters (\code{min.matched.features.num} and
#'   \code{min.leading.edge.num}) are applied only to the significant modules of the
#'   final iteration, after the iterative analysis is complete. They do not affect the
#'   enrichment calculation, the annotation re-weighting, the permutation test or the
#'   FDR estimation. Modules removed by these filters are recorded in
#'   \code{process_info$post_hoc_filter$removed_modules}.
#'
#' @author Xiaotao Shen \email{xiaotao.shen@@outlook.com}
#' @author Yijiang Liu \email{ejoliu@@outlook.com}
#'
#' @importFrom dplyr left_join select rename %>%
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
    min.matched.features.num = 15,
    min.leading.edge.num = 2,
    max.iter.num = 3,          
    verbose = TRUE
) {
  
  for (arg in c("min.matched.features.num", "min.leading.edge.num")) {
    val <- get(arg)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val < 0) {
      stop(sprintf("%s must be a single non-negative number.", arg))
    }
  }
  
  # --- Step 1: Database preparation ---
  if (!isS4(pathway_database) || !"database_info" %in% slotNames(pathway_database)) {
    stop("pathway_database must be an S4 object with a 'database_info' slot.")
  }
  
  db_source <- pathway_database@database_info$source
  valid_sources <- c("KEGG", "SMPDB", "IMETPD", "Reactome", "Wikipathway")
  
  if (!db_source %in% valid_sources) {
    stop(sprintf("Invalid database source: '%s'", db_source))
  }
  
  pathway_df <- switch(db_source,
                       "SMPDB"       = convert_hmdb2fmsea(pathway_database),
                       "KEGG"        = convert_kegg2fmsea(pathway_database),
                       "IMETPD"      = convert_imetpd2fmsea(pathway_database), 
                       "Reactome"    = convert_reactome2fmsea(pathway_database),
                       "Wikipathway" = convert_wikipathway2fmsea(pathway_database))
  
  pathway_df <- pathway_df %>%
    dplyr::rename(MFM_id = pathway_id, MFM_name = pathway_name, MFM_description = pathway_description)
  
  # --- Step 2: Weight preparation ---
  # Store raw ranking table and create absolute-weight table for calculation
  ranking_table_raw <- ranking_table
  ranking_table_calc <- ranking_table %>%
    dplyr::mutate(ranking_weight = abs(ranking_weight))
  
  # --- Step 3: Iterative enrichment analysis ---
  if (verbose) message("Performing featureMSEA analysis...")
  
  score_annotation_table <- annotation_table
  prev_feature_metabolite_count <- NULL
  converged <- FALSE
  iter_used <- 0L
  last_res_list <- NULL
  last_annotation_table_weighting <- NULL
  last_feature_metabolite_count <- NULL
  last_significant_mfm <- NULL      
  
  for (iter in seq_len(max.iter.num)) {
    iter_used <- iter
    
    # Run core enrichment using absolute weights.
    # Global normalization is applied inside annotation_long_fast_base()
    # to avoid materializing another full-size matrix in memory.
    last_res_list <- parallel_computing_pathways_indexed_fast(
      pathway_dataset = pathway_df,
      annotation_table = score_annotation_table,
      ranking_table = ranking_table_calc,
      min_compounds = min.compounds.num, 
      max_compounds = max.compounds.num, 
      id_col = id.col,                    
      n_perm = perm.num,                  
      seed = seed,
      n_cores = threads,
      verbose = verbose
    )
    
    last_significant_mfm <- get_significant_mfm(last_res_list, fdr_threshold = fdr.thr)
    
    if (nrow(last_significant_mfm) == 0) break
    
    # Update weighting based on significant modules
    last_feature_metabolite_count <- get_fm_long_table(
      last_significant_mfm,
      last_res_list,
      id_col = id.col,
      fdr_threshold = fdr.thr
    )
    last_annotation_table_weighting <- get_weighting_annotation_table_fast(
      annotation_table,
      last_feature_metabolite_count,
      id_col = id.col
    )
    
    # Convergence check
    if (!is.null(prev_feature_metabolite_count) &&
        counts_equal(prev_feature_metabolite_count, last_feature_metabolite_count, id_col = id.col)) {
      converged <- TRUE
      break
    }
    
    prev_feature_metabolite_count <- last_feature_metabolite_count
    score_annotation_table <- last_annotation_table_weighting
  }
  
  # --- Step 4: Post hoc filtering by feature numbers ---
  # Applied only to the final significant modules; the iterations above,
  # including re-weighting and FDR, are not affected.
  removed_mfm <- data.frame(MFM_id = character())
  if (!is.null(last_significant_mfm) && nrow(last_significant_mfm) > 0) {
    filtered <- filter_mfm_by_feature_number(
      last_significant_mfm,
      last_res_list,
      min_matched = min.matched.features.num,
      min_leading_edge = min.leading.edge.num
    )
    last_significant_mfm <- filtered$kept
    removed_mfm <- filtered$removed
    if (verbose && nrow(removed_mfm) > 0) {
      message(sprintf(
        "%d significant module(s) removed by the post hoc filters (matched features < %s or leading-edge features < %s): %s",
        nrow(removed_mfm), min.matched.features.num, min.leading.edge.num,
        paste(removed_mfm$MFM_id, collapse = ", ")
      ))
    }
  }
  
  # --- Step 5: Formatting and Output ---
  if (is.null(last_feature_metabolite_count)) last_feature_metabolite_count <- data.frame() 
  if (is.null(last_annotation_table_weighting)) last_annotation_table_weighting <- data.frame()
  if (is.null(last_significant_mfm)) last_significant_mfm <- data.frame(MFM_id = character())
  
  pathway_cols <- c("MFM_id", "MFM_name", "MFM_description", "pathway_class_all")
  if (id.col %in% colnames(pathway_df)) {
    pathway_cols <- c(pathway_cols, id.col)
  }
  
  significant_modules_final <- last_significant_mfm %>%
    dplyr::left_join(pathway_df[, pathway_cols], by = "MFM_id") %>%
    dplyr::rename(pathway_id = MFM_id, pathway_name = MFM_name, pathway_description = MFM_description) %>%
    dplyr::select(dplyr::all_of(pathway_cols[pathway_cols != "MFM_id" & pathway_cols != "MFM_name" & pathway_cols != "MFM_description"]), ES, NES, p_value, FDR, dplyr::any_of(c("n_matched_features", "n_leading_edge_features")), pathway_id, pathway_name, pathway_description) %>%
    dplyr::select(pathway_id, pathway_name, pathway_description, dplyr::everything())
  
  result_object <- new(
    Class = "featuremsea_object",
    ranking_table              = ranking_table_raw,
    feature_metabolite_count   = last_feature_metabolite_count,
    annotation_table_weighting = last_annotation_table_weighting,
    significant_modules        = significant_modules_final,   
    res_list                   = last_res_list,
    converged                  = converged,
    iterations_used            = iter_used,
    process_info               = list(
      creation_date = as.character(Sys.time()),
      post_hoc_filter = list(
        min.matched.features.num = min.matched.features.num,
        min.leading.edge.num = min.leading.edge.num,
        removed_modules = removed_mfm
      )
    )
  )
  
  return(result_object)
}
