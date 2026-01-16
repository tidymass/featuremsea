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
#' @param max.iter.num Integer. Maximum number of iterations for the algorithm. Default is 10.
#' @param verbose Logical. Whether to print progress messages. Default is `TRUE`.
#'
#' @return A \code{featuremsea_object} containing the analysis results. 
#'   The \code{significant_modules} slot will contain columns: 
#'   \code{pathway_id}, \code{pathway_name}, \code{pathway_description}, etc.
#'
#' @author Xiaotao Shen \email{xiaotao.shen@@outlook.com}
#' @author Yijiang Liu \email{ejoliu@@outlook.com}
#'
#' @importFrom dplyr left_join select rename %>%
#' @export
perform_fmsea_analysis_new2 <- function(
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
  
  # --- STEP 0: Database Conversion ---
  if (verbose) message("Checking pathway database source...")
  
  # Ensure the object has the required slot
  if (!isS4(pathway_database) || !"database_info" %in% slotNames(pathway_database)) {
    stop("pathway_database must be an S4 object with a 'database_info' slot.")
  }
  
  db_source <- pathway_database@database_info$source
  
  if (is.null(db_source)) {
    stop("Could not determine database source. pathway_database@database_info$source is NULL.")
  }
  
  # Strict validation of source names (Updated HMDB -> SMPDB)
  valid_sources <- c("KEGG", "SMPDB", "IMETPD", "Reactome", "Wikipathway")
  
  if (!db_source %in% valid_sources) {
    stop(sprintf(
      "Invalid database source: '%s'.\nAllowed sources are: %s", 
      db_source, 
      paste(valid_sources, collapse = ", ")
    ))
  }
  
  # Select the correct conversion function
  # Note: SMPDB source uses convert_hmdb2fmsea
  pathway_df <- switch(db_source,
                       "SMPDB"       = convert_hmdb2fmsea(pathway_database),
                       "KEGG"        = convert_kegg2fmsea(pathway_database),
                       "IMETPD"      = convert_imetpd2fmsea(pathway_database), 
                       "Reactome"    = convert_reactome2fmsea(pathway_database),
                       "Wikipathway" = convert_wikipathway2fmsea(pathway_database)
  )
  
  # --- 修改开始：处理显示名称 ---
  # 如果内部识别为 SMPDB，则显示名称改为 HMDB；否则保持原样
  display_source_name <- if (db_source == "SMPDB") "HMDB" else db_source
  
  if (verbose) message(sprintf("Successfully converted database from source: %s", display_source_name))
  # --- 修改结束 ---
  
  # INTERNAL RENAMING: 
  # Rename columns to MFM_* for internal processing compatibility
  pathway_df <- pathway_df %>%
    dplyr::rename(
      MFM_id = pathway_id,
      MFM_name = pathway_name,
      MFM_description = pathway_description
    )
  
  # --- END STEP 0 ---
  
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
    annotation_table_norm <- normalize_score_annotation_global(score_annotation = score_annotation_table)
    
    # 2) compute results
    last_res_list <- parallel_computing_pathways_indexed_fast(
      pathway_dataset = pathway_df,   # Use the converted internal data frame
      annotation_table = annotation_table_norm,
      ranking_table = ranking_table,
      min_compounds = min.compounds.num, 
      max_compounds = max.compounds.num, 
      id_col = id.col,                    
      n_perm = perm.num,                  
      seed = seed,
      return_perm = FALSE,
      n_cores = threads,
      verbose = verbose
    )
    
    # 3) select significant modules
    last_significant_mfm <- get_significant_mfm(
      last_res_list,
      fdr_threshold = fdr.thr 
    )
    
    if (nrow(last_significant_mfm) == 0) {
      if (verbose) message(sprintf("No significant modules found under the current permutation number and FDR threshold. Stopping iteration."))
      break
    }
    
    # 4) count feature–metabolite pairs
    last_feature_metabolite_count <- get_fm_long_table(
      last_significant_mfm,
      last_res_list,
      id_col = id.col 
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
  
  # --- Handle NULLs if loop broke early (No significant modules) ---
  if (is.null(last_feature_metabolite_count)) {
    last_feature_metabolite_count <- data.frame() 
  }
  
  if (is.null(last_annotation_table_weighting)) {
    last_annotation_table_weighting <- data.frame()
  }
  
  if (is.null(last_significant_mfm)) {
    last_significant_mfm <- data.frame(MFM_id = character())
  }
  
  # --- FINAL OUTPUT FORMATTING ---
  # Join with pathway info and RENAME back to pathway_*
  significant_modules_final <- last_significant_mfm %>%
    dplyr::left_join(
      pathway_df[, c("MFM_id", "MFM_name", "MFM_description", "pathway_class_all")], 
      by = c("MFM_id" = "MFM_id")
    ) %>%
    dplyr::rename(
      pathway_id = MFM_id,
      pathway_name = MFM_name,
      pathway_description = MFM_description
    ) %>%
    dplyr::select(
      pathway_id, pathway_name, pathway_description, pathway_class_all,
      ES, NES, p_value, FDR
    )
  
  # create and return featuremsea_object
  result_object <- new(
    Class = "featuremsea_object",
    feature_metabolite_count   = last_feature_metabolite_count,
    annotation_table_weighting = last_annotation_table_weighting,
    significant_modules        = significant_modules_final,   
    res_list                   = last_res_list,
    converged                  = converged,
    iterations_used            = iter_used,
    process_info               = list(creation_date = as.character(Sys.time()))
  )
  
  return(result_object)
}