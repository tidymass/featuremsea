#' Perform Iterative Feature-based Metabolite Set Enrichment Analysis (fMSEA)
#'
#' This function performs an iterative fMSEA analysis. It iteratively updates the
#' weighting of annotations based on the significance of metabolic feature modules (MFMs)
#' until convergence or a maximum number of iterations is reached.
#' 
#' @description
#' The function supports automatic column detection for multiple databases.
#' Internally it normalizes columns to "MFM_id", but the final output will 
#' restore the original column names (e.g., "KEGG_id").
#'
#' @param pathway_database A data frame containing pathway information. 
#'   Supports standard columns (`MFM_id`, `MFM_name`, `MFM_description`) 
#'   or database-specific columns (e.g., `KEGG_id`, `KEGG_name`, `KEGG_description`).
#'   Must also contain `pathway_class_all`.
#' @param annotation_table A data frame containing the initial annotation weights.
#' @param ranking_table A data frame containing the ranking of features.
#' @param threads Integer. Number of threads to use for parallel computing.
#' @param min.compounds.num Integer. Minimum number of compounds required in a pathway.
#' @param max.compounds.num Integer. Maximum number of compounds allowed in a pathway.
#' @param id.col Character. The column name for compound IDs (e.g., "KEGG_ID").
#' @param perm.num Integer. Number of permutations for significance testing.
#' @param seed Integer. Random seed for reproducibility.
#' @param fdr.thr Numeric. FDR threshold for selecting significant modules.
#' @param max.iter.num Integer. Maximum number of iterations for the algorithm.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A \code{featuremsea_object} containing the analysis results.
#'
#' @author Xiaotao Shen \email{xiaotao.shen@@outlook.com}
#' @author Yijiang Liu \email{ejoliu@@outlook.com}
#'
#' @importFrom dplyr left_join select rename %>%
#' @importFrom rlang :=
#' @export
perform_fmsea_analysis_new <- function(
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
  
  # ==============================================================================
  # 1. Database Column Detection & Normalization
  # ==============================================================================
  
  # Define mappings for supported databases
  db_mappings <- list(
    "KEGG"        = c(id = "KEGG_id",        name = "KEGG_name",        desc = "KEGG_description"),
    "HMDB"        = c(id = "HMDB_id",        name = "HMDB_name",        desc = "HMDB_description"),
    "REACTOME"    = c(id = "REACTOME_id",    name = "REACTOME_name",    desc = "REACTOME_description"),
    "WIKIPATHWAY" = c(id = "WIKIPATHWAY_id", name = "WIKIPATHWAY_name", desc = "WIKIPATHWAY_description")
  )
  
  input_cols <- colnames(pathway_database)
  is_standardized <- FALSE
  
  # Variables to store the FINAL output column names (Default to MFM)
  # We will update these if we detect a specific database like KEGG
  final_col_names <- list(id = "MFM_id", name = "MFM_name", desc = "MFM_description")
  
  # Check 1: Is it already standard/IMETPD format?
  if (all(c("MFM_id", "MFM_name", "MFM_description") %in% input_cols)) {
    is_standardized <- TRUE
    if (verbose) message("  - Detected standard/IMETPD database format.")
  } else {
    # Check 2: Loop through other DBs
    for (db_name in names(db_mappings)) {
      map <- db_mappings[[db_name]]
      
      if (all(map %in% input_cols)) {
        if (verbose) message(sprintf("Detected %s database format. Normalizing internally...", db_name))
        
        # KEY STEP: Save the original names for the final output
        final_col_names <- as.list(map)
        
        # Rename to internal standard (MFM_*) for calculation
        # pathway_database <- pathway_database %>%
        #   dplyr::rename(
        #     MFM_id = !!map["id"],
        #     MFM_name = !!map["name"],
        #     MFM_description = !!map["desc"]
        #   )
        
        pathway_database <- pathway_database %>%
          dplyr::rename(
            MFM_id = !!sym(map["id"]),
            MFM_name = !!sym(map["name"]),
            MFM_description = !!sym(map["desc"])
          )
        
        is_standardized <- TRUE
        break
      }
    }
  }
  
  # Check 3: Check for pathway_class_all
  if (!"pathway_class_all" %in% colnames(pathway_database)) {
    if (verbose) message("  - Warning: 'pathway_class_all' column not found. Creating NA column.")
    pathway_database$pathway_class_all <- NA
  }
  
  if (!is_standardized) {
    stop("Error: The 'pathway_database' columns could not be identified. Please ensure columns match one of the supported formats.")
  }
  
  # ==============================================================================
  # 2. Iterative Calculation (Uses MFM_* internally)
  # ==============================================================================
  
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
    if (verbose) message(sprintf("Iteration %d: computing enrichment...", iter))
    
    # Normalize
    annotation_table_norm <- normalize_score_annotation_global(score_annotation = score_annotation_table)
    
    # Compute
    last_res_list <- parallel_computing_pathways_indexed_fast(
      pathway_dataset = pathway_database, # This now has MFM columns
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
    
    # Select significant
    last_significant_mfm <- get_significant_mfm(
      last_res_list,
      fdr_threshold = fdr.thr 
    )
    
    if (nrow(last_significant_mfm) == 0) {
      if (verbose) message("No significant modules found. Stopping.")
      break
    }
    
    # Count pairs
    last_feature_metabolite_count <- get_fm_long_table(
      last_significant_mfm,
      last_res_list,
      id_col = id.col 
    )
    
    # New weights
    last_annotation_table_weighting <- get_weighting_annotation_table_fast(
      annotation_table,
      last_feature_metabolite_count
    )
    
    # Convergence
    if (!is.null(prev_feature_metabolite_count) &&
        counts_equal(prev_feature_metabolite_count, last_feature_metabolite_count)) {
      converged <- TRUE
      if (verbose) message(sprintf("Converged at iteration %d.", iter))
      break
    }
    
    prev_feature_metabolite_count <- last_feature_metabolite_count
    score_annotation_table <- last_annotation_table_weighting
  }
  
  if (!converged && verbose) {
    message(sprintf("Reached max iteration (%d) without convergence.", max.iter.num))
  }
  
  # Handle NULLs
  if (is.null(last_feature_metabolite_count)) last_feature_metabolite_count <- data.frame() 
  if (is.null(last_annotation_table_weighting)) last_annotation_table_weighting <- data.frame()
  if (is.null(last_significant_mfm)) last_significant_mfm <- data.frame(MFM_id = character())
  
  # ==============================================================================
  # 3. Result Construction & Renaming Back (Output Restoration)
  # ==============================================================================
  
  # Step A: Join using the internal MFM_id (since pathway_database was temporarily renamed)
  significant_modules_temp <- last_significant_mfm %>%
    dplyr::left_join(
      pathway_database[, c("MFM_id", "MFM_name", "MFM_description", "pathway_class_all")], 
      by = c("MFM_id" = "MFM_id")
    ) %>%
    dplyr::select(
      MFM_id, MFM_name, MFM_description, pathway_class_all,
      ES, NES, p_value, FDR
    )
  
  # Step B: Rename columns back to the original database format (e.g. KEGG_id)
  # using the 'final_col_names' we saved at the beginning
  significant_modules_final <- significant_modules_temp %>%
    dplyr::rename(
      !!final_col_names$id   := MFM_id,
      !!final_col_names$name := MFM_name,
      !!final_col_names$desc := MFM_description
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