#' Parallel mMSEA computation for multiple pathways (indexed + permutation)
#'
#' This function performs pathway-level mMSEA analysis in parallel using
#' the indexed fast implementation and C++-powered permutations.
#'
#' @param pathway_dataset A data frame containing pathway information.
#'   Must include a pathway identifier column \code{MFM_id} and at least one ID
#'   column among \code{KEGG_ID} or \code{HMDB_ID}. The \code{id_col} argument
#'   selects which one to use.
#'
#' @param annotation_table A feature-by-compound annotation matrix.
#'
#' @param ranking_table A data frame containing \code{variable_id} and
#'   \code{ranking_weight} columns for ranking features.
#'
#' @param min_compounds Minimum number of compounds required for a pathway
#'   to be retained.
#'
#' @param max_compounds Maximum number of compounds allowed for a pathway
#'   to be retained.
#'
#' @param id_col Name of the column in \code{pathway_dataset} containing
#'   compound IDs (in "{}" concatenated format).
#'
#' @param n_perm Number of permutations for each pathway.
#'
#' @param seed Random seed for reproducibility.
#'
#' @param return_perm Logical; whether to return permutation ES values.
#'
#' @param n_cores Number of CPU cores to use. If \code{NULL}, an automatic
#'   selection is applied.
#'
#' @param verbose Logical; print progress messages.
#'
#' @return
#' A named list of mMSEA results for each pathway. Each element contains:
#' \itemize{
#'   \item \code{ES} – observed enrichment score
#'   \item \code{p_value} – permutation p-value
#'   \item \code{NES} – normalized enrichment score
#'   \item \code{leading_edge} – leading-edge subset
#'   \item \code{steps} – step-by-step enrichment curve
#'   \item \code{perm_ES} – (optional) permutation ES values
#' }
#'
#' @importFrom foreach %dopar%
#' @importFrom dplyr %>%
#' @export
#'
parallel_computing_pathways_indexed_fast <- function(
    pathway_dataset,
    annotation_table,
    ranking_table,
    min_compounds = 5,
    max_compounds = 300,
    id_col = "KEGG_ID",
    n_perm = 1000,
    seed = 123,
    return_perm = FALSE,
    n_cores = NULL,
    verbose = TRUE
){
  # ---- Basic input checks ----
  if (!all(c("variable_id", "ranking_weight") %in% colnames(ranking_table))) {
    stop("ranking_table must contain 'variable_id' and 'ranking_weight' columns.")
  }
  if (!any(c("KEGG_ID", "HMDB_ID") %in% colnames(pathway_dataset))) {
    stop("pathway_dataset must contain at least one of KEGG_ID or HMDB_ID.")
  }
  if (!"MFM_id" %in% colnames(pathway_dataset)) {
    stop("pathway_dataset must contain the MFM_id column.")
  }
  
  # ---- Filter pathways ----
  pathway_dataset <- pathway_dataset %>%
    dplyr::mutate(
      id_count = ifelse(
        is.na(.data[[id_col]]), 0L,
        stringr::str_count(.data[[id_col]], "\\{\\}") + 1L
      )
    ) %>%
    dplyr::filter(
      id_count >= min_compounds & id_count <= max_compounds
    )
  
  if (nrow(pathway_dataset) == 0L) {
    if (verbose) message("No pathways left after filtering by compound count.")
    return(list())
  }
  
  # ---- Preprocessing shared across all pathways ----
  annotation_long <- annotation_long_fast_base(annotation_table, id_col = id_col)
  
  base_rk <- ranking_table %>%
    dplyr::distinct(variable_id, ranking_weight) %>%
    dplyr::arrange(dplyr::desc(ranking_weight))
  
  rk_idx <- build_ranking_index(base_rk)
  
  # ---- Parallel configuration ----
  if (is.null(n_cores)) {
    n_cores <- max(1L, min(parallel::detectCores() - 1L, 4L))
  }
  
  if (verbose) {
    cat("Total pathways to analyze:", nrow(pathway_dataset), "\n")
    cat("Using", n_cores, "CPU cores\n")
    cat("Starting parallel computation...\n")
  }
  
  doParallel::registerDoParallel(cores = n_cores)
  
  # ---- Main parallel mMSEA loop ----
  res_list <- foreach::foreach(
    i = seq_len(nrow(pathway_dataset)),
    .packages = c("dplyr", "stringr"),
    .export   = c(
      "split_ids",
      "annotation_long_fast_base",
      "build_ranking_index",
      "precompute_mMSEA_static",
      "compute_ES_from_mapping",
      "get_mMSEA_results_indexed_fast",
      "compute_permutations_cpp"
    )
  ) %dopar% {
    pathway_ids_vec <- split_ids(pathway_dataset[[id_col]][[i]])
    
    get_mMSEA_results_indexed_fast(
      pathway_ids_vec = pathway_ids_vec,
      annotation_long = annotation_long,
      rk_idx          = rk_idx,
      id_col          = id_col,
      n_perm          = n_perm,
      seed            = seed + i,
      return_perm     = return_perm
    )
  }
  
  doParallel::stopImplicitCluster()
  
  names(res_list) <- pathway_dataset$MFM_id
  
  if (verbose) cat("Parallel analysis completed.\n")
  
  res_list
}
