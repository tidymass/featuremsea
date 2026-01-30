#' Process Metabolite Annotation Table
#'
#' Reshapes the annotation table into a score matrix and extracts ranking weights.
#'
#' @param annotation_table_final2 Data frame containing annotation results.
#' @param database_type Character string, either "HMDB" or "KEGG".
#'
#' @return A list containing:
#'   \item{original_score_annotation}{Data frame of score matrix (variables rows x metabolites columns).}
#'   \item{ranking_table}{Data frame with variable_id and ranking_weight columns.}
#'
#' @importFrom data.table as.data.table setnames
#' @importFrom Matrix sparseMatrix
#' @importFrom dplyr rename filter distinct arrange desc all_of
#' @importFrom rlang sym .data
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#' @export
process_annotation_table <- function(annotation_table_final2, 
                                     database_type = "HMDB") {
  
  # Validate database_type
  if (!database_type %in% c("HMDB", "KEGG")) {
    stop("database_type must be either 'HMDB' or 'KEGG'")
  }
  
  # Set ID column name based on database type
  id_col_original <- if (database_type == "HMDB") "HMDB.ID" else "KEGG.ID"
  id_col_new <- if (database_type == "HMDB") "HMDB_ID" else "KEGG_ID"
  
  # Check if original column exists before proceeding
  if (!id_col_original %in% names(annotation_table_final2)) {
    stop(paste0("Column '", id_col_original, "' not found in input data."))
  }
  
  # ========== Part 1: Create original_score_annotation (using data.table) ==========
  # Ensure we don't modify the original input by reference
  dt <- data.table::as.data.table(annotation_table_final2)
  
  # Rename ID column
  data.table::setnames(dt, old = id_col_original, new = id_col_new)
  
  # Remove rows with missing IDs
  # Note: usage of get() is correct for dynamic column names in data.table
  dt <- dt[!is.na(get(id_col_new))]
  
  # Check if data is empty after filtering
  if (nrow(dt) == 0) {
    warning("No valid data after removing NA values.")
    return(list(
      original_score_annotation = data.frame(),
      ranking_table = data.frame()
    ))
  }
  
  # Take maximum score for each (variable_id, metabolite_ID) pair
  # .SD, .BY are not needed here, standard aggregation is fine
  dt_max <- dt[, .(score = max(score, na.rm = TRUE)), 
               by = c("variable_id", id_col_new)]
  
  # Construct sparse matrix
  row_fac <- factor(dt_max$variable_id)
  col_fac <- factor(dt_max[[id_col_new]])
  
  ann_sparse <- Matrix::sparseMatrix(
    i        = as.integer(row_fac),
    j        = as.integer(col_fac),
    x        = dt_max$score,
    dims     = c(nlevels(row_fac), nlevels(col_fac)),
    dimnames = list(levels(row_fac), levels(col_fac))
  )
  
  # Convert to regular matrix then data frame
  ann_mat <- as.matrix(ann_sparse)
  original_score_annotation <- as.data.frame(ann_mat)
  
  # ========== Part 2: Create ranking_table (using dplyr) ==========
  # We start from original input again to be safe
  
  # Define dynamic renaming
  rename_list <- stats::setNames(id_col_original, id_col_new)
  
  ranking_table <- annotation_table_final2 %>%
    dplyr::rename(!!rlang::sym(id_col_new) := dplyr::all_of(id_col_original)) %>%
    dplyr::filter(!is.na(.data[[id_col_new]])) %>%
    # Assuming 'condition' exists in the input data
    dplyr::distinct(variable_id, ranking_weight = abs(condition)) %>%
    dplyr::arrange(dplyr::desc(ranking_weight))
  
  # ========== Return results ==========
  return(list(
    original_score_annotation = original_score_annotation,
    ranking_table = ranking_table
  ))
}





