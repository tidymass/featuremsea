#' @title Internal Helper Functions
#' @description
#' These functions are internal utilities for processing MFM results,
#' handling annotation tables, and calculating scores.
#' They are not intended to be called directly by the user.
#'
#' @name internal_utils
#' @keywords internal
#'
#' @importFrom dplyr mutate filter pull count group_by summarise n
#' @importFrom purrr imap_dfr map_dfr
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom rlang %||%
#' @importFrom stats p.adjust cor
#' @importFrom matrixStats rowMaxs
NULL

# ==============================================================================
# Internal Utility Functions
# ==============================================================================

#' Get significant MFM
#'
#' Filters results based on FDR threshold.
#'
#' @param res_list List of results.
#' @param fdr_threshold Numeric threshold for FDR.
#'
#' @return A data frame of significant MFMs.
#' @noRd
get_significant_mfm <- function(res_list, fdr_threshold = 0.05) {
  res_df <- imap_dfr(
    res_list,
    ~ tibble(
      MFM_id  = .y,
      ES      = .x$ES %||% NA_real_,
      NES     = .x$NES %||% NA_real_,
      p_value = .x$p_value %||% NA_real_
    )
  )
  
  significant_mfm <- res_df %>%
    mutate(
      FDR = {
        pv <- p_value
        adj <- rep(NA_real_, length(pv))
        idx <- is.finite(pv) & !is.na(pv)
        if (any(idx)) adj[idx] <- p.adjust(pv[idx], method = "BH")
        adj
      }
    ) %>%
    filter(!is.na(FDR), FDR < fdr_threshold)
  
  return(significant_mfm)
}

#' Get feature metabolite long table
#'
#' @param significant_mfm Data frame from get_significant_mfm.
#' @param res_list Original result list.
#' @noRd
get_fm_long_table <- function(significant_mfm, res_list, id_col) {
  selected_ids <- significant_mfm %>%
    filter(!is.na(FDR), FDR < 0.05) %>%
    pull(MFM_id)
  
  leading_edge_df <- map_dfr(
    selected_ids,
    function(mid) {
      le <- res_list[[mid]]$leading_edge
      if (!is.null(le) && nrow(le) > 0) {
        le %>% mutate(MFM_id = mid)
      } else {
        NULL
      }
    }
  )
  
  # feature_metabolite_count <- leading_edge_df %>%
  #   dplyr::count(variable_id, KEGG_ID, name = "count")
  
  feature_metabolite_count <- leading_edge_df %>%
    dplyr::count(variable_id, .data[[id_col]], name = "count")
  
  return(feature_metabolite_count)
}

#' Calculate weighted annotation table (Fast Version)
#'
#' Optimized version using matrix operations.
#' 
#' @param original_annotation_table The base annotation matrix/df.
#' @param feature_metabolite_count Data frame with counts.
#' @noRd
get_weighting_annotation_table_fast <- function(original_annotation_table,
                                                feature_metabolite_count) {
  # --- Step 0: Data Preparation ---
  orig_df <- as.data.frame(original_annotation_table)
  orig_mat <- as.matrix(orig_df)
  storage.mode(orig_mat) <- "double"
  
  # --- Step 1: Vectorized construction of the counts matrix ---
  counts_mat <- matrix(0, nrow = nrow(orig_mat), ncol = ncol(orig_mat),
                       dimnames = dimnames(orig_mat))
  
  if (!is.null(feature_metabolite_count) && nrow(feature_metabolite_count) > 0) {
    row_indices <- match(feature_metabolite_count$variable_id, rownames(orig_mat))
    col_indices <- match(feature_metabolite_count$KEGG_ID, colnames(orig_mat))
    
    valid_indices <- !is.na(row_indices) & !is.na(col_indices)
    index_mat <- cbind(row_indices[valid_indices], col_indices[valid_indices])
    
    counts_mat[index_mat] <- feature_metabolite_count$count[valid_indices]
  }
  
  # --- Step 2: High-performance calculation of row statistics ---
  total_count_row <- rowSums(counts_mat, na.rm = TRUE)
  
  # Ensure matrixStats is in Imports
  max_score_row <- matrixStats::rowMaxs(orig_mat, na.rm = TRUE)
  max_score_row[is.infinite(max_score_row)] <- 0
  
  # --- Steps 3, 4, 5: Calculation ---
  denom <- ifelse(total_count_row == 0, 1, total_count_row)
  prop_mat <- sweep(counts_mat, 1, denom, "/")
  
  bonus_mat <- sweep(prop_mat, 1, max_score_row, "*")
  bonus_mat[counts_mat == 0] <- 0
  
  weighted_mat <- orig_mat + bonus_mat
  
  # --- Step 6: Return result ---
  weighted_df <- as.data.frame(weighted_mat, stringsAsFactors = FALSE)
  return(weighted_df)
}

#' Calculate Pearson score between two tables
#' 
#' @noRd
get_pearson_score <- function(annotation_table1, annotation_table2) {
  vec1 <- as.vector(as.matrix(annotation_table1))
  vec2 <- as.vector(as.matrix(annotation_table2))
  
  cor_score <- cor(vec1, vec2, method = "pearson", use = "pairwise.complete.obs")
  return(cor_score)
}

#' Canonicalize counts dataframe
#' Helper for equality check
#' @noRd
canon_counts <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(
      variable_id = character(0),
      KEGG_ID     = character(0),
      count       = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  cols <- intersect(c("variable_id", "KEGG_ID", "count"), colnames(df))
  out <- df[, cols, drop = FALSE]
  if (!"variable_id" %in% colnames(out)) out$variable_id <- NA_character_
  if (!"KEGG_ID"     %in% colnames(out)) out$KEGG_ID     <- NA_character_
  if (!"count"       %in% colnames(out)) out$count       <- NA_real_
  
  out$variable_id <- as.character(out$variable_id)
  out$KEGG_ID     <- as.character(out$KEGG_ID)
  out$count       <- as.numeric(out$count)
  
  out <- out[order(out$variable_id, out$KEGG_ID, out$count), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Check if counts are equal
#' @noRd
counts_equal <- function(a, b) {
  identical(canon_counts(a), canon_counts(b))
}