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
get_fm_long_table <- function(significant_mfm, res_list, id_col, fdr_threshold = 0.05) {
  selected_ids <- significant_mfm %>%
    filter(!is.na(FDR), FDR < fdr_threshold) %>%
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
                                                feature_metabolite_count,
                                                id_col = "KEGG_ID") {
  # --- Step 0: Data Preparation ---
  weighted_df <- as.data.frame(original_annotation_table, stringsAsFactors = FALSE)
  n_row <- nrow(weighted_df)
  n_col <- ncol(weighted_df)
  
  # --- Step 1: Row-wise max score (streaming over columns, no dense copy) ---
  max_score_row <- rep(-Inf, n_row)
  for (j in seq_len(n_col)) {
    col_j <- weighted_df[[j]]
    if (!is.numeric(col_j)) {
      col_j <- suppressWarnings(as.numeric(col_j))
    }
    col_for_max <- col_j
    col_for_max[is.na(col_for_max)] <- -Inf
    max_score_row <- pmax(max_score_row, col_for_max)
    weighted_df[[j]] <- col_j
  }
  max_score_row[!is.finite(max_score_row)] <- 0
  
  # Early return if no updates are needed
  if (is.null(feature_metabolite_count) || nrow(feature_metabolite_count) == 0) {
    return(weighted_df)
  }
  
  if (!id_col %in% colnames(feature_metabolite_count)) {
    stop("feature_metabolite_count must contain id column: ", id_col)
  }
  
  # --- Step 2: Sparse index matching ---
  row_indices <- match(feature_metabolite_count$variable_id, rownames(weighted_df))
  col_indices <- match(feature_metabolite_count[[id_col]], colnames(weighted_df))
  counts_vec <- as.numeric(feature_metabolite_count$count)
  
  valid <- !is.na(row_indices) & !is.na(col_indices) & is.finite(counts_vec) & counts_vec > 0
  if (!any(valid)) {
    return(weighted_df)
  }
  
  row_indices <- row_indices[valid]
  col_indices <- col_indices[valid]
  counts_vec <- counts_vec[valid]
  
  # --- Step 3: Row-wise total count for normalization ---
  total_count_row <- numeric(n_row)
  row_sum_tbl <- rowsum(counts_vec, group = row_indices, reorder = FALSE)
  total_count_row[as.integer(rownames(row_sum_tbl))] <- as.numeric(row_sum_tbl[, 1])
  denom <- ifelse(total_count_row == 0, 1, total_count_row)
  
  # --- Step 4: Compute sparse bonuses for hit cells only ---
  bonus_vec <- counts_vec / denom[row_indices] * max_score_row[row_indices]
  
  # --- Step 5: Apply sparse updates grouped by column ---
  split_idx <- split(seq_along(col_indices), col_indices)
  for (col_key in names(split_idx)) {
    idx <- split_idx[[col_key]]
    cj <- as.integer(col_key)
    rows <- row_indices[idx]
    addv <- bonus_vec[idx]
    col_j <- weighted_df[[cj]]
    col_j[rows] <- col_j[rows] + addv
    weighted_df[[cj]] <- col_j
  }
  
  weighted_df
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
canon_counts <- function(df, id_col = "KEGG_ID") {
  if (is.null(df) || nrow(df) == 0) {
    out <- data.frame(
      variable_id = character(0),
      count       = numeric(0),
      stringsAsFactors = FALSE
    )
    out[[id_col]] <- character(0)
    out <- out[, c("variable_id", id_col, "count"), drop = FALSE]
    return(out)
  }
  
  cols <- intersect(c("variable_id", id_col, "count"), colnames(df))
  out <- df[, cols, drop = FALSE]
  if (!"variable_id" %in% colnames(out)) out$variable_id <- NA_character_
  if (!id_col %in% colnames(out)) out[[id_col]] <- NA_character_
  if (!"count"       %in% colnames(out)) out$count       <- NA_real_
  
  out$variable_id <- as.character(out$variable_id)
  out[[id_col]]   <- as.character(out[[id_col]])
  out$count       <- as.numeric(out$count)
  
  out <- out[, c("variable_id", id_col, "count"), drop = FALSE]
  out <- out[order(out$variable_id, out[[id_col]], out$count), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Check if counts are equal
#' @noRd
counts_equal <- function(a, b, id_col = "KEGG_ID") {
  identical(canon_counts(a, id_col = id_col), canon_counts(b, id_col = id_col))
}
