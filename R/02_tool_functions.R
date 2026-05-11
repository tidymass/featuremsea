#' @title Internal Helper Functions
#' @description
#' These functions are internal utilities for processing MFM results,
#' handling annotation tables, and calculating scores.
#' They are not intended to be called directly by the user.
#'
#' @name internal_utils
#' @keywords internal
#'
#' @importFrom dplyr mutate filter pull count group_by summarise n select distinct rename arrange
#' @importFrom purrr imap_dfr map_dfr
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom rlang %||% sym
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
#' Optimized version using matrix operations. Applies bonus weights based on
#' feature-metabolite counts, with upper bound constraint to prevent scores
#' from exceeding the row-wise maximum original score.
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

  # --- Step 1: Global max score (streaming over columns, no dense copy) ---
  max_score_global <- -Inf
  for (j in seq_len(n_col)) {
    col_j <- weighted_df[[j]]
    if (!is.numeric(col_j)) {
      col_j <- suppressWarnings(as.numeric(col_j))
    }
    col_for_max <- col_j
    col_for_max[is.na(col_for_max)] <- -Inf
    max_score_global <- max(max_score_global, col_for_max, na.rm = TRUE)
    weighted_df[[j]] <- col_j
  }
  max_score_global <- if (!is.finite(max_score_global)) 0 else max_score_global
  
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
  bonus_vec <- counts_vec / denom[row_indices] * max_score_global

  # --- Step 5: Apply sparse updates grouped by column with upper bound constraint ---
  split_idx <- split(seq_along(col_indices), col_indices)
  for (col_key in names(split_idx)) {
    idx <- split_idx[[col_key]]
    cj <- as.integer(col_key)
    rows <- row_indices[idx]
    addv <- bonus_vec[idx]
    col_j <- weighted_df[[cj]]

    # Apply bonus but constrain by global maximum score
    new_values <- col_j[rows] + addv
    col_j[rows] <- pmin(new_values, max_score_global)

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

#' Extract leading edge feature-metabolite pairs with weighted scores
#'
#' This function extracts all feature-metabolite pairs from the leading edge of
#' significant pathways and retrieves their corresponding weighted scores from
#' the annotation table.
#'
#' @param result_object A featuremsea_object containing analysis results.
#' @param id_col Character. Column name for metabolite ID (default: "KEGG_ID").
#' @param fdr_threshold Numeric. FDR threshold for significance (default: 0.05).
#'
#' @return A data frame with three columns:
#'   \item{feature}{Character. Feature IDs (variable_id)}
#'   \item{metabolite}{Character. Metabolite IDs}
#'   \item{score}{Numeric. Weighted scores from annotation_table_weighting}
#'
#' @examples
#' \dontrun{
#' # Extract leading edge pairs with scores
#' leading_edge_scores <- get_leading_edge_scores(result_object)
#'
#' # With custom parameters
#' leading_edge_scores <- get_leading_edge_scores(
#'   result_object,
#'   id_col = "CompoundID",
#'   fdr_threshold = 0.01
#' )
#' }
#' @export
get_leading_edge_scores <- function(result_object, id_col = "KEGG_ID", fdr_threshold = 0.05) {

  # 检查输入
  if (!inherits(result_object, "featuremsea_object")) {
    stop("result_object must be a featuremsea_object")
  }

  # 获取significant modules
  significant_modules <- result_object@significant_modules
  if (nrow(significant_modules) == 0) {
    warning("No significant modules found")
    return(data.frame(
      feature = character(0),
      metabolite = character(0),
      score = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  # 获取significant pathway IDs (根据FDR阈值进一步筛选)
  sig_pathways <- significant_modules %>%
    filter(!is.na(FDR), FDR < fdr_threshold)

  if (nrow(sig_pathways) == 0) {
    warning("No pathways pass the FDR threshold")
    return(data.frame(
      feature = character(0),
      metabolite = character(0),
      score = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  # 获取res_list和annotation_table_weighting
  res_list <- result_object@res_list
  annotation_table_weighting <- result_object@annotation_table_weighting

  # 提取所有significant pathway的leading edge
  all_leading_edge <- map_dfr(
    sig_pathways$pathway_id,  # 使用pathway_id作为key
    function(pathway_id) {
      le <- res_list[[pathway_id]]$leading_edge
      if (!is.null(le) && nrow(le) > 0) {
        # 选择需要的列：variable_id和代谢物ID列
        le_subset <- le[, c("variable_id", id_col), drop = FALSE]
        le_subset$pathway_id <- pathway_id  # 保留pathway信息便于调试
        return(le_subset)
      } else {
        return(NULL)
      }
    }
  )

  if (nrow(all_leading_edge) == 0) {
    warning("No leading edge data found in significant pathways")
    return(data.frame(
      feature = character(0),
      metabolite = character(0),
      score = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  # 去重（同一个feature-metabolite对可能在多个pathway的leading edge中出现）
  unique_pairs <- all_leading_edge %>%
    select(variable_id, !!sym(id_col)) %>%
    distinct()

  # 从annotation_table_weighting中提取对应的score
  # annotation_table_weighting的行名是variable_id，列名是metabolite_id

  # 初始化结果dataframe
  result_df <- unique_pairs
  result_df$score <- NA_real_

  # 获取有效的行名和列名
  valid_features <- rownames(annotation_table_weighting)
  valid_metabolites <- colnames(annotation_table_weighting)

  # 遍历每一行，提取对应的score
  for (i in seq_len(nrow(result_df))) {
    feat <- result_df$variable_id[i]
    metab <- result_df[[id_col]][i]

    if (!is.na(feat) && !is.na(metab) &&
        feat %in% valid_features && metab %in% valid_metabolites) {
      result_df$score[i] <- annotation_table_weighting[feat, metab]
    }
  }

  # 重命名列并清理数据
  result_df <- result_df %>%
    rename(
      feature = variable_id,
      metabolite = !!sym(id_col)
    ) %>%
    # 移除缺失值
    filter(!is.na(score)) %>%
    # score取整（四舍五入）
    mutate(score = round(score)) %>%
    # 按score降序排列
    arrange(desc(score))

  return(result_df)
}
