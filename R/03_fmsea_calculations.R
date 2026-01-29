# ============================================================
# Core calculation functions for featureMSEA
# ============================================================

# ------------------------------------------------------------
# Normalize the annotation matrix
# ------------------------------------------------------------
# normalize_score_annotation_new <- function(score_annotation) {
#   num_mat <- as.matrix(score_annotation)
#   storage.mode(num_mat) <- "double"
#   
#   Sj <- colSums(num_mat, na.rm = TRUE)
#   Sj[!is.finite(Sj) | Sj == 0] <- 1
#   
#   sweep(num_mat, 2, Sj, "/")
# }

normalize_score_annotation_global <- function(score_annotation) {
  # 1. Convert input to a numeric matrix
  # Ensure storage mode is double to prevent integer division issues
  num_mat <- as.matrix(score_annotation)
  storage.mode(num_mat) <- "double"
  
  # 2. Calculate the global maximum value across the entire matrix
  max_val <- max(num_mat, na.rm = TRUE)
  
  # 3. Safety check: prevent division by zero or infinity
  # If max_val is 0 (e.g., matrix of all zeros) or non-finite, set denominator to 1
  # This ensures the original values are returned without generating NaNs or Infs
  if (!is.finite(max_val) || max_val == 0) {
    max_val <- 1
  }
  
  # 4. Perform global normalization
  # Divide every element in the matrix by the global maximum value
  num_mat / max_val
}

# ------------------------------------------------------------
# Convert annotation matrix to long format
# ------------------------------------------------------------
annotation_long_fast_base <- function(annotation_table, id_col = "KEGG_ID") {
  rn <- rownames(annotation_table)
  cn <- colnames(annotation_table)
  
  if (is.null(rn) || is.null(cn)) {
    stop("annotation_table must have both row and column names")
  }
  
  m <- as.matrix(annotation_table)
  storage.mode(m) <- "double"
  
  rs_keep <- rowSums(m > 0, na.rm = TRUE) > 0
  cs_keep <- colSums(m > 0, na.rm = TRUE) > 0
  
  m  <- m[rs_keep, cs_keep, drop = FALSE]
  rn <- rn[rs_keep]
  cn <- cn[cs_keep]
  
  idx <- which(m > 0, arr.ind = TRUE)
  if (nrow(idx) == 0L) {
    out <- data.frame(
      variable_id = character(),
      tmp         = character(),
      a_ij        = numeric(),
      stringsAsFactors = FALSE
    )
    names(out)[2] <- id_col
    return(out)
  }
  
  out <- data.frame(
    variable_id = rn[idx[, 1]],
    tmp         = cn[idx[, 2]],
    a_ij        = m[idx],
    stringsAsFactors = FALSE
  )
  names(out)[2] <- id_col
  out
}

# ------------------------------------------------------------
# Split pathway IDs separated by "{}"
# ------------------------------------------------------------
split_ids <- function(x) {
  if (is.na(x) || x == "") {
    return(character())
  }
  res <- stringr::str_split(x, "\\{\\}")[[1]]
  res <- unique(res)
  res[!is.na(res) & res != ""]
}

# ------------------------------------------------------------
# Build ranking index from ranking table
# ranking_table must contain columns: variable_id, ranking_weight
# ------------------------------------------------------------
build_ranking_index <- function(ranking_table) {
  # Keep only needed columns
  if (!all(c("variable_id", "ranking_weight") %in% colnames(ranking_table))) {
    stop("ranking_table must contain 'variable_id' and 'ranking_weight' columns.")
  }
  
  rt <- ranking_table[, c("variable_id", "ranking_weight")]
  
  # Remove duplicate rows (distinct variable_id, ranking_weight)
  base_rk <- rt[!duplicated(rt), , drop = FALSE]
  
  # Order by ranking_weight descending
  ord <- order(base_rk[["ranking_weight"]], decreasing = TRUE)
  base_rk <- base_rk[ord, , drop = FALSE]
  
  var_ids <- base_rk[["variable_id"]]
  weights <- base_rk[["ranking_weight"]]
  V       <- length(var_ids)
  
  var2pos <- stats::setNames(seq_len(V), var_ids)
  
  list(
    V             = V,
    var_ids       = var_ids,
    weights_by_pos = weights,
    var2pos       = var2pos
  )
}

# # ------------------------------------------------------------
# # Precompute static information for mMSEA
# # ------------------------------------------------------------
# precompute_mMSEA_static <- function(pathway_ids_vec,
#                                     annotation_long,
#                                     rk_idx,
#                                     id_col = "KEGG_ID") {
#   id_vec  <- annotation_long[[id_col]]
#   vid_all <- annotation_long[["variable_id"]]
#   a_all   <- annotation_long[["a_ij"]]
#   
#   hits_ids <- intersect(pathway_ids_vec, unique(id_vec))
#   if (length(hits_ids) == 0L) {
#     return(NULL)
#   }
#   
#   in_hits <- match(id_vec, hits_ids, nomatch = 0L) > 0L
#   in_rank <- !is.na(match(vid_all, rk_idx$var_ids))
#   keep    <- in_hits & in_rank & !is.na(a_all) & (a_all > 0)
#   
#   if (!any(keep)) {
#     return(NULL)
#   }
#   
#   vid_pairs <- vid_all[keep]
#   id_pairs  <- id_vec[keep]
#   a_pairs   <- a_all[keep]
#   
#   n_pairs_tab <- table(vid_pairs)
#   n_pairs <- integer(length(rk_idx$var_ids))
#   names(n_pairs) <- rk_idx$var_ids
#   n_pairs[names(n_pairs_tab)] <- as.integer(n_pairs_tab)
#   
#   slot_count <- ifelse(n_pairs > 0L, n_pairs, 1L)
#   total_slots <- sum(slot_count)
#   
#   ends   <- cumsum(slot_count)
#   starts <- ends - slot_count + 1L
#   
#   pair_positions     <- integer(sum(n_pairs))
#   miss_positions     <- integer(sum(slot_count == 1L & n_pairs == 0L))
#   miss_variable_id   <- rep(NA_character_, total_slots)
#   
#   pp <- 1L
#   mp <- 1L
#   
#   for (i in seq_along(rk_idx$var_ids)) {
#     s <- starts[i]
#     if (n_pairs[i] > 0L) {
#       len <- n_pairs[i]
#       idx <- s:(s + len - 1L)
#       pair_positions[pp:(pp + len - 1L)] <- idx
#       pp <- pp + len
#     } else {
#       miss_positions[mp]   <- s
#       miss_variable_id[s]  <- rk_idx$var_ids[i]
#       mp <- mp + 1L
#     }
#   }
#   
#   list(
#     vid_pairs       = vid_pairs,
#     id_pairs        = id_pairs,
#     a_pairs         = a_pairs,
#     total_slots     = total_slots,
#     pair_positions  = pair_positions,
#     miss_positions  = miss_positions,
#     miss_variable_id = miss_variable_id,
#     V               = rk_idx$V,
#     var_ids         = rk_idx$var_ids,
#     n_pairs         = n_pairs,
#     slot_count      = slot_count
#   )
# }

# ------------------------------------------------------------
# Precompute static information for mMSEA
# Double Weighting Strategy is PERMANENTLY ENABLED
# ------------------------------------------------------------
precompute_mMSEA_static <- function(pathway_ids_vec,
                                    annotation_long,
                                    rk_idx,
                                    id_col = "KEGG_ID") { # [移除] double_weighting 参数
  id_vec  <- annotation_long[[id_col]]
  vid_all <- annotation_long[["variable_id"]]
  a_all   <- annotation_long[["a_ij"]]
  
  hits_ids <- intersect(pathway_ids_vec, unique(id_vec))
  if (length(hits_ids) == 0L) {
    return(NULL)
  }
  
  in_hits <- match(id_vec, hits_ids, nomatch = 0L) > 0L
  in_rank <- !is.na(match(vid_all, rk_idx$var_ids))
  keep    <- in_hits & in_rank & !is.na(a_all) & (a_all > 0)
  
  if (!any(keep)) {
    return(NULL)
  }
  
  vid_pairs <- vid_all[keep]
  id_pairs  <- id_vec[keep]
  a_pairs   <- a_all[keep]
  
  # ============================================================
  # Double Weighting Strategy (Always Executed)
  # ============================================================
  # 1. 计算全局代谢物冗余度 (即每个代谢物对应多少个 Feature)
  # 使用整个 annotation 表的 id_vec 进行统计，确保反映真实的"一对多"情况
  met_counts_all <- table(id_vec)
  
  # 2. 获取当前筛选出的 pairs 对应的代谢物计数
  # 使用 character indexing 快速匹配
  counts_for_pairs <- as.numeric(met_counts_all[id_pairs])
  
  # 3. 计算惩罚权重 w_j = 1 / sqrt(N_j)
  # 既惩罚了冗余，又保留了 Feature 内部的相对打分差异
  w_j <- 1 / sqrt(counts_for_pairs)
  
  # 4. 更新注释分数
  # a_pairs 在这里被直接修改，后续所有计算(ES/Permutation)均基于此加权值
  a_pairs <- a_pairs * w_j
  # ============================================================
  
  n_pairs_tab <- table(vid_pairs)
  n_pairs <- integer(length(rk_idx$var_ids))
  names(n_pairs) <- rk_idx$var_ids
  n_pairs[names(n_pairs_tab)] <- as.integer(n_pairs_tab)
  
  slot_count <- ifelse(n_pairs > 0L, n_pairs, 1L)
  total_slots <- sum(slot_count)
  
  ends   <- cumsum(slot_count)
  starts <- ends - slot_count + 1L
  
  pair_positions      <- integer(sum(n_pairs))
  miss_positions      <- integer(sum(slot_count == 1L & n_pairs == 0L))
  miss_variable_id    <- rep(NA_character_, total_slots)
  
  pp <- 1L
  mp <- 1L
  
  for (i in seq_along(rk_idx$var_ids)) {
    s <- starts[i]
    if (n_pairs[i] > 0L) {
      len <- n_pairs[i]
      idx <- s:(s + len - 1L)
      pair_positions[pp:(pp + len - 1L)] <- idx
      pp <- pp + len
    } else {
      miss_positions[mp]   <- s
      miss_variable_id[s]  <- rk_idx$var_ids[i]
      mp <- mp + 1L
    }
  }
  
  list(
    vid_pairs       = vid_pairs,
    id_pairs        = id_pairs,
    a_pairs         = a_pairs, # 始终包含二次加权逻辑
    total_slots     = total_slots,
    pair_positions  = pair_positions,
    miss_positions  = miss_positions,
    miss_variable_id = miss_variable_id,
    V               = rk_idx$V,
    var_ids         = rk_idx$var_ids,
    n_pairs         = n_pairs,
    slot_count      = slot_count
  )
}


# # ------------------------------------------------------------
# # Compute ES given a ranking mapping
# # (with leading-edge truncation at ES peak)
# # ------------------------------------------------------------
# compute_ES_from_mapping <- function(static,
#                                     weights_by_pos,
#                                     var2pos,
#                                     id_col = "KEGG_ID",
#                                     return_details = FALSE,
#                                     pair_positions_override = NULL,
#                                     miss_positions_override = NULL,
#                                     miss_variable_id_override = NULL) {
#   pos <- unname(var2pos[static$vid_pairs])
#   wgt <- weights_by_pos[pos]
#   
#   pw   <- static$a_pairs * wgt
#   keep <- !is.na(pw) & (pw > 0)
#   
#   if (!any(keep)) {
#     out <- list(ES = NA_real_)
#     if (return_details) {
#       out[["leading_edge"]] <- NULL
#       out[["steps"]]        <- NULL
#     }
#     return(out)
#   }
#   
#   pw  <- pw[keep]
#   vhk <- static$vid_pairs[keep]
#   idk <- static$id_pairs[keep]
#   
#   ord         <- order(pw, decreasing = TRUE)
#   top_idx_all <- ord
#   
#   pair_pos <- if (is.null(pair_positions_override)) {
#     static$pair_positions
#   } else {
#     pair_positions_override
#   }
#   
#   miss_pos <- if (is.null(miss_positions_override)) {
#     static$miss_positions
#   } else {
#     miss_positions_override
#   }
#   
#   total_slots <- length(pair_pos) + length(miss_pos)
#   
#   n_fill <- min(length(pair_pos), length(ord))
#   step  <- numeric(total_slots)
#   
#   if (n_fill == 0L) {
#     if (length(miss_pos) > 0L) {
#       step[miss_pos] <- -1 / length(miss_pos)
#     }
#     ES_cum <- cumsum(step)
#     return(list(ES = max(ES_cum)))
#   }
#   
#   top_idx <- top_idx_all[seq_len(n_fill)]
#   top_pw  <- pw[top_idx]
#   
#   hit_total <- sum(top_pw)
#   if (!is.finite(hit_total) || hit_total <= 0) {
#     return(list(ES = NA_real_))
#   }
#   
#   step[pair_pos[seq_len(n_fill)]] <- top_pw / hit_total
#   
#   N_miss <- length(miss_pos)
#   if (N_miss > 0L) {
#     step[miss_pos] <- -1 / N_miss
#   }
#   
#   ES_cum <- cumsum(step)
#   ES_val <- max(ES_cum)
#   
#   if (!return_details) {
#     return(list(ES = ES_val))
#   }
#   
#   # Leading edge
#   max_idx <- which.max(ES_cum)
#   positions_assigned <- pair_pos[seq_len(n_fill)]
#   keep_le <- positions_assigned <= max_idx
#   
#   leading_edge <- data.frame(
#     variable_id    = vhk[top_idx][keep_le],
#     tmp_id         = idk[top_idx][keep_le],
#     pair_weight    = top_pw[keep_le],
#     a_ij           = static$a_pairs[keep][top_idx][keep_le],
#     ranking_weight = wgt[keep][top_idx][keep_le],
#     stringsAsFactors = FALSE
#   )
#   names(leading_edge)[names(leading_edge) == "tmp_id"] <- id_col
#   
#   type <- rep("miss", total_slots)
#   type[pair_pos] <- "pair_slot"
#   
#   variable_id_miss <- if (is.null(miss_variable_id_override)) {
#     static$miss_variable_id
#   } else {
#     miss_variable_id_override
#   }
#   
#   variable_id_hit <- rep(NA_character_, total_slots)
#   id_hit          <- rep(NA_character_, total_slots)
#   
#   variable_id_hit[pair_pos[seq_len(n_fill)]] <- vhk[top_idx]
#   id_hit[pair_pos[seq_len(n_fill)]]          <- idk[top_idx]
#   
#   steps <- data.frame(
#     order_pos        = seq_len(total_slots),
#     type             = type,
#     variable_id_miss = variable_id_miss,
#     variable_id_hit  = variable_id_hit,
#     tmp_id           = id_hit,
#     step             = step,
#     ES_cum           = ES_cum,
#     stringsAsFactors = FALSE
#   )
#   names(steps)[names(steps) == "tmp_id"] <- id_col
#   
#   list(
#     ES           = ES_val,
#     leading_edge = leading_edge,
#     steps        = steps
#   )
# }

# ------------------------------------------------------------
# Compute ES given a ranking mapping
# (Modified: Feature only contributes score on first occurrence)
# ------------------------------------------------------------
compute_ES_from_mapping <- function(static,
                                    weights_by_pos,
                                    var2pos,
                                    id_col = "KEGG_ID",
                                    return_details = FALSE,
                                    pair_positions_override = NULL,
                                    miss_positions_override = NULL,
                                    miss_variable_id_override = NULL) {
  pos <- unname(var2pos[static$vid_pairs])
  wgt <- weights_by_pos[pos]
  
  # 原始 Pair Weight
  pw   <- static$a_pairs * wgt
  keep <- !is.na(pw) & (pw > 0)
  
  if (!any(keep)) {
    out <- list(ES = NA_real_)
    if (return_details) {
      out[["leading_edge"]] <- NULL
      out[["steps"]]        <- NULL
    }
    return(out)
  }
  
  pw  <- pw[keep]
  vhk <- static$vid_pairs[keep]
  idk <- static$id_pairs[keep]
  
  # 1. 按照权重从大到小排序
  ord         <- order(pw, decreasing = TRUE)
  top_idx_all <- ord
  
  pair_pos <- if (is.null(pair_positions_override)) {
    static$pair_positions
  } else {
    pair_positions_override
  }
  
  miss_pos <- if (is.null(miss_positions_override)) {
    static$miss_positions
  } else {
    miss_positions_override
  }
  
  total_slots <- length(pair_pos) + length(miss_pos)
  n_fill <- min(length(pair_pos), length(ord))
  step   <- numeric(total_slots)
  
  if (n_fill == 0L) {
    if (length(miss_pos) > 0L) {
      step[miss_pos] <- -1 / length(miss_pos)
    }
    ES_cum <- cumsum(step)
    return(list(ES = max(ES_cum)))
  }
  
  # --- [新增/修改逻辑开始] ---
  
  # 获取排序后的 Feature ID
  top_idx     <- top_idx_all[seq_len(n_fill)]
  sorted_vids <- vhk[top_idx]
  
  # 识别重复的 Feature (sorted_vids 中第二次及以后出现的标记为 TRUE)
  is_dup <- duplicated(sorted_vids)
  
  # 获取原始排序后的权重
  top_pw_raw <- pw[top_idx]
  
  # 创建"有效权重" (effective weights)
  # 如果是重复出现的 Feature，将其有效权重设为 0
  top_pw_effective <- top_pw_raw
  top_pw_effective[is_dup] <- 0
  
  # 计算总 Hit 分数时，只使用有效权重（即去重后的总分）
  # 这样分母更合理，避免被重复项稀释
  hit_total <- sum(top_pw_effective)
  
  if (!is.finite(hit_total) || hit_total <= 0) {
    return(list(ES = NA_real_))
  }
  
  # 计算步长：使用有效权重计算
  # 重复的 Feature 这里分子是 0，步长也是 0，曲线会走平路
  step[pair_pos[seq_len(n_fill)]] <- top_pw_effective / hit_total
  
  # --- [新增/修改逻辑结束] ---
  
  N_miss <- length(miss_pos)
  if (N_miss > 0L) {
    step[miss_pos] <- -1 / N_miss
  }
  
  ES_cum <- cumsum(step)
  ES_val <- max(ES_cum)
  
  if (!return_details) {
    return(list(ES = ES_val))
  }
  
  # Leading edge
  max_idx <- which.max(ES_cum)
  positions_assigned <- pair_pos[seq_len(n_fill)]
  keep_le <- positions_assigned <= max_idx
  
  leading_edge <- data.frame(
    variable_id    = vhk[top_idx][keep_le],
    tmp_id         = idk[top_idx][keep_le],
    
    # 这里输出 effective weight，这样用户能直观看到重复项的贡献是 0
    # 如果想保留原始计算权重供参考，可以增加一列 original_weight
    pair_weight    = top_pw_effective[keep_le], 
    original_weight = top_pw_raw[keep_le],     # 新增：原始权重供参考
    
    a_ij           = static$a_pairs[keep][top_idx][keep_le],
    ranking_weight = wgt[keep][top_idx][keep_le],
    is_duplicate   = is_dup[keep_le],          # 新增：标记是否因为重复被归零
    stringsAsFactors = FALSE
  )
  names(leading_edge)[names(leading_edge) == "tmp_id"] <- id_col
  
  type <- rep("miss", total_slots)
  type[pair_pos] <- "pair_slot"
  
  variable_id_miss <- if (is.null(miss_variable_id_override)) {
    static$miss_variable_id
  } else {
    miss_variable_id_override
  }
  
  variable_id_hit <- rep(NA_character_, total_slots)
  id_hit          <- rep(NA_character_, total_slots)
  
  variable_id_hit[pair_pos[seq_len(n_fill)]] <- vhk[top_idx]
  id_hit[pair_pos[seq_len(n_fill)]]          <- idk[top_idx]
  
  steps <- data.frame(
    order_pos        = seq_len(total_slots),
    type             = type,
    variable_id_miss = variable_id_miss,
    variable_id_hit  = variable_id_hit,
    tmp_id           = id_hit,
    step             = step,
    ES_cum           = ES_cum,
    stringsAsFactors = FALSE
  )
  names(steps)[names(steps) == "tmp_id"] <- id_col
  
  list(
    ES           = ES_val,
    leading_edge = leading_edge,
    steps        = steps
  )
}



# ------------------------------------------------------------
# Main function: mMSEA with C++ permutation back-end
# ------------------------------------------------------------
get_mMSEA_results_indexed_fast <- function(
    pathway_ids_vec,
    annotation_long,
    rk_idx,
    id_col    = "KEGG_ID",
    n_perm    = 1000L,
    seed      = 123L,
    return_perm = FALSE
) {
  set.seed(seed)
  
  static <- precompute_mMSEA_static(
    pathway_ids_vec = pathway_ids_vec,
    annotation_long = annotation_long,
    rk_idx          = rk_idx,
    id_col          = id_col
  )
  
  if (is.null(static)) {
    out <- list(
      ES          = NA_real_,
      leading_edge = NULL,
      p_value     = NA_real_,
      NES         = NA_real_,
      steps       = NULL
    )
    if (return_perm) {
      out[["perm_ES"]] <- numeric(0L)
    }
    return(out)
  }
  
  # Observed ES with fixed ranking
  obs <- compute_ES_from_mapping(
    static          = static,
    weights_by_pos  = rk_idx$weights_by_pos,
    var2pos         = rk_idx$var2pos,
    id_col          = id_col,
    return_details  = TRUE
  )
  
  ES_obs <- obs$ES
  if (is.na(ES_obs)) {
    out <- list(
      ES           = NA_real_,
      leading_edge = NULL,
      p_value      = NA_real_,
      NES          = NA_real_,
      steps        = NULL
    )
    if (return_perm) {
      out[["perm_ES"]] <- numeric(0L)
    }
    return(out)
  }
  
  # Use C++ to compute permutation ES values
  vid_pairs_idx <- match(static$vid_pairs, rk_idx$var_ids) - 1L
  
  perm_ES <- compute_permutations_cpp(
    vid_pairs_idx = as.integer(vid_pairs_idx),
    a_pairs       = static$a_pairs,
    n_pairs       = as.integer(static$n_pairs),
    slot_count    = as.integer(static$slot_count),
    base_weights  = rk_idx$weights_by_pos,
    V             = as.integer(static$V),
    n_perm        = as.integer(n_perm),
    seed          = as.integer(seed)
  )
  
  valid <- !is.na(perm_ES)
  p_val <- (sum(perm_ES[valid] >= ES_obs) + 1) / (sum(valid) + 1)
  NES   <- ES_obs / mean(abs(perm_ES[valid]), na.rm = TRUE)
  
  res <- list(
    ES           = ES_obs,
    leading_edge = obs$leading_edge,
    p_value      = p_val,
    NES          = NES,
    steps        = obs$steps
  )
  
  if (return_perm) {
    res[["perm_ES"]] <- perm_ES
  }
  
  res
}