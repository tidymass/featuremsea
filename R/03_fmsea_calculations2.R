# ============================================================
# Core calculation functions for featureMSEA / mMSEA
# Modified Strategy:
# 1. Feature Level: Max Score (a_ij) per Feature.
# 2. Metabolite Level: Top Rank (First occurrence) per Metabolite.
# 3. Information: Retain all features mapped to a metabolite in 'all_mapped_features'.
# ============================================================

# ------------------------------------------------------------
# 1. Normalize the annotation matrix
# ------------------------------------------------------------
normalize_score_annotation_global <- function(score_annotation) {
  # 1. Convert input to a numeric matrix
  # Ensure storage mode is double to prevent integer division issues
  num_mat <- as.matrix(score_annotation)
  storage.mode(num_mat) <- "double"
  
  # 2. Calculate the global maximum value across the entire matrix
  max_val <- max(num_mat, na.rm = TRUE)
  
  # 3. Safety check: prevent division by zero or infinity
  if (!is.finite(max_val) || max_val == 0) {
    max_val <- 1
  }
  
  # 4. Perform global normalization
  num_mat / max_val
}

# ------------------------------------------------------------
# 2. Convert annotation matrix to long format
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
# 3. Split pathway IDs separated by "{}"
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
# 4. Build ranking index from ranking table
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
    V              = V,
    var_ids        = var_ids,
    weights_by_pos = weights,
    var2pos        = var2pos
  )
}

# ------------------------------------------------------------
# 5. Precompute static information for mMSEA
# [MODIFIED] Strategy:
# 1. Feature Deduplication: Keep max a_ij per feature.
# 2. Metabolite Deduplication: Keep best ranked feature per metabolite.
# 3. Aggregation: Store all mapped features in 'all_mapped_features'.
# ------------------------------------------------------------
precompute_mMSEA_static <- function(pathway_ids_vec,
                                    annotation_long,
                                    rk_idx,
                                    id_col = "KEGG_ID") {
  id_vec  <- annotation_long[[id_col]]
  vid_all <- annotation_long[["variable_id"]]
  a_all   <- annotation_long[["a_ij"]]
  
  # --- 1. Basic Filtering ---
  hits_ids <- intersect(pathway_ids_vec, unique(id_vec))
  if (length(hits_ids) == 0L) return(NULL)
  
  in_hits <- match(id_vec, hits_ids, nomatch = 0L) > 0L
  in_rank <- !is.na(match(vid_all, rk_idx$var_ids))
  keep    <- in_hits & in_rank & !is.na(a_all) & (a_all > 0)
  
  if (!any(keep)) return(NULL)
  
  vid_pairs <- vid_all[keep]
  id_pairs  <- id_vec[keep]
  a_pairs   <- a_all[keep]
  
  # --- 2. Deduplication Strategy ---
  
  # Create temporary DF for processing
  temp_df <- data.frame(
    vid = vid_pairs,
    id  = id_pairs,
    a   = a_pairs,
    stringsAsFactors = FALSE
  )
  
  # A. Feature Level Strategy: Keep Highest Score
  # Sort by score descending, then remove duplicate vids
  temp_df <- temp_df[order(temp_df$a, decreasing = TRUE), ]
  temp_df <- temp_df[!duplicated(temp_df$vid), ]
  
  # B. Preparation for Metabolite Strategy
  # Get Rank Position for each surviving feature (Smaller index = Better rank)
  temp_df$rank_pos <- rk_idx$var2pos[temp_df$vid]
  
  # Sort by Rank Position Ascending (Best rank first)
  temp_df <- temp_df[order(temp_df$rank_pos, decreasing = FALSE), ]
  
  # C. Aggregate Info (Mapping all valid features to the metabolite)
  # Aggregating based on 'id'. Since temp_df is filtered, this lists all features
  # that "won" the feature-level step for this metabolite.
  agg_res <- aggregate(vid ~ id, data = temp_df, FUN = function(x) {
    paste(unique(x), collapse = "; ")
  })
  colnames(agg_res)[2] <- "all_mapped_features"
  
  # D. Metabolite Level Strategy: Keep First Occurrence (Best Rank)
  # Since we sorted by rank_pos, duplicated() keeps the first one.
  temp_df_clean <- temp_df[!duplicated(temp_df$id), ]
  
  # Merge aggregated info back
  temp_df_clean <- merge(temp_df_clean, agg_res, by = "id", sort = FALSE)
  
  # Re-sort to ensure correct GSEA order (Merge might disturb order)
  temp_df_clean <- temp_df_clean[order(temp_df_clean$rank_pos, decreasing = FALSE), ]
  
  # Update main vectors
  vid_pairs           <- temp_df_clean$vid
  id_pairs            <- temp_df_clean$id
  a_pairs             <- temp_df_clean$a
  all_mapped_features <- temp_df_clean$all_mapped_features
  
  # --- 3. Build Structural Vectors for Permutation/ES ---
  
  # Count pairs (After deduplication, most n_pairs will be 1, rest 0)
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
  sorted_all_features <- rep(NA_character_, sum(n_pairs)) # New vector for storage
  
  # Map lookup for fast filling
  map_lookup <- setNames(all_mapped_features, vid_pairs)
  
  pp <- 1L
  mp <- 1L
  
  for (i in seq_along(rk_idx$var_ids)) {
    s <- starts[i]
    curr_vid <- rk_idx$var_ids[i]
    
    if (n_pairs[i] > 0L) {
      len <- n_pairs[i]
      # Fill hit positions
      pair_positions[pp:(pp + len - 1L)] <- s:(s + len - 1L)
      # Fill aggregated feature info
      sorted_all_features[pp:(pp + len - 1L)] <- map_lookup[curr_vid]
      
      pp <- pp + len
    } else {
      # Fill miss positions
      miss_positions[mp]   <- s
      miss_variable_id[s]  <- curr_vid
      mp <- mp + 1L
    }
  }
  
  list(
    vid_pairs           = vid_pairs,
    id_pairs            = id_pairs,
    a_pairs             = a_pairs,
    sorted_all_features = sorted_all_features, # Passed to compute_ES
    total_slots         = total_slots,
    pair_positions      = pair_positions,
    miss_positions      = miss_positions,
    miss_variable_id    = miss_variable_id,
    V                   = rk_idx$V,
    var_ids             = rk_idx$var_ids,
    n_pairs             = n_pairs,
    slot_count          = slot_count
  )
}


# ------------------------------------------------------------
# 6. Compute ES given a ranking mapping
# [MODIFIED] Adds 'all_mapped_features' to leading_edge output
# ------------------------------------------------------------
compute_ES_from_mapping <- function(static,
                                    weights_by_pos,
                                    var2pos,
                                    id_col = "KEGG_ID",
                                    return_details = FALSE,
                                    pair_positions_override = NULL,
                                    miss_positions_override = NULL,
                                    miss_variable_id_override = NULL) {
  
  # Extract Rank Weights
  pos <- unname(var2pos[static$vid_pairs])
  wgt <- weights_by_pos[pos]
  
  # Calculate Pair Weight (a_ij * rank_weight)
  # Note: a_pairs no longer has redundancy penalty, it's just raw score
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
  
  # Filter Data
  pw          <- pw[keep]
  vhk         <- static$vid_pairs[keep]
  idk         <- static$id_pairs[keep]
  all_feats_k <- static$sorted_all_features[keep] # Get stored feature strings
  
  # Sort by Contribution (Standard GSEA logic)
  ord         <- order(pw, decreasing = TRUE)
  top_idx_all <- ord
  
  # Determine Positions
  pair_pos <- if (is.null(pair_positions_override)) static$pair_positions else pair_positions_override
  miss_pos <- if (is.null(miss_positions_override)) static$miss_positions else miss_positions_override
  
  total_slots <- length(pair_pos) + length(miss_pos)
  n_fill      <- min(length(pair_pos), length(ord))
  step        <- numeric(total_slots)
  
  # Edge Case: No hits filled
  if (n_fill == 0L) {
    if (length(miss_pos) > 0L) {
      step[miss_pos] <- -1 / length(miss_pos)
    }
    ES_cum <- cumsum(step)
    return(list(ES = max(ES_cum)))
  }
  
  top_idx <- top_idx_all[seq_len(n_fill)]
  top_pw  <- pw[top_idx]
  
  # Calculate Hit Step
  hit_total <- sum(top_pw)
  if (!is.finite(hit_total) || hit_total <= 0) {
    return(list(ES = NA_real_))
  }
  step[pair_pos[seq_len(n_fill)]] <- top_pw / hit_total
  
  # Calculate Miss Step
  N_miss <- length(miss_pos)
  if (N_miss > 0L) {
    step[miss_pos] <- -1 / N_miss
  }
  
  ES_cum <- cumsum(step)
  ES_val <- max(ES_cum)
  
  if (!return_details) {
    return(list(ES = ES_val))
  }
  
  # --- Details: Leading Edge ---
  max_idx <- which.max(ES_cum)
  positions_assigned <- pair_pos[seq_len(n_fill)]
  keep_le <- positions_assigned <= max_idx
  
  leading_edge <- data.frame(
    variable_id         = vhk[top_idx][keep_le],        # The specific feature that contributed
    tmp_id              = idk[top_idx][keep_le],        # The metabolite ID
    all_mapped_features = all_feats_k[top_idx][keep_le],# [NEW] All features mapping to this met
    pair_weight         = top_pw[keep_le],
    a_ij                = static$a_pairs[keep][top_idx][keep_le],
    ranking_weight      = wgt[keep][top_idx][keep_le],
    stringsAsFactors    = FALSE
  )
  names(leading_edge)[names(leading_edge) == "tmp_id"] <- id_col
  
  # --- Details: Steps Table ---
  type <- rep("miss", total_slots)
  type[pair_pos] <- "pair_slot"
  
  variable_id_miss <- if (is.null(miss_variable_id_override)) static$miss_variable_id else miss_variable_id_override
  variable_id_hit  <- rep(NA_character_, total_slots)
  id_hit           <- rep(NA_character_, total_slots)
  
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
# 7. Main function: mMSEA with C++ permutation back-end
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
  # NOTE: Expects 'compute_permutations_cpp' to be available in environment/package
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