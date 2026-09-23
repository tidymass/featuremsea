#' Iteration-Aware Permutation Test for Iterative fMSEA (Experimental)
#'
#' Computes pathway p-values that account for the data-dependent iterative
#' re-weighting of the annotation matrix performed by
#' \code{\link{perform_fmsea_analysis}}.
#'
#' In \code{perform_fmsea_analysis}, the reported p-values come from a
#' permutation test that treats the annotation matrix of the final iteration as
#' fixed. That matrix, however, is derived from the observed ranking (FDR
#' selection + leading-edge re-weighting in earlier iterations), so the
#' reported p-values are conditional on a data-dependent matrix. This function
#' instead re-runs the \emph{entire} iterative procedure (starting from the
#' original annotation matrix A0) for every permuted ranking, and compares the
#' observed final-iteration ES of each pathway with the null distribution of
#' final-iteration ES values.
#'
#' The procedure for one run (observed or permuted) is identical to
#' \code{perform_fmsea_analysis}: the same pathway filtering, global
#' normalization, inner permutation test (used only for FDR-based selection),
#' re-weighting and convergence rule. With the same \code{seed} and
#' \code{inner.perm.num == perm.num}, the observed run reproduces the results of
#' \code{perform_fmsea_analysis}.
#'
#' The outer permutation shuffles \code{ranking_weight} across features, which
#' is the same feature-label null used by the inner test.
#'
#' @section Computational cost:
#' Every outer permutation costs about as much as one full
#' \code{perform_fmsea_analysis} run. Total cost is roughly
#' \code{(outer.perm.num + 1) x iterations x pathways x inner.perm.num} ES
#' evaluations. The smallest attainable iteration-aware p-value is
#' \code{1 / (outer.perm.num + 1)}.
#'
#' @section Choosing inner.perm.num:
#' Keep \code{inner.perm.num} equal to the \code{perm.num} used in
#' \code{perform_fmsea_analysis}. The inner p-values drive which pathways are
#' selected for re-weighting; if \code{inner.perm.num} is too small, the minimum
#' inner p-value \code{1 / (inner.perm.num + 1)} may never pass the FDR
#' threshold, re-weighting never happens, and the test no longer reflects the
#' procedure it is meant to evaluate.
#'
#' @param pathway_database An S4 pathway database object (same as in
#'   \code{perform_fmsea_analysis}).
#' @param annotation_table The original annotation score table (A0).
#' @param ranking_table Data frame with \code{variable_id} and
#'   \code{ranking_weight} columns.
#' @param outer.perm.num Integer. Number of outer permutations, each re-running
#'   the full iterative procedure. Default is 100.
#' @param inner.perm.num Integer. Number of inner permutations per pathway used
#'   for FDR-based selection within each iteration. Default is 1000.
#' @param threads Integer. Number of cores. The observed run is parallelized
#'   over pathways. Null runs are parallelized over outer permutations when
#'   \code{outer.perm.num >= threads}; otherwise they run one at a time, each
#'   parallelized over pathways. Results do not depend on this choice. Default
#'   is \code{NULL} (all cores minus one).
#' @param min.compounds.num,max.compounds.num Pathway size limits.
#' @param id.col Character. Compound ID column. Default is "KEGG_ID".
#' @param seed Integer. Random seed. Default is 123.
#' @param fdr.thr Numeric. FDR threshold used for selection inside the
#'   iterations. Default is 0.05.
#' @param max.iter.num Integer. Maximum number of iterations. Default is 3.
#' @param return.null Logical. Whether to return the null ES matrix
#'   (outer permutations x pathways). Default is \code{FALSE}.
#' @param verbose Logical. Whether to print progress messages. Default is \code{TRUE}.
#'
#' @return A list with:
#'   \item{results}{Data frame, one row per pathway, with columns
#'     \code{pathway_id}, \code{pathway_name}, \code{pathway_description},
#'     \code{ES}, \code{NES} (final iteration, as in \code{perform_fmsea_analysis}),
#'     \code{p_value_conditional}, \code{FDR_conditional} (conditional on the
#'     data-derived final matrix, as reported by \code{perform_fmsea_analysis}),
#'     \code{p_value_A0}, \code{FDR_A0} (first iteration, original matrix A0),
#'     \code{p_value_iter_aware}, \code{FDR_iter_aware}, \code{NES_iter_aware}
#'     (iteration-aware null), and \code{n_null_valid}.}
#'   \item{observed}{List with \code{iterations_used}, \code{converged},
#'     \code{n_sig_iter1}, \code{n_sig_final} and \code{res_list} of the
#'     observed run.}
#'   \item{null_diagnostics}{Data frame with, per outer permutation,
#'     \code{iterations_used}, \code{converged}, \code{n_sig_iter1} and
#'     \code{n_sig_final}.}
#'   \item{null_ES}{Matrix of null ES values (only if \code{return.null = TRUE}).}
#'   \item{params}{List of the parameters used.}
#'
#' @author Yijiang Liu \email{ejoliu@@outlook.com}
#'
#' @importFrom foreach %dopar%
#' @importFrom dplyr rename distinct arrange desc mutate filter %>%
#' @export
perform_fmsea_iterative_permutation <- function(
    pathway_database,
    annotation_table,
    ranking_table,
    outer.perm.num = 100,
    inner.perm.num = 1000,
    threads = NULL,
    min.compounds.num = 15,
    max.compounds.num = 300,
    id.col = "KEGG_ID",
    seed = 123,
    fdr.thr = 0.05,
    max.iter.num = 3,
    return.null = FALSE,
    verbose = TRUE
) {
  outer.perm.num <- as.integer(outer.perm.num)
  if (is.na(outer.perm.num) || outer.perm.num < 1L) {
    stop("outer.perm.num must be a positive integer.")
  }
  if (!all(c("variable_id", "ranking_weight") %in% colnames(ranking_table))) {
    stop("ranking_table must contain 'variable_id' and 'ranking_weight' columns.")
  }

  # --- Step 1: Pathway preparation (same filtering as the main function) ---
  pw <- .itperm_prepare_pathways(
    pathway_database,
    id_col        = id.col,
    min_compounds = min.compounds.num,
    max_compounds = max.compounds.num
  )
  n_pw <- length(pw$ids_list)
  if (n_pw == 0L) {
    stop("No pathways left after filtering by compound count.")
  }

  # --- Step 2: Ranking preparation (absolute weights, as in the main function) ---
  ranking_table_calc <- ranking_table %>%
    dplyr::mutate(ranking_weight = abs(ranking_weight)) %>%
    dplyr::distinct(variable_id, ranking_weight)

  if (anyDuplicated(ranking_table_calc$variable_id) > 0L) {
    warning("ranking_table contains variable_ids with multiple ranking_weight values; ",
            "the permutation shuffles rows, not unique features.")
  }

  # --- Step 3: Parallel configuration ---
  total_cores <- parallel::detectCores()
  if (is.null(threads) || threads >= total_cores) {
    threads <- max(1L, total_cores - 1L)
  }
  threads <- as.integer(threads)

  # Parallelize over outer permutations when there are enough of them to keep
  # all cores busy; otherwise run them one by one and parallelize over pathways
  # inside each run (as perform_fmsea_analysis does).
  outer_parallel <- threads > 1L && outer.perm.num >= threads

  doParallel::registerDoParallel(cores = threads)
  on.exit(doParallel::stopImplicitCluster(), add = TRUE)

  # --- Step 4: Observed run (pathway-level parallel) ---
  if (verbose) {
    message(sprintf("Running observed iterative fMSEA (%d pathways, %d core(s))...",
                    n_pw, threads))
  }
  t_obs <- system.time(
    obs <- .itperm_run_iterations(
      ids_list         = pw$ids_list,
      mfm_ids          = pw$mfm_ids,
      annotation_table = annotation_table,
      ranking_table    = ranking_table_calc,
      id_col           = id.col,
      n_perm           = inner.perm.num,
      seed             = seed,
      fdr_thr          = fdr.thr,
      max_iter         = max.iter.num,
      keep_first       = TRUE,
      n_cores          = threads,
      verbose          = verbose
    )
  )[["elapsed"]]

  ES_obs <- .itperm_extract(obs$res_list, "ES")

  # --- Step 5: Pre-generate outer permutations (independent of core count) ---
  set.seed(seed)
  n_rank <- nrow(ranking_table_calc)
  perm_idx <- lapply(seq_len(outer.perm.num), function(b) sample.int(n_rank))

  if (verbose) {
    message(sprintf(
      "Observed run: %d iteration(s), %d significant pathway(s) in iteration 1, %.1f s.",
      obs$iterations_used, obs$n_sig_iter1, t_obs
    ))
    if (outer_parallel) {
      # A serial run costs roughly threads x the parallel observed run.
      est <- t_obs * threads * ceiling(outer.perm.num / threads)
      mode_txt <- sprintf("parallel over permutations on %d core(s)", threads)
    } else {
      est <- t_obs * outer.perm.num
      mode_txt <- sprintf("one at a time, parallel over pathways on %d core(s)", threads)
    }
    message(sprintf(
      "Running %d outer permutations (%s); upper-bound estimate ~%.1f min.",
      outer.perm.num, mode_txt, est / 60
    ))
  }

  # --- Step 6: Null runs (full iterative procedure per permuted ranking) ---
  run_null <- function(b, n_cores) {
    rk_perm <- ranking_table_calc
    rk_perm$ranking_weight <- ranking_table_calc$ranking_weight[perm_idx[[b]]]

    run_b <- .itperm_run_iterations(
      ids_list         = pw$ids_list,
      mfm_ids          = pw$mfm_ids,
      annotation_table = annotation_table,
      ranking_table    = rk_perm,
      id_col           = id.col,
      n_perm           = inner.perm.num,
      seed             = seed + b * n_pw,
      fdr_thr          = fdr.thr,
      max_iter         = max.iter.num,
      keep_first       = FALSE,
      n_cores          = n_cores,
      verbose          = FALSE
    )

    list(
      ES              = .itperm_extract(run_b$res_list, "ES"),
      iterations_used = run_b$iterations_used,
      converged       = run_b$converged,
      n_sig_iter1     = run_b$n_sig_iter1,
      n_sig_final     = run_b$n_sig_final
    )
  }

  if (outer_parallel) {
    null_out <- foreach::foreach(
      b = seq_len(outer.perm.num),
      .packages = c("dplyr", "stringr"),
      .export   = .itperm_export_names()
    ) %dopar% {
      run_null(b, n_cores = 1L)
    }
  } else {
    null_out <- vector("list", outer.perm.num)
    for (b in seq_len(outer.perm.num)) {
      t_b <- system.time(null_out[[b]] <- run_null(b, n_cores = threads))[["elapsed"]]
      if (verbose) {
        message(sprintf("  Outer permutation %d/%d: %d iteration(s), %.1f s.",
                        b, outer.perm.num, null_out[[b]]$iterations_used, t_b))
      }
    }
  }

  null_ES <- do.call(rbind, lapply(null_out, `[[`, "ES"))
  null_ES <- matrix(null_ES, nrow = outer.perm.num, ncol = n_pw,
                    dimnames = list(NULL, pw$mfm_ids))

  null_diagnostics <- data.frame(
    perm            = seq_len(outer.perm.num),
    iterations_used = vapply(null_out, `[[`, integer(1), "iterations_used"),
    converged       = vapply(null_out, `[[`, logical(1), "converged"),
    n_sig_iter1     = vapply(null_out, `[[`, integer(1), "n_sig_iter1"),
    n_sig_final     = vapply(null_out, `[[`, integer(1), "n_sig_final"),
    stringsAsFactors = FALSE
  )

  # --- Step 7: Iteration-aware p-values ---
  n_null_valid <- colSums(!is.na(null_ES))
  p_iter <- vapply(seq_len(n_pw), function(k) {
    if (is.na(ES_obs[k]) || n_null_valid[k] == 0L) return(NA_real_)
    nk <- null_ES[, k]
    nk <- nk[!is.na(nk)]
    (sum(nk >= ES_obs[k]) + 1) / (length(nk) + 1)
  }, numeric(1))
  NES_iter <- vapply(seq_len(n_pw), function(k) {
    nk <- null_ES[, k]
    nk <- nk[!is.na(nk)]
    if (is.na(ES_obs[k]) || length(nk) == 0L) return(NA_real_)
    ES_obs[k] / mean(abs(nk))
  }, numeric(1))

  p_cond <- .itperm_extract(obs$res_list, "p_value")
  p_A0   <- .itperm_extract(obs$first_res_list, "p_value")

  results <- data.frame(
    pathway_id          = pw$mfm_ids,
    pathway_name        = pw$mfm_names,
    pathway_description = pw$mfm_descriptions,
    ES                  = ES_obs,
    NES                 = .itperm_extract(obs$res_list, "NES"),
    p_value_conditional = p_cond,
    FDR_conditional     = .itperm_bh(p_cond),
    p_value_A0          = p_A0,
    FDR_A0              = .itperm_bh(p_A0),
    p_value_iter_aware  = p_iter,
    FDR_iter_aware      = .itperm_bh(p_iter),
    NES_iter_aware      = NES_iter,
    n_null_valid        = as.integer(n_null_valid),
    stringsAsFactors    = FALSE
  )
  results <- results[order(results$p_value_iter_aware, -results$ES, na.last = TRUE), , drop = FALSE]
  rownames(results) <- NULL

  if (verbose) {
    message(sprintf(
      "Done. Re-weighting triggered in %d / %d null runs (>= 1 significant pathway in iteration 1).",
      sum(null_diagnostics$n_sig_iter1 > 0L), outer.perm.num
    ))
    message(sprintf(
      "Significant at FDR < %.2f: conditional = %d, A0 = %d, iteration-aware = %d (min attainable p = %.4g).",
      fdr.thr,
      sum(results$FDR_conditional < fdr.thr, na.rm = TRUE),
      sum(results$FDR_A0 < fdr.thr, na.rm = TRUE),
      sum(results$FDR_iter_aware < fdr.thr, na.rm = TRUE),
      1 / (outer.perm.num + 1)
    ))
  }

  out <- list(
    results = results,
    observed = list(
      iterations_used = obs$iterations_used,
      converged       = obs$converged,
      n_sig_iter1     = obs$n_sig_iter1,
      n_sig_final     = obs$n_sig_final,
      res_list        = obs$res_list
    ),
    null_diagnostics = null_diagnostics,
    params = list(
      outer.perm.num    = outer.perm.num,
      inner.perm.num    = inner.perm.num,
      min.compounds.num = min.compounds.num,
      max.compounds.num = max.compounds.num,
      id.col            = id.col,
      seed              = seed,
      fdr.thr           = fdr.thr,
      max.iter.num      = max.iter.num
    )
  )
  if (return.null) out[["null_ES"]] <- null_ES
  out
}

# ------------------------------------------------------------
# Internal helpers for perform_fmsea_iterative_permutation
# ------------------------------------------------------------

# Convert the pathway database and apply the same compound-count filter as
# parallel_computing_pathways_indexed_fast(). Pathway order (and therefore the
# per-pathway seeds seed + i) matches the main function.
.itperm_prepare_pathways <- function(pathway_database,
                                     id_col,
                                     min_compounds,
                                     max_compounds) {
  if (!isS4(pathway_database) || !"database_info" %in% methods::slotNames(pathway_database)) {
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

  if (!id_col %in% colnames(pathway_df)) {
    stop("pathway_database does not provide the id column: ", id_col)
  }

  id_str <- pathway_df[[id_col]]
  id_count <- ifelse(is.na(id_str), 0L, stringr::str_count(id_str, "\\{\\}") + 1L)
  pathway_df <- pathway_df[id_count >= min_compounds & id_count <= max_compounds, , drop = FALSE]

  list(
    mfm_ids          = as.character(pathway_df$pathway_id),
    mfm_names        = as.character(pathway_df$pathway_name),
    mfm_descriptions = as.character(pathway_df$pathway_description),
    ids_list         = lapply(pathway_df[[id_col]], split_ids)
  )
}

# One full iterative fMSEA run (serial over pathways), mirroring the loop in
# perform_fmsea_analysis(). Returns the res_list of the last iteration run,
# i.e. the one whose statistics perform_fmsea_analysis() reports.
.itperm_run_iterations <- function(ids_list,
                                   mfm_ids,
                                   annotation_table,
                                   ranking_table,
                                   id_col,
                                   n_perm,
                                   seed,
                                   fdr_thr,
                                   max_iter,
                                   keep_first = FALSE,
                                   n_cores = 1L,
                                   verbose = FALSE) {
  base_rk <- ranking_table %>%
    dplyr::distinct(variable_id, ranking_weight) %>%
    dplyr::arrange(dplyr::desc(ranking_weight))
  rk_idx <- build_ranking_index(base_rk)

  score_annotation_table <- annotation_table
  prev_count     <- NULL
  converged      <- FALSE
  iter_used      <- 0L
  res_list       <- NULL
  first_res_list <- NULL
  n_sig_iter1    <- 0L
  n_sig_final    <- 0L

  for (iter in seq_len(max_iter)) {
    iter_used <- iter

    annotation_long <- annotation_long_fast_base(
      score_annotation_table,
      id_col = id_col,
      normalize_global = TRUE
    )

    run_one <- function(i) {
      get_mMSEA_results_indexed_fast(
        pathway_ids_vec = ids_list[[i]],
        annotation_long = annotation_long,
        rk_idx          = rk_idx,
        id_col          = id_col,
        n_perm          = n_perm,
        seed            = seed + i,
        return_perm     = FALSE
      )
    }

    # Pathway-level parallelism uses the doParallel backend registered by the
    # caller; per-pathway seeds keep results identical to the serial path.
    res_list <- if (n_cores > 1L) {
      foreach::foreach(
        i = seq_along(ids_list),
        .packages = c("dplyr", "stringr"),
        .export   = .itperm_export_names()
      ) %dopar% {
        run_one(i)
      }
    } else {
      lapply(seq_along(ids_list), run_one)
    }
    names(res_list) <- mfm_ids

    if (iter == 1L && keep_first) first_res_list <- res_list

    significant_mfm <- get_significant_mfm(res_list, fdr_threshold = fdr_thr)
    n_sig_final <- nrow(significant_mfm)
    if (iter == 1L) n_sig_iter1 <- n_sig_final
    if (verbose) {
      message(sprintf("  Iteration %d: %d significant pathway(s).", iter, n_sig_final))
    }

    if (n_sig_final == 0L) break

    fm_count <- get_fm_long_table(
      significant_mfm,
      res_list,
      id_col = id_col,
      fdr_threshold = fdr_thr
    )
    weighted_table <- get_weighting_annotation_table_fast(
      annotation_table,
      fm_count,
      id_col = id_col
    )

    if (!is.null(prev_count) && counts_equal(prev_count, fm_count, id_col = id_col)) {
      converged <- TRUE
      break
    }

    prev_count <- fm_count
    score_annotation_table <- weighted_table
  }

  list(
    res_list        = res_list,
    first_res_list  = first_res_list,
    iterations_used = as.integer(iter_used),
    converged       = converged,
    n_sig_iter1     = as.integer(n_sig_iter1),
    n_sig_final     = as.integer(n_sig_final)
  )
}

# Internal functions needed on PSOCK workers (no-op for forked workers).
.itperm_export_names <- function() {
  c(
    ".itperm_run_iterations",
    ".itperm_extract",
    ".itperm_export_names",
    "annotation_long_fast_base",
    "build_ranking_index",
    "precompute_mMSEA_static",
    "compute_ES_from_mapping",
    "get_mMSEA_results_indexed_fast",
    "compute_permutations_cpp",
    "get_significant_mfm",
    "get_fm_long_table",
    "get_weighting_annotation_table_fast",
    "counts_equal",
    "canon_counts"
  )
}

# Extract a scalar field from every element of a res_list.
.itperm_extract <- function(res_list, field) {
  vapply(res_list, function(x) {
    v <- x[[field]]
    if (is.null(v) || length(v) == 0L) NA_real_ else as.numeric(v[[1]])
  }, numeric(1), USE.NAMES = FALSE)
}

# BH adjustment over finite p-values only (same convention as get_significant_mfm).
.itperm_bh <- function(p) {
  adj <- rep(NA_real_, length(p))
  idx <- is.finite(p)
  if (any(idx)) adj[idx] <- stats::p.adjust(p[idx], method = "BH")
  adj
}
