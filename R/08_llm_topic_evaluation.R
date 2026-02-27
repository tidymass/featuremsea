#' @title Analyze Literature Relevance for Metabolic Pathways
#' @description
#' Evaluate the literature support linking significant metabolic pathways to a
#' specific research topic (e.g., cancer, pregnancy, diabetes) by querying PubMed.
#' For each pathway, the function searches PubMed for co-occurrence with the
#' research topic and returns supporting PMIDs. If no literature is found,
#' GPT-4 provides a plausible biological explanation for the association.
#'
#' @param results Object containing significant_modules slot with pathway information,
#'   must be an S4 object with required columns
#'   (pathway_id, pathway_name, pathway_description) in `significant_modules`.
#' @param research_topic Character string describing the research topic (e.g., "cancer",
#'   "pregnancy", "type 2 diabetes", "Parkinson's disease")
#' @param api_key OpenAI API key (required for biological explanations when no literature found)
#' @param pubmed_api_key Optional NCBI API key for higher rate limits on PubMed queries.
#'   Without it, requests are limited to 3/second; with it, 10/second.
#' @param max_pmids Deprecated and kept for backward compatibility. Fuzzy matching
#'   PMIDs are always truncated to top 10.
#' @param model OpenAI model to use for biological explanations. Default: "gpt-4o"
#' @param temperature Sampling temperature for model. Default: 0.2
#' @param max_tokens Maximum tokens in response. Default: 8000
#' @param rate_delay Delay in seconds between PubMed requests to avoid rate limiting.
#'   Default: 0.35
#'
#' @return Updated input S4 object with `significant_modules` augmented by:
#' \itemize{
#'   \item \code{literature_pmids_exact}: PMIDs from exact query, separated by \code{{}}
#'   \item \code{literature_pmids_fuzzy}: PMIDs from fuzzy query, separated by \code{{}} (top 10 only)
#'   \item \code{biological_explanation}: AI-generated explanation when no literature found, NA otherwise
#'   \item \code{research_topic}: The research topic used for evaluation
#' }
#'
#' @examples
#' \dontrun{
#' result_obj <- analyze_literature_relevance(
#'   results = my_results,
#'   research_topic = "breast cancer",
#'   api_key = "your-openai-api-key"
#' )
#' }
#'
#' @importFrom httr2 request req_url_query req_headers req_perform req_timeout
#'   resp_body_json req_auth_bearer_token req_body_json
#' @importFrom jsonlite toJSON fromJSON
#' @importFrom dplyr mutate transmute select distinct left_join %>% bind_rows filter pull
#' @importFrom tibble as_tibble tibble
#' @importFrom stringr str_squish str_replace_all
#' @importFrom glue glue
#' @importFrom methods slotNames is
#'
#' @export
analyze_literature_relevance <- function(results,
                                         research_topic,
                                         api_key,
                                         pubmed_api_key = NULL,
                                         max_pmids = 10,
                                         model = "gpt-4.1",
                                         temperature = 0.2,
                                         max_tokens = 8000,
                                         rate_delay = 0.35) {
  
  # -------- Validate Inputs --------
  if (missing(api_key) || is.null(api_key) || nchar(api_key) == 0) {
    stop("api_key is required. Please provide a valid OpenAI API key.")
  }
  if (missing(research_topic) || is.null(research_topic) || nchar(trimws(research_topic)) == 0) {
    stop("research_topic is required. Please provide a research topic (e.g., 'cancer', 'pregnancy').")
  }
  research_topic <- trimws(research_topic)
  
  # -------- Validate results is S4 with significant_modules --------
  if (!isS4(results) || !"significant_modules" %in% methods::slotNames(results)) {
    stop(
      "results must be an S4 object with a 'significant_modules' slot."
    )
  }
  
  req_cols <- c("pathway_id", "pathway_name", "pathway_description")
  significant_modules <- methods::slot(results, "significant_modules")
  if (!all(req_cols %in% names(significant_modules))) {
    stop("significant_modules must contain columns: ", paste(req_cols, collapse = ", "))
  }
  
  # -------- Internal: Build PubMed Search Query --------
  .build_pubmed_query <- function(pathway_name, research_topic) {
    clean_name <- pathway_name
    clean_name <- stringr::str_replace_all(clean_name, "\\s*\\(.*?\\)\\s*", " ")
    clean_name <- stringr::str_replace_all(clean_name, ",", "")
    clean_name <- stringr::str_squish(clean_name)
    query <- glue::glue('("{clean_name}"[Title/Abstract]) AND ("{research_topic}"[Title/Abstract])')
    return(as.character(query))
  }
  
  # -------- Internal: Search PubMed via E-utilities --------
  .search_pubmed <- function(query, retmax, retstart = 0L, pubmed_api_key) {
    base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    
    req <- httr2::request(base_url) |>
      httr2::req_url_query(
        db = "pubmed",
        term = query,
        retmax = as.integer(retmax),
        retstart = as.integer(retstart),
        retmode = "json",
        sort = "relevance"
      ) |>
      httr2::req_timeout(seconds = 30)
    
    if (!is.null(pubmed_api_key) && nchar(pubmed_api_key) > 0) {
      req <- req |> httr2::req_url_query(api_key = pubmed_api_key)
    }
    
    resp <- tryCatch(
      httr2::req_perform(req),
      error = function(e) {
        warning("PubMed search failed for query: ", query, "\n  Error: ", conditionMessage(e))
        return(NULL)
      }
    )
    
    if (is.null(resp)) return(list(count = 0L, pmids = character(0)))
    
    body <- tryCatch(httr2::resp_body_json(resp), error = function(e) NULL)
    
    if (is.null(body) || is.null(body$esearchresult)) {
      return(list(count = 0L, pmids = character(0)))
    }
    
    count <- as.integer(body$esearchresult$count %||% 0)
    pmids <- unlist(body$esearchresult$idlist) %||% character(0)
    
    return(list(count = count, pmids = as.character(pmids)))
  }
  
  # -------- Internal: Exact Search (retrieve all PMIDs with pagination) --------
  .search_pubmed_exact_all <- function(query, pubmed_api_key, page_size = 10000L) {
    first <- .search_pubmed(query = query, retmax = 0L, retstart = 0L, pubmed_api_key = pubmed_api_key)
    total <- as.integer(first$count %||% 0L)
    if (!is.finite(total) || total <= 0L) {
      return(list(count = 0L, pmids = character(0)))
    }
    
    ids <- character(0)
    starts <- seq.int(0L, total - 1L, by = as.integer(page_size))
    for (st in starts) {
      take_n <- min(as.integer(page_size), total - st)
      chunk <- .search_pubmed(query = query, retmax = take_n, retstart = st, pubmed_api_key = pubmed_api_key)
      ids <- c(ids, chunk$pmids)
    }
    
    ids <- unique(as.character(ids))
    list(count = length(ids), pmids = ids)
  }
  
  # -------- Internal: Broader Fallback Search --------
  .search_pubmed_broad <- function(pathway_name, research_topic, pubmed_api_key) {
    clean_name <- stringr::str_replace_all(pathway_name, "\\s*\\(.*?\\)\\s*", " ")
    clean_name <- stringr::str_replace_all(clean_name, ",", "")
    clean_name <- stringr::str_squish(clean_name)
    
    terms <- unlist(strsplit(clean_name, "\\s+"))
    terms <- terms[nchar(terms) > 3]
    if (length(terms) > 3) terms <- terms[1:3]
    if (length(terms) == 0) return(list(count = 0, pmids = character(0)))
    
    broad_query <- glue::glue(
      '({paste(terms, collapse = " AND ")}) AND ("{research_topic}"[Title/Abstract])'
    )
    out <- .search_pubmed(as.character(broad_query), retmax = 10L, retstart = 0L, pubmed_api_key = pubmed_api_key)
    out$pmids <- head(out$pmids, 10)
    out$count <- length(out$pmids)
    out
  }
  
  # -------- Internal: Get Biological Explanations via OpenAI --------
  .get_biological_explanations <- function(pathways_without_lit, research_topic, api_key,
                                           model, temperature, max_tokens) {
    if (nrow(pathways_without_lit) == 0) {
      return(tibble::tibble(pathway_id = character(0), biological_explanation = character(0)))
    }
    
    system_prompt <- glue::glue(
      "You are an expert in metabolomics, biochemistry, and molecular biology.

Task: For each metabolic pathway listed below, provide a concise biological explanation
for why this pathway *might* (or might not) be relevant to the research topic: \"{research_topic}\".

Consider:
- Known biochemical links between the pathway and the disease/condition
- Potential upstream/downstream regulatory connections
- Shared metabolites or cofactors
- Published hypotheses or emerging evidence

If the pathway is likely irrelevant to the research topic, say so clearly and explain why.

Output STRICTLY compact JSON:
{{
  \"explanations\": [
    {{
      \"pathway_id\": string,
      \"biological_explanation\": string (<= 3 sentences)
    }},
    ...
  ]
}}"
    )
    
    pathway_info <- pathways_without_lit %>%
      dplyr::transmute(
        pathway_id = as.character(pathway_id),
        pathway_name = as.character(pathway_name),
        pathway_description = stringr::str_squish(as.character(pathway_description))
      )
    
    user_content <- paste0(
      "Research topic: ", research_topic, "\n",
      "Provide biological explanations for the following pathways (no PubMed literature found):\n",
      jsonlite::toJSON(pathway_info, auto_unbox = TRUE, pretty = FALSE)
    )
    
    messages <- list(
      list(role = "system", content = system_prompt),
      list(role = "user", content = user_content)
    )
    
    api_base <- "https://api.openai.com/v1"
    req <- httr2::request(paste0(api_base, "/chat/completions")) |>
      httr2::req_auth_bearer_token(api_key) |>
      httr2::req_headers(`Content-Type` = "application/json") |>
      httr2::req_timeout(seconds = 300)
    
    body <- list(
      model = model,
      temperature = temperature,
      max_tokens = max_tokens,
      messages = messages,
      response_format = list(type = "json_object")
    )
    
    resp <- tryCatch(
      req |> httr2::req_body_json(body) |> httr2::req_perform(),
      error = function(e) {
        warning("OpenAI API request failed: ", conditionMessage(e))
        return(NULL)
      }
    )
    
    if (is.null(resp)) {
      return(tibble::tibble(
        pathway_id = pathways_without_lit$pathway_id,
        biological_explanation = rep("Failed to generate explanation (API error).",
                                     nrow(pathways_without_lit))
      ))
    }
    
    resp_json <- tryCatch(httr2::resp_body_json(resp), error = function(e) NULL)
    content <- if (!is.null(resp_json) && length(resp_json$choices) > 0) {
      resp_json$choices[[1]]$message$content
    } else {
      ""
    }
    
    parsed <- tryCatch(jsonlite::fromJSON(content), error = function(e) NULL)
    
    if (is.null(parsed) || is.null(parsed$explanations)) {
      return(tibble::tibble(
        pathway_id = pathways_without_lit$pathway_id,
        biological_explanation = rep("Failed to parse explanation from model.",
                                     nrow(pathways_without_lit))
      ))
    }
    
    tibble::as_tibble(parsed$explanations) %>%
      dplyr::transmute(
        pathway_id = as.character(pathway_id),
        biological_explanation = as.character(biological_explanation)
      )
  }
  
  # -------- Main Logic --------
  
  df <- significant_modules %>%
    dplyr::distinct(pathway_id, .keep_all = TRUE) %>%
    dplyr::mutate(
      pathway_id = as.character(pathway_id),
      pathway_name = as.character(pathway_name),
      pathway_description = as.character(pathway_description)
    )
  
  if (nrow(df) == 0) {
    warning("No pathways found in input data.")
    updated_modules <- significant_modules %>%
      dplyr::mutate(
        literature_pmids_exact = NA_character_,
        literature_pmids_fuzzy = NA_character_,
        biological_explanation = NA_character_,
        research_topic = research_topic
      )
    methods::slot(results, "significant_modules") <- updated_modules
    return(results)
  }
  
  message("Searching PubMed for ", nrow(df), " pathways related to '", research_topic, "'...")
  
  # -------- PubMed Search Loop --------
  pubmed_results <- list()
  
  for (i in seq_len(nrow(df))) {
    pathway_id <- df$pathway_id[i]
    pathway_name <- df$pathway_name[i]
    
    message("  [", i, "/", nrow(df), "] Searching: ", pathway_name)
    
    # Exact search (keep all matched PMIDs)
    query <- .build_pubmed_query(pathway_name, research_topic)
    exact_result <- .search_pubmed_exact_all(query, pubmed_api_key)
    
    # Fuzzy search only when exact search has no hit (top 10 only)
    if (length(exact_result$pmids) == 0) {
      message("    -> No exact-match results. Trying fuzzy search...")
      fuzzy_result <- .search_pubmed_broad(pathway_name, research_topic, pubmed_api_key)
    } else {
      fuzzy_result <- list(count = 0L, pmids = character(0))
    }
    
    # Format PMIDs with {} separator, or NA
    pubmed_results[[pathway_id]] <- tibble::tibble(
      pathway_id = pathway_id,
      literature_pmids_exact = if (length(exact_result$pmids) > 0) {
        paste(exact_result$pmids, collapse = "{}")
      } else {
        NA_character_
      },
      literature_pmids_fuzzy = if (length(fuzzy_result$pmids) > 0) {
        paste(head(fuzzy_result$pmids, 10), collapse = "{}")
      } else {
        NA_character_
      }
    )
    
    Sys.sleep(rate_delay)
  }
  
  pubmed_df <- dplyr::bind_rows(pubmed_results)
  
  # -------- Biological Explanations for Pathways Without Literature --------
  no_literature_ids <- pubmed_df %>%
    dplyr::filter(is.na(literature_pmids_exact) & is.na(literature_pmids_fuzzy)) %>%
    dplyr::pull(pathway_id)
  
  pathways_needing_explanation <- df %>%
    dplyr::filter(pathway_id %in% no_literature_ids)
  
  if (nrow(pathways_needing_explanation) > 0) {
    message("Generating biological explanations for ", nrow(pathways_needing_explanation),
            " pathways without literature support...")
    
    explanations <- .get_biological_explanations(
      pathways_without_lit = pathways_needing_explanation,
      research_topic = research_topic,
      api_key = api_key,
      model = model,
      temperature = temperature,
      max_tokens = max_tokens
    )
  } else {
    explanations <- tibble::tibble(pathway_id = character(0), biological_explanation = character(0))
  }
  
  # -------- Merge Results --------
  updated_modules <- significant_modules %>%
    dplyr::mutate(pathway_id = as.character(pathway_id)) %>%
    dplyr::left_join(pubmed_df, by = "pathway_id") %>%
    dplyr::left_join(explanations, by = "pathway_id") %>%
    dplyr::mutate(
      research_topic = research_topic,
      biological_explanation = ifelse(
        !is.na(literature_pmids_exact) | !is.na(literature_pmids_fuzzy),
        NA_character_,
        biological_explanation
      )
    )
  
  n_with_lit <- sum(!is.na(updated_modules$literature_pmids_exact) | !is.na(updated_modules$literature_pmids_fuzzy))
  n_without <- sum(is.na(updated_modules$literature_pmids_exact) & is.na(updated_modules$literature_pmids_fuzzy))
  message("Analysis complete. ", n_with_lit, " pathways with literature support, ", n_without, " with biological explanations.")
  
  methods::slot(results, "significant_modules") <- updated_modules
  return(results)
}
