#' @title Analyze Matrix Confidence for Metabolic Pathways
#' @description
#' Evaluate the confidence of inferring metabolic pathway activity from
#' metabolites detected in a given sample matrix (e.g., urine, plasma, feces) using GPT-4 or Qwen.
#'
#' @param results Object containing significant_modules slot with pathway information,
#'   or a data.frame with required columns (pathway_id, pathway_name, pathway_description)
#' @param sample_source Sample matrix type. Options: "urine", "plasma", "serum", "blood", "feces"
#' @param api_key API key (OpenAI or SiliconFlow, depending on provider)
#' @param provider API provider. Required. Options: "openai", "siliconflow"
#' @param model Model to use. If NULL, uses provider default: "gpt-4o" for OpenAI, "Qwen3-32B" for SiliconFlow
#' @param temperature Sampling temperature for model. Default: 0.2
#' @param max_tokens Maximum tokens in response. Default: 8000
#'
#' @return A data frame based on significant_modules with three additional columns:
#' \itemize{
#'   \item \code{matrix_confidence_score}: Integer score (0/25/50/75/100)
#'   \item \code{matrix_confidence_reason}: Brief explanation of the score
#'   \item \code{matrix_source}: The sample source used for evaluation
#' }
#'
#' @details
#' The function uses LLM to assess the confidence of inferring pathway activity from
#' metabolites in the given matrix with the following scoring:
#' \itemize{
#'   \item 100: High confidence - metabolites in this matrix are reliable indicators of pathway activity
#'   \item 75: Good confidence - metabolites likely reflect pathway activity with some certainty
#'   \item 50: Moderate confidence - metabolites may indicate pathway activity but with uncertainty
#'   \item 25: Low confidence - metabolites provide weak evidence of pathway activity
#'   \item 0: No confidence - metabolites in this matrix do not reliably indicate pathway activity
#' }
#'
#' @examples
#' \dontrun{
#' # Using OpenAI
#' result_df <- analyze_matrix_relevance(
#'   results = my_results,
#'   sample_source = "urine",
#'   api_key = "sk-openai-key",
#'   provider = "openai"
#' )
#'
#' # Using SiliconFlow with Qwen
#' result_df <- analyze_matrix_relevance(
#'   results = my_results,
#'   sample_source = "urine",
#'   api_key = "sk-siliconflow-key",
#'   provider = "siliconflow"
#' )
#'
#' # Custom model specification
#' result_df <- analyze_matrix_relevance(
#'   results = my_results,
#'   sample_source = "plasma",
#'   api_key = "sk-siliconflow-key",
#'   provider = "siliconflow",
#'   model = "Qwen3-72B"
#' )
#' }
#'
#' @importFrom httr2 request req_auth_bearer_token req_headers req_body_json req_perform req_timeout resp_body_json
#' @importFrom jsonlite toJSON fromJSON
#' @importFrom dplyr mutate transmute select distinct left_join %>%
#' @importFrom tibble as_tibble
#' @importFrom stringr str_squish
#' @importFrom glue glue
#' @importFrom methods slotNames is
#'
#' @export
analyze_matrix_relevance <- function(results,
                                     sample_source,
                                     api_key,
                                     provider,
                                     model = NULL,
                                     temperature = 0.2,
                                     max_tokens = 8000) {

  # -------- Validate API Key --------
  if (missing(api_key) || is.null(api_key) || nchar(api_key) == 0) {
    stop("api_key is required. Please provide a valid API key.")
  }

  # -------- Validate Provider --------
  if (missing(provider) || is.null(provider) || nchar(trimws(provider)) == 0) {
    stop("provider is required. Please choose one of: 'openai' or 'siliconflow'")
  }

  provider <- tolower(trimws(provider))
  valid_providers <- c("openai", "siliconflow")
  if (!provider %in% valid_providers) {
    stop("Invalid provider '", provider, "'. Please choose one of: 'openai' or 'siliconflow'")
  }

  # -------- Configure API based on provider --------
  if (provider == "openai") {
    api_base <- "https://api.openai.com/v1"
    default_model <- "gpt-4o"
  } else if (provider == "siliconflow") {
    api_base <- "https://api.siliconflow.cn/v1"
    default_model <- "Qwen3-32B"
  }

  # Use provided model or default
  if (is.null(model)) {
    model <- default_model
  }
  
  # -------- Input Validation --------
  valid_sources <- c("urine", "plasma", "serum", "blood", "feces")
  if (!tolower(sample_source) %in% valid_sources) {
    warning(
      "sample_source '", sample_source, "' is not in the standard list: ",
      paste(valid_sources, collapse = ", "),
      ". Proceeding anyway."
    )
  }
  sample_source <- tolower(sample_source)
  
  # -------- Validate results is S4 with significant_modules --------
  if (!isS4(results) || !"significant_modules" %in% methods::slotNames(results)) {
    stop(
      "results must be an S4 object with a 'significant_modules' slot. ",
      "If you have a data.frame, wrap it in the appropriate S4 object first."
    )
  }
  
  req_cols <- c("pathway_id", "pathway_name", "pathway_description")
  significant_modules <- methods::slot(results, "significant_modules")
  if (!all(req_cols %in% names(significant_modules))) {
    stop("significant_modules must contain columns: ", paste(req_cols, collapse = ", "))
  }
  
  # -------- Internal: Build System Prompt --------
  .build_system_prompt <- function() {
    glue::glue(
      "You're an expert in LC-MS untargeted metabolomics pathway interpretation.
Task: For each pathway, assess the CONFIDENCE of inferring pathway activity from metabolites detected in the given Sample Source (matrix). Note that pathways do not exist directly in biological matrices, but their activity can be inferred from the presence and levels of related metabolites.
Use DISCRETE ANCHORS:
matrix_confidence (0/25/50/75/100):
- 100: High confidence - metabolites in this matrix reliably indicate pathway activity.
- 75: Good confidence - metabolites likely reflect pathway activity with reasonable certainty.
- 50: Moderate confidence - metabolites may indicate pathway activity but interpretation requires caution.
- 25: Low confidence - metabolites provide weak evidence of pathway activity due to confounding factors.
- 0: No confidence - metabolites in this matrix are unreliable indicators of pathway activity.
Matrix interpretation guidelines (STRICT CONSTRAINTS):
- Urine: Reflects excretion and filtration products; good for catabolic end-products.
  * STRICT CONSTRAINT: **Score matrix_confidence strictly as 0** for
    'Aminoacyl-tRNA biosynthesis' and 'Amino acid biosynthesis'.
    Amino acids in urine typically result from filtration/pathological leaks
    rather than reflecting active biosynthesis, making inference unreliable.
- Feces: Dominated by microbial metabolism and dietary remnants.
  * STRICT CONSTRAINT: **Score matrix_confidence strictly as 0** for
    'Lipid biosynthesis' and host 'Amino acid metabolism'.
    Fecal metabolites are primarily microbial or dietary in origin,
    not reliable indicators of host metabolic pathway activity.
- Plasma/Serum/Blood: Reflects systemic circulation and active metabolism.
  * STRICT CONSTRAINT: **Score matrix_confidence strictly as 0** for
    'Amino acid biosynthesis'. Blood amino acids reflect consumption/turnover
    rather than active biosynthesis, making inference of biosynthetic
    pathway activity unreliable.
Output STRICTLY compact JSON:
{{
  \"results\": [
    {{
      \"pathway_id\": string,
      \"pathway_name\": string,
      \"matrix_confidence\": integer (0,25,50,75/100),
      \"reason\": string (<= 2 sentences; explain inference confidence)
    }},
    ...
  ]
}}
"
    )
  }
  
  # -------- Internal: Build User Prompt --------
  .build_user_prompt <- function(df, sample_source) {
    items <- df %>%
      dplyr::transmute(
        pathway_id = as.character(pathway_id),
        pathway_name = as.character(pathway_name),
        pathway_description = stringr::str_squish(as.character(pathway_description))
      ) %>%
      dplyr::mutate(pathway_description = substr(pathway_description, 1, 800))
    
    list(
      sample_source = sample_source,
      pathways = items
    ) %>%
      jsonlite::toJSON(auto_unbox = TRUE, pretty = FALSE)
  }
  
  # -------- Internal: Call LLM API --------
  .call_llm_api <- function(messages, model, temperature, max_tokens, api_key, api_base) {
    
    req <- httr2::request(paste0(api_base, "/chat/completions")) |>
      httr2::req_auth_bearer_token(api_key) |>
      httr2::req_headers("Content-Type" = "application/json") |>
      httr2::req_timeout(seconds = 300)
    
    body <- list(
      model = model,
      temperature = temperature,
      max_tokens = max_tokens,
      messages = messages,
      response_format = list(type = "json_object")
    )
    
    resp <- tryCatch(
      {
        req |>
          httr2::req_body_json(body) |>
          httr2::req_perform()
      },
      error = function(e) {
        stop("OpenAI API request failed: ", conditionMessage(e))
      }
    )
    
    return(resp)
  }
  
  # -------- Internal: Parse Response and Merge --------
  .parse_and_merge <- function(txt, original_df, sample_source) {
    out <- tryCatch(jsonlite::fromJSON(txt), error = function(e) NULL)
    
    if (is.null(out) || is.null(out$results)) {
      warning("Model response could not be parsed as expected JSON.")
      return(
        original_df %>%
          dplyr::mutate(
            matrix_confidence_score = NA_integer_,
            matrix_confidence_reason = "Model response could not be parsed as expected JSON.",
            matrix_source = sample_source
          )
      )
    }
    
    res <- tibble::as_tibble(out$results) %>%
      dplyr::mutate(
        pathway_id = as.character(pathway_id),
        matrix_confidence_score = suppressWarnings(as.integer(matrix_confidence)),
        matrix_confidence_reason = as.character(reason)
      ) %>%
      dplyr::select(pathway_id, matrix_confidence_score, matrix_confidence_reason)
    
    clamp_to_set <- function(x) {
      allowed <- c(0L, 25L, 50L, 75L, 100L)
      ifelse(x %in% allowed, x, NA_integer_)
    }
    
    res <- res %>%
      dplyr::mutate(matrix_confidence_score = clamp_to_set(matrix_confidence_score))
    
    original_df %>%
      dplyr::mutate(pathway_id = as.character(pathway_id)) %>%
      dplyr::left_join(res, by = "pathway_id") %>%
      dplyr::mutate(matrix_source = sample_source)
  }
  
  # -------- Main Logic --------
  
  # Prepare data
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
        matrix_confidence_score = NA_integer_,
        matrix_confidence_reason = "No pathways to analyze.",
        matrix_source = sample_source
      )
    methods::slot(results, "significant_modules") <- updated_modules
    return(results)
  }
  
  message("Analyzing ", nrow(df), " unique pathways for matrix confidence in '", sample_source, "'...")
  
  # Build prompts
  system_prompt <- .build_system_prompt()
  user_payload_json <- .build_user_prompt(df, sample_source)
  
  messages <- list(
    list(role = "system", content = system_prompt),
    list(
      role = "user",
      content = paste0(
        "Sample Source: ", sample_source, "\n",
        "Evaluate the matrix confidence for inferring pathway activity from the following pathways (as JSON). Each reason <= 2 sentences.\n",
        user_payload_json
      )
    )
  )
  
  # Call API
  resp <- .call_llm_api(
    messages = messages,
    model = model,
    temperature = temperature,
    max_tokens = max_tokens,
    api_key = api_key,
    api_base = api_base
  )
  
  # Parse response
  resp_json <- tryCatch(
    httr2::resp_body_json(resp),
    error = function(e) {
      warning("Failed to parse API response: ", conditionMessage(e))
      NULL
    }
  )
  
  content <- if (!is.null(resp_json) && length(resp_json$choices) > 0) {
    resp_json$choices[[1]]$message$content
  } else {
    ""
  }
  
  updated_modules <- .parse_and_merge(content, original_df = significant_modules, sample_source = sample_source)
  
  # Update the S4 object's significant_modules slot
  methods::slot(results, "significant_modules") <- updated_modules
  
  message("Analysis complete. ", sum(!is.na(updated_modules$matrix_confidence_score)), " pathways scored.")
  
  return(results)
}