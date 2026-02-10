#' @title Analyze Matrix Relevance for Metabolic Pathways
#' @description
#' Evaluate the detectability and biological meaningfulness of metabolic pathways
#' in a given sample matrix (e.g., urine, plasma, feces) using GPT-4.
#'
#' @param results Object containing significant_modules slot with pathway information,
#'   or a data.frame with required columns (pathway_id, pathway_name, pathway_description)
#' @param sample_source Sample matrix type. Options: "urine", "plasma", "serum", "blood", "feces"
#' @param api_key OpenAI API key (required)
#' @param model OpenAI model to use. Default: "gpt-4o"
#' @param temperature Sampling temperature for model. Default: 0.2
#' @param max_tokens Maximum tokens in response. Default: 8000
#'
#' @return A data frame based on significant_modules with three additional columns:
#' \itemize{
#'   \item \code{matrix_relevance_score}: Integer score (0/25/50/75/100)
#'   \item \code{matrix_relevance_reason}: Brief explanation of the score
#'   \item \code{matrix_source}: The sample source used for evaluation
#' }
#'
#' @details
#' The function uses GPT-4 to assess pathway detectability with the following scoring:
#' \itemize{
#'   \item 100: Strong, reliable evidence expected in this matrix
#'   \item 75: Likely detectable; several expected analytes
#'   \item 50: Possibly detectable; proxies exist but uncertain
#'   \item 25: Unlikely to be a valid readout
#'   \item 0: Not observable or completely misleading in this matrix
#' }
#'
#' @examples
#' \dontrun{
#' result_df <- analyze_matrix_relevance(
#'   results = my_results,
#'   sample_source = "urine",
#'   api_key = "your-api-key"
#' )
#'
#' # Using a data.frame directly
#' df <- data.frame(
#'   pathway_id = c("path1", "path2"),
#'   pathway_name = c("Glycolysis", "TCA Cycle"),
#'   pathway_description = c("Glucose metabolism", "Central carbon metabolism")
#' )
#' result_df <- analyze_matrix_relevance(
#'   results = df,
#'   sample_source = "plasma",
#'   api_key = "your-api-key"
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
                                     model = "gpt-4.1",
                                     temperature = 0.2,
                                     max_tokens = 8000) {
  
  # -------- Validate API Key --------
  if (missing(api_key) || is.null(api_key) || nchar(api_key) == 0) {
    stop("api_key is required. Please provide a valid OpenAI API key.")
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
  
  # -------- Internal: Build System Prompt --------
  .build_system_prompt <- function() {
    glue::glue(
      "You're an expert in LC-MS untargeted metabolomics pathway interpretation.

Task: For each pathway, assess the likelihood that its small-molecule evidence is detectable and *biologically meaningful* in the given Sample Source (matrix).

Use DISCRETE ANCHORS:
matrix_detectability (0/25/50/75/100):
- 100: Strong, reliable evidence expected in this matrix.
- 75: Likely detectable; several expected analytes.
- 50: Possibly detectable; proxies exist but uncertain.
- 25: Unlikely to be a valid readout (e.g., metabolites are artifacts, dietary waste, or pathological leaks).
- 0: Essentially not observable or completely misleading in this matrix.

Matrix heuristics (STRICT CONSTRAINTS):
- Urine: Enriched for small, excreted metabolites.
  * STRICT CONSTRAINT: **Score matrix_detectability strictly as 0** for 
    'Aminoacyl-tRNA biosynthesis' and 'Amino acid biosynthesis'. 
    Detection in urine implies filtration/leak rather than systemic 
    biosynthesis, making it a misleading marker for this pathway.

- Feces: Microbiome metabolism, SCFAs, bile acids.
  * STRICT CONSTRAINT: **Score matrix_detectability strictly as 0** for 
    'Lipid biosynthesis' and 'Amino acid metabolism' if the topic implies 
    host physiology. Fecal levels are dominated by microbial activity or 
    diet, rendering them invalid as host pathway readouts.

- Plasma/Serum/Blood: Systemic metabolism, generally high detectability.
  * STRICT CONSTRAINT: **Score matrix_detectability strictly as 0** for 
    'Amino acid biosynthesis'. Detection in blood implies consumption/
    degradation rather than active biosynthesis, making it a misleading 
    marker for this pathway.

Output STRICTLY compact JSON:
{{
  \"results\": [
    {{
      \"pathway_id\": string,
      \"pathway_name\": string,
      \"matrix_detectability\": integer (0,25,50,75/100),
      \"reason\": string (<= 2 sentences; explain matrix validity)
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
  
  # -------- Internal: Call OpenAI API --------
  .call_openai <- function(messages, model, temperature, max_tokens, api_key) {
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
            matrix_relevance_score = NA_integer_,
            matrix_relevance_reason = "Model response could not be parsed as expected JSON.",
            matrix_source = sample_source
          )
      )
    }
    
    res <- tibble::as_tibble(out$results) %>%
      dplyr::mutate(
        pathway_id = as.character(pathway_id),
        matrix_relevance_score = suppressWarnings(as.integer(matrix_detectability)),
        matrix_relevance_reason = as.character(reason)
      ) %>%
      dplyr::select(pathway_id, matrix_relevance_score, matrix_relevance_reason)
    
    # Clamp to allowed discrete values
    clamp_to_set <- function(x) {
      allowed <- c(0L, 25L, 50L, 75L, 100L)
      ifelse(x %in% allowed, x, NA_integer_)
    }
    
    res <- res %>%
      dplyr::mutate(matrix_relevance_score = clamp_to_set(matrix_relevance_score))
    
    # Join with original data and add matrix_source
    original_df %>%
      dplyr::mutate(pathway_id = as.character(pathway_id)) %>%
      dplyr::left_join(res, by = "pathway_id") %>%
      dplyr::mutate(matrix_source = sample_source)
  }
  
  # -------- Extract Data from Input --------
  .extract_data <- function(results) {
    req_cols <- c("pathway_id", "pathway_name", "pathway_description")
    
    # Case 1: S4 object with significant_modules slot (use isS4() for detection)
    if (isS4(results) && "significant_modules" %in% methods::slotNames(results)) {
      df <- methods::slot(results, "significant_modules")
      if (!all(req_cols %in% names(df))) {
        stop("significant_modules must contain columns: ", paste(req_cols, collapse = ", "))
      }
      return(df)
    }
    
    # Case 2: Direct data.frame input
    if (is.data.frame(results)) {
      if (!all(req_cols %in% names(results))) {
        stop("Input data.frame must contain columns: ", paste(req_cols, collapse = ", "))
      }
      return(results)
    }
    
    # Case 3: List with significant_modules element
    if (is.list(results) && "significant_modules" %in% names(results)) {
      df <- results$significant_modules
      if (!all(req_cols %in% names(df))) {
        stop("significant_modules must contain columns: ", paste(req_cols, collapse = ", "))
      }
      return(df)
    }
    
    stop(
      "results must be one of:\n",
      "  1. An S4 object with 'significant_modules' slot\n",
      "  2. A data.frame with columns: ", paste(req_cols, collapse = ", "), "\n",
      "  3. A list with 'significant_modules' element"
    )
  }
  
  # -------- Main Logic --------
  
  # Extract significant_modules from input
  significant_modules <- .extract_data(results)
  
  # Prepare data
  df <- significant_modules %>%
    dplyr::distinct(pathway_id, .keep_all = TRUE) %>%
    dplyr::mutate(
      pathway_id = as.character(pathway_id),
      pathway_name = as.character(pathway_name),
      pathway_description = as.character(pathway_description)
    )
  
  # Check if data is empty
  if (nrow(df) == 0) {
    warning("No pathways found in input data.")
    return(
      significant_modules %>%
        dplyr::mutate(
          matrix_relevance_score = NA_integer_,
          matrix_relevance_reason = "No pathways to analyze.",
          matrix_source = sample_source
        )
    )
  }
  
  message("Analyzing ", nrow(df), " unique pathways for matrix relevance in '", sample_source, "'...")
  
  # Build prompts
  system_prompt <- .build_system_prompt()
  user_payload_json <- .build_user_prompt(df, sample_source)
  
  messages <- list(
    list(role = "system", content = system_prompt),
    list(
      role = "user",
      content = paste0(
        "Sample Source: ", sample_source, "\n",
        "Evaluate the matrix detectability for the following pathways (as JSON). Each reason <= 2 sentences.\n",
        user_payload_json
      )
    )
  )
  
  # Call API
  resp <- .call_openai(
    messages = messages,
    model = model,
    temperature = temperature,
    max_tokens = max_tokens,
    api_key = api_key
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
  
  result <- .parse_and_merge(content, original_df = significant_modules, sample_source = sample_source)
  
  message("Analysis complete. ", sum(!is.na(result$matrix_relevance_score)), " pathways scored.")
  
  return(result)
}
