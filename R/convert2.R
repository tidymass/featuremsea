#' Convert metpath KEGG object to fmsea database format
#'
#' @param kegg_hsa_pathway A KEGG pathway object from the metpath package
#' @return A data frame formatted for featuremsea/fmsea
#' @export
convert_kegg2fmsea <- function(kegg_hsa_pathway) {
  
  # Ensure required packages are available
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("stringr is required.")
  
  # 1. Extract Description
  # Note: Accessing @describtion slot (preserving original typo if present in S4 class)
  description_vec <- vapply(
    kegg_hsa_pathway@describtion,
    function(x) {
      if (length(x) == 0) {
        NA_character_
      } else {
        as.character(x[[1]])
      }
    },
    FUN.VALUE = character(1)
  )
  
  # 2. Extract Pathway Class
  pathway_class_vec <- vapply(
    kegg_hsa_pathway@pathway_class,
    function(x) {
      if (length(x) == 0) {
        NA_character_
      } else {
        as.character(x[[1]])
      }
    },
    FUN.VALUE = character(1)
  )
  
  # 3. Extract Compound IDs (collapsed with "{}")
  compound_list_vec <- vapply(
    kegg_hsa_pathway@compound_list,
    function(x) {
      if (length(x$KEGG.ID) == 0) {
        NA_character_
      } else {
        paste(x$KEGG.ID, collapse = "{}")
      }
    },
    FUN.VALUE = character(1)
  )
  
  # 4. Extract Compound Names (collapsed with "{}")
  compound_name_vec <- vapply(
    kegg_hsa_pathway@compound_list,
    function(x) {
      if (length(x$Compound.name) == 0) {
        NA_character_
      } else {
        paste(x$Compound.name, collapse = "{}")
      }
    },
    FUN.VALUE = character(1)
  )
  
  # 5. Construct initial data frame
  kegg_database <- data.frame(
    pathway_id    = kegg_hsa_pathway@pathway_id,
    pathway_name  = kegg_hsa_pathway@pathway_name,
    description   = description_vec,
    compound_name = compound_name_vec,
    compound_list = compound_list_vec,
    pathway_class = pathway_class_vec,
    database      = "KEGG",
    stringsAsFactors = FALSE
  )
  
  # 6. Final formatting and renaming
  kegg_database <- kegg_database %>%
    dplyr::mutate(
      pathway_class_all = stringr::str_replace_all(pathway_class, "; ", "{}")
    ) %>%
    dplyr::select(-pathway_class) %>%
    dplyr::rename(
      KEGG_id          = pathway_id,
      KEGG_name        = pathway_name,
      KEGG_description = description
    )
  
  return(kegg_database)
}