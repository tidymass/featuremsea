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
      pathway_description = description,
      KEGG_ID          = compound_list
    )
  
  return(kegg_database)
}



#' Convert metpath HMDB object to fmsea database format
#'
#' @param hmdb_pathway An HMDB pathway object from the metpath package
#' @return A data frame formatted for featuremsea/fmsea
#' @export
convert_hmdb2fmsea <- function(hmdb_pathway) {
  
  # Ensure required packages are available
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("stringr is required.")
  
  # 1. Extract Description
  # Note: Accessing @describtion slot (preserving original typo as seen in metpath S4 classes)
  description_vec <- vapply(
    hmdb_pathway@describtion,
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
    hmdb_pathway@pathway_class,
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
  # Note: Changing target column to HMDB.ID
  compound_list_vec <- vapply(
    hmdb_pathway@compound_list,
    function(x) {
      if (length(x$HMDB.ID) == 0) {
        NA_character_
      } else {
        paste(x$HMDB.ID, collapse = "{}")
      }
    },
    FUN.VALUE = character(1)
  )
  
  # 4. Extract Compound Names (collapsed with "{}")
  compound_name_vec <- vapply(
    hmdb_pathway@compound_list,
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
  hmdb_database <- data.frame(
    pathway_id    = hmdb_pathway@pathway_id,
    pathway_name  = hmdb_pathway@pathway_name,
    description   = description_vec,
    compound_name = compound_name_vec,
    compound_list = compound_list_vec,
    pathway_class = pathway_class_vec,
    database      = "HMDB",
    stringsAsFactors = FALSE
  )
  
  # 6. Final formatting and renaming
  hmdb_database <- hmdb_database %>%
    dplyr::mutate(
      pathway_class_all = stringr::str_replace_all(pathway_class, ";", "{}")
    ) %>%
    dplyr::select(-pathway_class) %>%
    dplyr::rename(
      pathway_description = description,
      HMDB_ID             = compound_list
    )
  
  return(hmdb_database)
}

#' Convert iMetPD Pathway Object to featuremsea Format
#'
#' This function extracts pathway information from an iMetPD S4 object and formats it
#' into a data frame suitable for the featuremsea package. It processes HMDB and KEGG IDs,
#' removing NAs and duplicates, and collapses them using "{}".
#'
#' @param imetpd_pathway An S4 object containing iMetPD pathway data. It is expected to have
#' slots for `describtion`, `pathway_class`, `compound_list` (containing HMDB.ID and KEGG.ID),
#' `pathway_id`, and `pathway_name`.
#'
#' @return A data frame containing the following columns:
#' \itemize{
#'   \item pathway_id
#'   \item pathway_name
#'   \item pathway_description
#'   \item HMDB_ID (Collapsed unique IDs)
#'   \item KEGG_ID (Collapsed unique IDs)
#'   \item database (Set to "IMETPD")
#'   \item pathway_class_all
#' }
#' @importFrom dplyr mutate select rename
#' @importFrom stringr str_replace_all
#' @export
convert_imetpd2fmsea <- function(imetpd_pathway) {
  
  # Ensure required packages are available
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("stringr is required.")
  
  # 1. Extract Description
  description_vec <- vapply(
    imetpd_pathway@describtion,
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
    imetpd_pathway@pathway_class,
    function(x) {
      if (length(x) == 0) {
        NA_character_
      } else {
        as.character(x[[1]])
      }
    },
    FUN.VALUE = character(1)
  )
  
  # 3. Extract HMDB IDs (Remove NA, Deduplicate, Collapse with "{}")
  hmdb_id_vec <- vapply(
    imetpd_pathway@compound_list,
    function(x) {
      # Extract HMDB.ID column
      ids <- as.character(x$HMDB.ID)
      # Remove NAs and empty strings
      ids <- ids[!is.na(ids) & ids != ""]
      # Deduplicate
      ids <- unique(ids)
      
      if (length(ids) == 0) {
        NA_character_
      } else {
        paste(ids, collapse = "{}")
      }
    },
    FUN.VALUE = character(1)
  )
  
  # 4. Extract KEGG IDs (Remove NA, Deduplicate, Collapse with "{}")
  kegg_id_vec <- vapply(
    imetpd_pathway@compound_list,
    function(x) {
      # Extract KEGG.ID column
      ids <- as.character(x$KEGG.ID)
      # Remove NAs and empty strings
      ids <- ids[!is.na(ids) & ids != ""]
      # Deduplicate
      ids <- unique(ids)
      
      if (length(ids) == 0) {
        NA_character_
      } else {
        paste(ids, collapse = "{}")
      }
    },
    FUN.VALUE = character(1)
  )
  
  # 5. Construct initial data frame
  imetpd_database <- data.frame(
    pathway_id    = imetpd_pathway@pathway_id,
    pathway_name  = imetpd_pathway@pathway_name,
    pathway_description = description_vec,
    pathway_class = pathway_class_vec,
    HMDB_ID       = hmdb_id_vec,
    KEGG_ID       = kegg_id_vec,
    database      = "IMETPD",
    stringsAsFactors = FALSE
  )
  
  # 6. Final formatting
  imetpd_database <- imetpd_database %>%
    dplyr::mutate(
      pathway_class_all = stringr::str_replace_all(pathway_class, ";", "{}")
    ) %>%
    dplyr::select(-pathway_class) 
  
  return(imetpd_database)
}


#' Convert Reactome Pathway Object to featuremsea Format
#'
#' This function extracts pathway information from a Reactome S4 object and formats it
#' into a data frame suitable for the featuremsea package. 
#' It assumes the compound list has already been mapped to HMDB and KEGG IDs.
#'
#' @param reactome_pathway An S4 object containing Reactome pathway data. It is expected to have
#' slots for `describtion`, `pathway_class`, `compound_list` (containing HMDB.ID and KEGG.ID),
#' `pathway_id`, and `pathway_name`.
#'
#' @return A data frame containing the following columns:
#' \itemize{
#'   \item pathway_id
#'   \item pathway_name
#'   \item pathway_description
#'   \item HMDB_ID (Collapsed unique IDs)
#'   \item KEGG_ID (Collapsed unique IDs)
#'   \item database (Set to "Reactome")
#'   \item pathway_class_all
#' }
#' @importFrom dplyr mutate select rename
#' @importFrom stringr str_replace_all
#' @export
convert_reactome2fmsea <- function(reactome_pathway) {
  
  # Ensure required packages are available
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("stringr is required.")
  
  # 1. Extract Description
  description_vec <- vapply(
    reactome_pathway@describtion,
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
    reactome_pathway@pathway_class,
    function(x) {
      if (length(x) == 0) {
        NA_character_
      } else {
        as.character(x[[1]])
      }
    },
    FUN.VALUE = character(1)
  )
  
  # 3. Extract HMDB IDs (Remove NA, Deduplicate, Collapse with "{}")
  hmdb_id_vec <- vapply(
    reactome_pathway@compound_list,
    function(x) {
      # Check if column exists, return NA if not
      if (!"HMDB.ID" %in% colnames(x)) return(NA_character_)
      
      # Extract HMDB.ID column
      ids <- as.character(x$HMDB.ID)
      # Remove NAs and empty strings
      ids <- ids[!is.na(ids) & ids != ""]
      # Deduplicate
      ids <- unique(ids)
      
      if (length(ids) == 0) {
        NA_character_
      } else {
        paste(ids, collapse = "{}")
      }
    },
    FUN.VALUE = character(1)
  )
  
  # 4. Extract KEGG IDs (Remove NA, Deduplicate, Collapse with "{}")
  kegg_id_vec <- vapply(
    reactome_pathway@compound_list,
    function(x) {
      # Check if column exists, return NA if not
      if (!"KEGG.ID" %in% colnames(x)) return(NA_character_)
      
      # Extract KEGG.ID column
      ids <- as.character(x$KEGG.ID)
      # Remove NAs and empty strings
      ids <- ids[!is.na(ids) & ids != ""]
      # Deduplicate
      ids <- unique(ids)
      
      if (length(ids) == 0) {
        NA_character_
      } else {
        paste(ids, collapse = "{}")
      }
    },
    FUN.VALUE = character(1)
  )
  
  # 5. Construct initial data frame
  reactome_database <- data.frame(
    pathway_id    = reactome_pathway@pathway_id,
    pathway_name  = reactome_pathway@pathway_name,
    pathway_description = description_vec,
    pathway_class = pathway_class_vec,
    HMDB_ID       = hmdb_id_vec,
    KEGG_ID       = kegg_id_vec,
    database      = "Reactome",
    stringsAsFactors = FALSE
  )
  
  # 6. Final formatting
  reactome_database <- reactome_database %>%
    dplyr::mutate(
      pathway_class_all = stringr::str_replace_all(pathway_class, ";", "{}")
    ) %>%
    dplyr::select(-pathway_class) 
  
  return(reactome_database)
}


#' Convert WikiPathways Object to featuremsea Format
#'
#' This function extracts pathway information from a WikiPathways S4 object and formats it
#' into a data frame suitable for the featuremsea package.
#' It assumes the compound list has already been mapped to HMDB and KEGG IDs.
#'
#' @param wiki_pathway An S4 object containing WikiPathways data. It is expected to have
#' slots for `describtion`, `pathway_class`, `compound_list` (containing HMDB.ID and KEGG.ID),
#' `pathway_id`, and `pathway_name`.
#'
#' @return A data frame containing the following columns:
#' \itemize{
#'   \item pathway_id
#'   \item pathway_name
#'   \item pathway_description
#'   \item HMDB_ID (Collapsed unique IDs)
#'   \item KEGG_ID (Collapsed unique IDs)
#'   \item database (Set to "WikiPathways")
#'   \item pathway_class_all
#' }
#' @importFrom dplyr mutate select rename
#' @importFrom stringr str_replace_all
#' @export
convert_wikipathway2fmsea <- function(wiki_pathway) {
  
  # Ensure required packages are available
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr is required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("stringr is required.")
  
  # 1. Extract Description
  description_vec <- vapply(
    wiki_pathway@describtion,
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
    wiki_pathway@pathway_class,
    function(x) {
      if (length(x) == 0) {
        NA_character_
      } else {
        as.character(x[[1]])
      }
    },
    FUN.VALUE = character(1)
  )
  
  # 3. Extract HMDB IDs (Remove NA, Deduplicate, Collapse with "{}")
  hmdb_id_vec <- vapply(
    wiki_pathway@compound_list,
    function(x) {
      # Check if column exists, return NA if not
      if (!"HMDB.ID" %in% colnames(x)) return(NA_character_)
      
      # Extract HMDB.ID column
      ids <- as.character(x$HMDB.ID)
      # Remove NAs and empty strings
      ids <- ids[!is.na(ids) & ids != ""]
      # Deduplicate
      ids <- unique(ids)
      
      if (length(ids) == 0) {
        NA_character_
      } else {
        paste(ids, collapse = "{}")
      }
    },
    FUN.VALUE = character(1)
  )
  
  # 4. Extract KEGG IDs (Remove NA, Deduplicate, Collapse with "{}")
  kegg_id_vec <- vapply(
    wiki_pathway@compound_list,
    function(x) {
      # Check if column exists, return NA if not
      if (!"KEGG.ID" %in% colnames(x)) return(NA_character_)
      
      # Extract KEGG.ID column
      ids <- as.character(x$KEGG.ID)
      # Remove NAs and empty strings
      ids <- ids[!is.na(ids) & ids != ""]
      # Deduplicate
      ids <- unique(ids)
      
      if (length(ids) == 0) {
        NA_character_
      } else {
        paste(ids, collapse = "{}")
      }
    },
    FUN.VALUE = character(1)
  )
  
  # 5. Construct initial data frame
  wiki_database <- data.frame(
    pathway_id    = wiki_pathway@pathway_id,
    pathway_name  = wiki_pathway@pathway_name,
    pathway_description = description_vec,
    pathway_class = pathway_class_vec,
    HMDB_ID       = hmdb_id_vec,
    KEGG_ID       = kegg_id_vec,
    database      = "Wikipathways",
    stringsAsFactors = FALSE
  )
  
  # 6. Final formatting
  wiki_database <- wiki_database %>%
    dplyr::mutate(
      pathway_class_all = stringr::str_replace_all(pathway_class, ";", "{}")
    ) %>%
    dplyr::select(-pathway_class) 
  
  return(wiki_database)
}