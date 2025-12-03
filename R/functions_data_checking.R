#' Check the Format and Integrity of a Feature Table
#'
#' Performs a series of checks on a feature table to ensure it meets
#' the expected structure and content required for downstream analysis.
#' The table must contain exactly six columns with specific names and types.
#'
#' @param feature_table A `data.frame` with exactly six columns named:
#'   `variable_id`, `mz`, `rt`, `condition`, `polarity`, and `mean_intensity`.
#'   Each column must conform to specific data types and constraints.
#'
#' @return Returns `TRUE` (invisibly) if all checks pass. Otherwise, the function
#'   stops with an informative error message.
#'
#' @details
#' The function checks:
#' \itemize{
#'   \item That `feature_table` is a `data.frame`
#'   \item That it has exactly six columns
#'   \item That column names and order match the expected names
#'   \item That data types are correct for each column
#'   \item That `variable_id` values are unique
#'   \item That `polarity` only contains "positive" or "negative"
#' }
#'
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @examples
#' # Example feature table
#' feature_table <- data.frame(
#'   variable_id = c("M1", "M2"),
#'   mz = c(100.1, 200.2),
#'   rt = c(300, 400),
#'   condition = c(1, 2),
#'   polarity = c("positive", "negative"),
#'   mean_intensity = c(123456, 654321)
#' )
#'
#' check_feature_table(feature_table)
#'
#' @export

check_feature_table <-
  function(feature_table) {
    # Check class
    if (!is.data.frame(feature_table)) {
      stop("feature_table must be a data.frame.")
    }
    
    # Check number of columns
    if (ncol(feature_table) != 6) {
      stop("feature_table must have exactly 5 columns.")
    }
    
    # Expected column names
    expected_names <- c("variable_id",
                        "mz",
                        "rt",
                        "condition",
                        "polarity",
                        "mean_intensity")
    if (!all(expected_names == colnames(feature_table))) {
      stop(paste0(
        "feature_table columns must be named: ",
        paste(expected_names, collapse = ", ")
      ))
    }
    
    # Check variable_id
    if (!is.character(feature_table$variable_id)) {
      stop("variable_id column must be character.")
    }
    if (anyDuplicated(feature_table$variable_id)) {
      stop("variable_id column contains duplicated values.")
    }
    
    # Check mz
    if (!is.numeric(feature_table$mz)) {
      stop("mz column must be numeric.")
    }
    
    # Check rt
    if (!is.numeric(feature_table$rt)) {
      stop("rt column must be numeric.")
    }
    
    # Check condition
    if (!is.numeric(feature_table$condition)) {
      stop("condition column must be numeric.")
    }
    
    # Check polarity
    if (!is.character(feature_table$polarity)) {
      stop("polarity column must be character.")
    }
    
    # Check mean_intensity
    if (!is.numeric(feature_table$mean_intensity)) {
      stop("polarity column must be numeric.")
    }
    
    allowed_polarities <- c("positive", "negative")
    if (!all(feature_table$polarity %in% allowed_polarities)) {
      stop("polarity column can only contain 'positive' or 'negative'.")
    }
    
    message("✅ feature_table passed all checks.")
    invisible(TRUE)
  }



#' Check the Format of an MS2 Annotation Table
#'
#' Validates that the MS2 annotation table meets required structural and type specifications.
#' The function checks for required columns and ensures proper data types for each.
#'
#' @param check_annotation_table_ms2 A `data.frame` containing MS2 annotation results.
#'   It must include the following columns: `variable_id`, `mz`, `rt`, `HMDB.ID`, and `KEGG.ID`.
#'
#' @return Returns `TRUE` (invisibly) if all checks pass. Otherwise, stops with an informative error message.
#'
#' @details The function performs the following checks:
#' \itemize{
#'   \item `check_annotation_table_ms2` must be a data frame
#'   \item Required columns: `variable_id`, `mz`, `rt`, `HMDB.ID`, `KEGG.ID`
#'   \item `variable_id`: character, unique
#'   \item `mz` and `rt`: numeric
#'   \item `HMDB.ID` and `KEGG.ID`: character (can include `NA`)
#' }
#'
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @examples
#' # Example annotation table
#' annotation_table <- data.frame(
#'   variable_id = c("M1", "M2"),
#'   mz = c(100.1, 200.2),
#'   rt = c(300, 400),
#'   HMDB.ID = c("HMDB00001", NA),
#'   KEGG.ID = c("C00031", NA),
#'   stringsAsFactors = FALSE
#' )
#'
#' check_annotation_table_ms2(annotation_table)
#'
#' @export

check_annotation_table_ms2 <-
  function(check_annotation_table_ms2) {
    # Check that it is a data.frame
    if (!is.data.frame(check_annotation_table_ms2)) {
      stop("check_annotation_table_ms2 must be a data.frame.")
    }
    
    # Required columns
    required_cols <- c("variable_id", "mz", "rt", "HMDB.ID", "KEGG.ID")
    missing_cols <- setdiff(required_cols, colnames(check_annotation_table_ms2))
    if (length(missing_cols) > 0) {
      stop(paste(
        "check_annotation_table_ms2 is missing required columns:",
        paste(missing_cols, collapse = ", ")
      ))
    }
    
    # Check variable_id
    if (!is.character(check_annotation_table_ms2$variable_id)) {
      stop("variable_id column must be character.")
    }
    if (anyDuplicated(check_annotation_table_ms2$variable_id)) {
      stop("variable_id column contains duplicated values.")
    }
    
    # Check mz
    if (!is.numeric(check_annotation_table_ms2$mz)) {
      stop("mz column must be numeric.")
    }
    
    # Check rt
    if (!is.numeric(check_annotation_table_ms2$rt)) {
      stop("rt column must be numeric.")
    }
    
    # Check HMDB.ID
    if (!is.character(check_annotation_table_ms2$HMDB.ID)) {
      stop(
        "HMDB.ID column must be character. It can contain NA but the column itself must be character type."
      )
    }
    
    # Check KEGG.ID
    if (!is.character(check_annotation_table_ms2$KEGG.ID)) {
      stop(
        "KEGG.ID column must be character. It can contain NA but the column itself must be character type."
      )
    }
    
    message("✅ check_annotation_table_ms2 passed all checks.")
    invisible(TRUE)
  }


#' #' @title Check data format for fMSEA
#' #' @description Check data format for fMSEA
#' #' @author Xiaotao Shen
#' #' \email{xiaotao.shen@@outlook.com}
#' #' @param feature_table object list, the items should be mass_dataset class from massdataset
#' #' @return Checking results
#' #' @export
#'
#' check4fmsea <- function(feature_table) {
#'   if (!is.list(feature_table)) {
#'     stop("feature_table should be a list.")
#'   }
#'
#'   check_is_mass_dataset <-
#'     feature_table %>%
#'     lapply(function(x) {
#'       is(x, class2 = "mass_dataset")
#'     }) %>%
#'     unlist()
#'
#'   if (any(!check_is_mass_dataset)) {
#'     stop(paste(which(!check_is_mass_dataset), collapse = ","),
#'          " in feature_table is (are) not mass_dataset.")
#'   }
#'
#'   column_polarity_check <-
#'     feature_table %>%
#'     purrr::map(
#'       .f = function(x) {
#'         variable_info <- slot(x, "variable_info")
#'         if (all(colnames(variable_info) != "column")) {
#'           error1 <- c("no column")
#'           column <- NA
#'         } else{
#'           error1 <- c("no")
#'           column <- unique(variable_info$column)
#'         }
#'
#'         if (all(colnames(variable_info) != "polarity")) {
#'           error2 <- c("no polarity")
#'           polarity <- NA
#'         } else{
#'           error2 <- c("no")
#'           polarity <- unique(variable_info$polarity)
#'         }
#'
#'         list(
#'           error = c(error1, error2),
#'           column = column,
#'           polarity = polarity
#'         )
#'       }
#'     )
#'
#'   error <-
#'     column_polarity_check %>%
#'     purrr::map(function(x) {
#'       x$error
#'     }) %>%
#'     do.call(rbind, .) %>%
#'     as.data.frame()
#'
#'   seq_len(nrow(error)) %>%
#'     lapply(function(i) {
#'       temp_error <-
#'         as.character(error[i, ])
#'
#'       if (any(temp_error != "no")) {
#'         temp_error[temp_error != "no"] %>%
#'           purrr::walk(function(x) {
#'             stop("object ", i , " has ", x)
#'           })
#'       }
#'     })
#'
#'   column <-
#'     column_polarity_check %>%
#'     purrr::map(function(x) {
#'       x$column
#'     })
#'
#'   seq_along(column) %>%
#'     purrr::walk(function(i) {
#'       if (length(column[[i]]) > 1) {
#'         stop("object ",
#'              i,
#'              " column are not unique: ",
#'              paste(column[[i]], collapse = ";"))
#'       }
#'       if (!column[[i]] %in% c("RPLC", "HILIC")) {
#'         stop("object ", i, " column should be RPLC or HILIC: ", column[[i]])
#'       }
#'     })
#'
#'
#'   polarity <-
#'     column_polarity_check %>%
#'     purrr::map(function(x) {
#'       x$polarity
#'     })
#'
#'   seq_along(polarity) %>%
#'     purrr::walk(function(i) {
#'       if (length(polarity[[i]]) > 1) {
#'         stop("object ",
#'              i,
#'              " polarity are not unique: ",
#'              paste(polarity[[i]], collapse = ";"))
#'       }
#'       if (!polarity[[i]] %in% c("positive", "negative")) {
#'         stop("object ",
#'              i,
#'              " polarity should be positive or negative: ",
#'              polarity[[i]])
#'       }
#'     })
#'
#' }
#'
#'
#'
#'
#'
#'
#'
