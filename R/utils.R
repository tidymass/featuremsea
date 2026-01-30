get_extension = function(file) {
  utils::tail(stringr::str_split(string = file, pattern = "\\.")[[1]], 1)
}

#' Read a Table from CSV or Excel Files
#'
#' Reads a data table from a file, automatically detecting the format
#' based on its extension (`csv`, `xlsx`, or `xls`). Supports additional
#' arguments passed to the underlying reading functions.
#'
#' @param file Character. The file path to read.
#' @param ... Additional parameters passed to `readr::read_csv()` or
#'   `readxl::read_xlsx()` / `readxl::read_xls()`.
#'
#' @return A data frame containing the table contents.
#'
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @examples
#' \dontrun{
#' # Read CSV file
#' df_csv <- readTable("data/example.csv")
#'
#' # Read Excel file (.xlsx)
#' df_xlsx <- readTable("data/example.xlsx")
#' }
readTable = function(file, ...) {
  extension <- get_extension(file = file)
  if (extension == "csv") {
    return(readr::read_csv(file = file, show_col_types = FALSE, ...))
  }
  
  if (extension == 'xlsx') {
    return(readxl::read_xlsx(path = file, ...))
  }
  
  if (extension == "xls") {
    return(readxl::read_xls(path = file, ...))
  }
  
  if (extension != "csv" &
      extension != "xlsx" &
      extension != "xls") {
    message(crayon::red("file is not csv, xlsx or xls."))
  }
}


msg <- function(..., startup = FALSE) {
  if (startup) {
    if (!isTRUE(getOption("fmsea.quiet"))) {
      packageStartupMessage(text_col(...))
    }
  } else {
    message(text_col(...))
  }
}

text_col <- function(x) {
  # If RStudio not available, messages already printed in black
  if (!rstudioapi::isAvailable()) {
    return(x)
  }
  
  if (!rstudioapi::hasFun("getThemeInfo")) {
    return(x)
  }
  
  theme <- rstudioapi::getThemeInfo()
  
  if (isTRUE(theme$dark))
    crayon::white(x)
  else
    crayon::black(x)
  
}

#' List all packages in the fmsea
#'
#' @param include_self Include fmsea in the list?
#' @export
#' @return fmsea_packages
#' @examples
#' fmsea_packages()
fmsea_packages <- function(include_self = TRUE) {
  raw <- utils::packageDescription("featuremsea")$Imports
  imports <- strsplit(raw, ",")[[1]]
  parsed <- gsub("^\\s+|\\s+$", "", imports)
  names <-
    vapply(strsplit(parsed, "\\s+"), "[[", 1, FUN.VALUE = character(1))
  
  if (include_self) {
    names <- c(names, "featuremsea")
  }
  
  names
}

invert <- function(x) {
  if (length(x) == 0)
    return()
  stacked <- utils::stack(x)
  tapply(as.character(stacked$ind), stacked$values, list)
}


style_grey <- function(level, ...) {
  crayon::style(paste0(...),
                crayon::make_style(grDevices::grey(level), grey = TRUE))
}
