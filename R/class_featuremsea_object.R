#' @title featuremsea_object Class
#' @description An S4 class to store the results of feature-based MSEA analysis.
#' @slot feature_metabolite_count A data frame of feature-metabolite counts.
#' @slot annotation_table_weighting A data frame of weighted annotations.
#' @slot significant_modules A data frame of significant modules.
#' @slot res_list A list containing raw results.
#' @slot converged Logical. Whether the algorithm converged.
#' @slot iterations_used Integer. Number of iterations used.
#' @slot process_info A list containing processing timestamp/info.
#' @export
setClass(
  Class = "featuremsea_object",
  slots = c(
    feature_metabolite_count = "data.frame",
    annotation_table_weighting = "data.frame",
    significant_modules = "data.frame",
    res_list = "list",
    converged = "logical",
    iterations_used = "integer",
    process_info = "list"
  )
)

#' @title Show method for featuremsea_object
#' @param object An object of class featuremsea_object.
#' @importFrom methods show
#' @importFrom utils packageVersion
#' @importFrom cli col_green style_bold
#' @export
setMethod(
  f = "show",
  signature = "featuremsea_object",
  definition = function(object) {
    
    # Get dimension information
    n_sig <- nrow(object@significant_modules)
    n_fm <- nrow(object@feature_metabolite_count)
    n_anno_r <- nrow(object@annotation_table_weighting)
    n_anno_c <- ncol(object@annotation_table_weighting)
    n_list <- length(object@res_list)
    
    # Get package version (default to 0.0.1 if package is not installed)
    pkg_ver <- tryCatch(as.character(utils::packageVersion("featuremsea")), 
                        error = function(e) "0.0.1")
    
    cat("--------------------\n")
    # Display version number (Green)
    cat(cli::col_green(paste0("featuremsea version: ", pkg_ver)), "\n")
    cat("--------------------\n")
    
    # 1. feature_metabolite_count (Green)
    cat(cli::style_bold(cli::col_green("1. feature_metabolite_count:")), 
        sprintf("[ %d pairs retained ]\n", n_fm))
    
    # 2. annotation_table_weighting (Green)
    cat(cli::style_bold(cli::col_green("2. annotation_table_weighting:")), 
        sprintf("[ %d x %d data.frame ]\n", n_anno_r, n_anno_c))
    
    # 3. significant_modules (Green)
    cat(cli::style_bold(cli::col_green("3. significant_modules:")), 
        sprintf("[ %d significant modules found ]\n", n_sig))
    
    if (n_sig > 0) {
      top_names <- head(object@significant_modules$MFM_name, 3)
      # Indent module names slightly for cleaner look
      cat(paste("   Top modules:", paste(top_names, collapse = ", ")))
      if (n_sig > 3) cat(" ...")
      cat("\n")
    }
    
    # 4. res_list (Green)
    cat(cli::style_bold(cli::col_green("4. res_list:")), 
        sprintf("[ List with %d elements ]\n", n_list))
    
    # 5. Algorithm status (Green)
    cat(cli::style_bold(cli::col_green("5. algorithm_status:")), "\n")
    
    status_text <- if (object@converged) "YES" else "NO"
    
    # Status text is also Green as requested
    cat(sprintf("   Converged: %s\n", cli::col_green(status_text)))
    cat(sprintf("   Iterations used: %d\n", object@iterations_used))
    
    cat("--------------------\n")
    
    # Processing information (Green)
    if(length(object@process_info) > 0){
      cat(cli::style_bold(cli::col_green("Processing information")), "\n")
      cat(sprintf("   Created at: %s\n", object@process_info$creation_date))
    }
  }
)