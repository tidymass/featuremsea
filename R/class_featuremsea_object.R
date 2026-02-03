#' @title featuremsea_object Class
#' @description An S4 class to store the results of feature-based MSEA analysis.
#' @slot ranking_table A data frame of the original feature rankings.
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
    ranking_table = "data.frame",            # 新增
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
#' @importFrom utils packageVersion head
#' @importFrom cli col_green style_bold
#' @export
setMethod(
  f = "show",
  signature = "featuremsea_object",
  definition = function(object) {
    
    # 获取维度信息
    n_ranking <- nrow(object@ranking_table)
    n_sig     <- nrow(object@significant_modules)
    n_fm      <- nrow(object@feature_metabolite_count)
    n_anno_r  <- nrow(object@annotation_table_weighting)
    n_anno_c  <- ncol(object@annotation_table_weighting)
    n_list    <- length(object@res_list)
    
    pkg_ver <- tryCatch(as.character(utils::packageVersion("featuremsea")), 
                        error = function(e) "0.0.1")
    
    cat("--------------------\n")
    cat(cli::col_green(paste0("featuremsea version: ", pkg_ver)), "\n")
    cat("--------------------\n")
    
    # 1. ranking_table (新增展示)
    cat(cli::style_bold(cli::col_green("1. ranking_table:")), 
        sprintf("[ %d features with weights ]\n", n_ranking))
    
    # 2. feature_metabolite_count
    cat(cli::style_bold(cli::col_green("2. feature_metabolite_count:")), 
        sprintf("[ %d pairs retained ]\n", n_fm))
    
    # 3. annotation_table_weighting
    cat(cli::style_bold(cli::col_green("3. annotation_table_weighting:")), 
        sprintf("[ %d x %d data.frame ]\n", n_anno_r, n_anno_c))
    
    # 4. significant_modules
    cat(cli::style_bold(cli::col_green("4. significant_modules:")), 
        sprintf("[ %d significant modules found ]\n", n_sig))
    
    if (n_sig > 0) {
      # 注意：此处列名已由 MFM_name 更改为 pathway_name
      top_names <- head(object@significant_modules$pathway_name, 3)
      cat(paste("    Top modules:", paste(top_names, collapse = ", ")))
      if (n_sig > 3) cat(" ...")
      cat("\n")
    }
    
    # 5. res_list
    cat(cli::style_bold(cli::col_green("5. res_list:")), 
        sprintf("[ List with %d elements ]\n", n_list))
    
    # 6. Algorithm status
    cat(cli::style_bold(cli::col_green("6. algorithm_status:")), "\n")
    status_text <- if (object@converged) "YES" else "NO"
    cat(sprintf("    Converged: %s\n", cli::col_green(status_text)))
    cat(sprintf("    Iterations used: %d\n", object@iterations_used))
    
    cat("--------------------\n")
    
    if(length(object@process_info) > 0){
      cat(cli::style_bold(cli::col_green("Processing information")), "\n")
      cat(sprintf("    Created at: %s\n", object@process_info$creation_date))
    }
  }
)