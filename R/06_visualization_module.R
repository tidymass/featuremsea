#' Plot fMSEA GSEA Enrichment Results
#'
#' @description
#' This function visualizes the Functional Metabolite Set Enrichment Analysis (fMSEA) 
#' results using the classic GSEA three-panel plot. It utilizes the pre-calculated 
#' `steps` data to ensure the visualization exactly matches the analytical results.
#'
#' @param fmsea_obj An fMSEA result object. Must contain a `res_list` slot 
#'   (with detailed step data) and a `significant_modules` slot (with summary statistics).
#' @param pathway_id Character. The specific pathway ID to be plotted (e.g., "hsa00140").
#' @param title Character. The main title of the plot. If \code{NULL} (default), 
#'   the pathway name is automatically extracted from the results.
#'
#' @return A \code{patchwork} combined \code{ggplot2} object consisting of:
#' \itemize{
#'   \item \strong{Top:} Running Enrichment Score (ES) curve.
#'   \item \strong{Middle:} Barcode plot showing the position of pathway hits.
#'   \item \strong{Bottom:} Distribution of the ranking metric/step values.
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_linerange geom_area labs theme_classic theme_void theme element_blank element_text element_rect
#' @importFrom patchwork plot_layout
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'results' is your fMSEA output object
#' p <- plot_fmsea_plot(results, pathway_id = "hsa00140")
#' print(p)
#' }
plot_fmsea_plot <- function(fmsea_obj, 
                            pathway_id, 
                            title = NULL) {
  
  # --- 1. Data Extraction ---
  
  # Verify if the pathway_id exists in the result list
  if (!pathway_id %in% names(fmsea_obj@res_list)) {
    stop(paste("Pathway ID", pathway_id, "not found in results@res_list."))
  }
  
  # Extract step-by-step calculation data for the specific pathway
  step_data <- fmsea_obj@res_list[[pathway_id]]$steps
  
  # --- 2. Statistical Metadata (NES and FDR) ---
  
  # Extract summary statistics for annotation
  stats <- fmsea_obj@significant_modules
  
  # Handle potential column name variations (MFM_id vs pathway_id)
  id_col <- if("MFM_id" %in% colnames(stats)) "MFM_id" else "pathway_id"
  pathway_stats <- stats[stats[[id_col]] == pathway_id, ]
  
  if (nrow(pathway_stats) > 0) {
    # Format NES and FDR values
    nes_val <- round(pathway_stats$NES, 3)
    fdr_val <- format(pathway_stats$FDR, scientific = TRUE, digits = 3)
    sub_title <- paste0("NES: ", nes_val, " | FDR: ", fdr_val)
    
    # Extract pathway name for the title
    pathway_name <- if("MFM_name" %in% colnames(pathway_stats)) {
      pathway_stats$MFM_name 
    } else {
      pathway_id
    }
  } else {
    sub_title <- "No significant stats found"
    pathway_name <- pathway_id
  }
  
  # Define final plot title
  plot_title <- if(is.null(title)) pathway_name else title
  
  # --- 3. Visualization Components ---
  
  # P1: Running Enrichment Score Curve
  p1 <- ggplot2::ggplot(step_data, ggplot2::aes(x = order_pos, y = ES_cum)) +
    ggplot2::geom_line(color = "#00AFBB", linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(title = plot_title, subtitle = sub_title, y = "Enrichment Score") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                   axis.text.x = ggplot2::element_blank(), 
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold", size = 12))
  
  # P2: Hit Barcode Plot
  # Filter only rows where a metabolite hit occurred ('pair_slot')
  hits_data <- step_data[step_data$type == "pair_slot", ]
  
  p2 <- ggplot2::ggplot(hits_data, ggplot2::aes(x = order_pos, y = 1)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = 0, ymax = 1), color = "black", alpha = 0.6) +
    ggplot2::labs(y = "") + 
    ggplot2::theme_void() +
    ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5))
  
  # P3: Ranking Metric / Step Contribution
  p3 <- ggplot2::ggplot(step_data, ggplot2::aes(x = order_pos, y = step)) +
    ggplot2::geom_area(fill = "grey80", alpha = 0.5) + 
    ggplot2::geom_line(color = "grey40", linewidth = 0.5) +
    ggplot2::labs(x = "Rank in Ordered Dataset", y = "Step Value") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey") + 
    ggplot2::theme_classic()
  
  # --- 4. Assembly ---
  
  # Combine plots vertically using patchwork
  combined_plot <- (p1 / p2 / p3) + 
    patchwork::plot_layout(heights = c(3, 0.5, 1.5))
  
  return(combined_plot)
}