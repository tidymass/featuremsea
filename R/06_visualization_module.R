#' Plot fMSEA GSEA Enrichment Results
#'
#' @description
#' Visualizes fMSEA results.
#' Updated visuals: P4, P5, P6 now use stick/segment plots for better clarity of sparse weights.
#' Fixed: X-axis alignment across all panels.
#'
#' @param fmsea_obj An fMSEA result object.
#' @param pathway_id Character. The specific pathway ID to be plotted.
#' @param title Character. The main title of the plot.
#'
#' @return A \code{patchwork} combined \code{ggplot2} object.
#' @export
plot_fmsea_plot <- function(fmsea_obj, 
                            pathway_id, 
                            title = NULL) {
  
  # --- 0. Configuration: Colors ---
  col_pos  <- "#E64B35"  # Red
  col_neg  <- "#3498DB"  # Blue
  col_line <- "#00A087"  # Teal
  
  # --- 1. Data Extraction ---
  if (!pathway_id %in% names(fmsea_obj@res_list)) {
    stop(paste("Pathway ID", pathway_id, "not found."))
  }
  
  # 1.1 Get Steps Data (X-axis source)
  step_data <- fmsea_obj@res_list[[pathway_id]]$steps
  
  # Preprocess P5/P6 columns: Replace NA with 0 immediately
  if("a_ij" %in% colnames(step_data)) step_data$a_ij[is.na(step_data$a_ij)] <- 0
  if("pair_weight" %in% colnames(step_data)) step_data$pair_weight[is.na(step_data$pair_weight)] <- 0
  
  # --- 2. Link with Ranking Table (for Direction) ---
  if (nrow(fmsea_obj@ranking_table) == 0) stop("ranking_table is empty.")
  
  # Reconstruct global order based on absolute weight
  rank_ref <- fmsea_obj@ranking_table %>%
    dplyr::arrange(dplyr::desc(abs(ranking_weight))) %>%
    dplyr::mutate(
      order_pos = dplyr::row_number(),
      real_ranking_weight = ranking_weight,
      norm_weight = abs(ranking_weight) / max(abs(ranking_weight))
    ) %>%
    dplyr::select(order_pos, real_ranking_weight, norm_weight)
  
  # Merge to get signed weights for colors
  # distinct() prevents potential row duplication if ranking table has issues
  plot_data <- step_data %>%
    dplyr::left_join(rank_ref, by = "order_pos")
  
  # --- 3. Sub-plot Data Prep ---
  
  # P2: Hits
  hits_data <- plot_data %>%
    dplyr::filter(!is.na(variable_id_hit) & variable_id_hit != "") %>%
    dplyr::mutate(direction = ifelse(real_ranking_weight > 0, "Positive", "Negative"))
  
  # P4: Global Weight Background
  # Filter rank_ref to match the plotting range
  min_pos <- min(plot_data$order_pos)
  max_pos <- max(plot_data$order_pos)
  weight_dist_data <- rank_ref %>% 
    dplyr::filter(order_pos >= min_pos & order_pos <= max_pos)
  
  # --- 4. Stats & Table ---
  stats_df <- fmsea_obj@significant_modules
  id_col <- if("MFM_id" %in% colnames(stats_df)) "MFM_id" else "pathway_id"
  p_stats <- stats_df[stats_df[[id_col]] == pathway_id, ]
  
  # Helper to safe extract
  safe_val <- function(x, digits=3, sci=FALSE) {
    if(length(x)==0 || is.na(x)) return("-")
    if(is.numeric(x)) {
      if(sci) format(x, scientific=TRUE, digits=digits) else round(x, digits)
    } else as.character(x)
  }
  
  # Metadata
  p_name <- if(nrow(p_stats)>0) (if("pathway_name" %in% names(p_stats)) p_stats$pathway_name else p_stats$MFM_name) else pathway_id
  le_pos <- plot_data$order_pos[which.max(abs(plot_data$ES_cum))]
  
  # Table Data
  table_df <- data.frame(
    Metric = c("ID", "Name", "LE Pos", "P-val", "FDR", "ES", "NES"),
    Value = c(
      pathway_id,
      ifelse(nchar(p_name) > 25, paste0(substr(p_name, 1, 22), "..."), p_name),
      safe_val(le_pos, 0),
      safe_val(p_stats$p_value, sci=TRUE),
      safe_val(p_stats$FDR, sci=TRUE),
      safe_val(p_stats$ES, 4),
      safe_val(p_stats$NES, 4)
    ), stringsAsFactors = FALSE
  )
  
  # Table Grob
  tbl_grob <- gridExtra::tableGrob(table_df, rows=NULL, theme=gridExtra::ttheme_minimal(
    core=list(fg_params=list(hjust=0, x=0.1, fontsize=5), bg_params=list(fill=c("white","grey95"), col="grey60", lwd=0.2)),
    colhead=list(fg_params=list(fontsize=5, fontface="bold"), bg_params=list(fill="grey85", col="grey40"))
  ))
  tbl_grob <- gtable::gtable_add_grob(tbl_grob, grid::rectGrob(gp=grid::gpar(fill=NA, col="grey30")), t=1, l=1, b=nrow(tbl_grob), r=ncol(tbl_grob), z=0)
  
  # --- 5. Visualization Panels ---
  # FIXED: Ensure x-axis alignment across all panels
  x_lim <- c(min_pos, max_pos)
  
  # Unified x-axis scale for all panels
  common_x_scale <- ggplot2::scale_x_continuous(
    limits = x_lim,
    expand = c(0, 0)
  )
  
  # Common theme
  theme_strip <- ggplot2::theme_classic() + 
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(), 
      axis.text.x = ggplot2::element_blank(), 
      axis.ticks.x = ggplot2::element_blank(), 
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(1, 5, 1, 1),
      axis.title.y = ggplot2::element_text(size = 9, angle = 90, vjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 7)
    )
  
  # P1: ES Curve
  y_rng <- range(plot_data$ES_cum)
  tbl_y <- if(abs(max(y_rng)) >= abs(min(y_rng))) c(max(y_rng)*0.5, max(y_rng)*0.95) else c(min(y_rng)*0.1, min(y_rng)*0.5)
  
  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x=order_pos, y=ES_cum)) +
    ggplot2::geom_vline(xintercept=le_pos, linetype="dashed", color=col_pos, alpha=0.8) +
    ggplot2::geom_line(color=col_line, linewidth=1) +
    ggplot2::geom_hline(yintercept=0, color="grey60", linewidth=0.3) +
    ggplot2::annotation_custom(tbl_grob, xmin=x_lim[2]*0.6, xmax=x_lim[2]*1.02, ymin=tbl_y[1], ymax=tbl_y[2]) +
    ggplot2::labs(title=if(is.null(title)) p_name else title, y="ES") +
    common_x_scale +
    theme_strip + 
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold", size=12), axis.title.y=ggplot2::element_text(size=11))
  
  # P2: Hits
  p2 <- ggplot2::ggplot(hits_data, ggplot2::aes(x=order_pos, y=1)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin=0, ymax=1, color=direction)) +
    ggplot2::scale_color_manual(values=c("Positive"=col_pos, "Negative"=col_neg), guide="none") +
    ggplot2::labs(y="Hit") +
    common_x_scale +
    theme_strip +
    ggplot2::theme(axis.line.x=ggplot2::element_blank(), axis.line.y=ggplot2::element_blank(),
                   axis.ticks.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(),
                   panel.border=ggplot2::element_rect(color="black", fill=NA))
  
  # P3: Step (Running metric) - Keep as area/line as it is continuous
  p3 <- ggplot2::ggplot(plot_data, ggplot2::aes(x=order_pos, y=step)) +
    ggplot2::geom_area(fill="grey85", alpha=0.5) +
    ggplot2::geom_line(color="grey40", linewidth=0.4) +
    ggplot2::labs(y="Step") +
    common_x_scale +
    theme_strip
  
  # --- UPDATED: P4, P5, P6 use line plots with original values (no normalization) ---
  
  # P4: Ranking Weight (Abs) - Line plot with original absolute weight
  p4 <- ggplot2::ggplot(weight_dist_data, ggplot2::aes(x=order_pos, y=abs(real_ranking_weight))) +
    ggplot2::geom_line(color="#C0392B", linewidth=0.6) +
    ggplot2::labs(y="Weight") +
    common_x_scale +
    ggplot2::scale_y_continuous(expand=c(0, 0)) +
    theme_strip +
    ggplot2::theme(panel.border=ggplot2::element_rect(color="black", fill=NA))
  
  # P5: a_ij (Adjacency) - Line plot with original a_ij values
  p5 <- ggplot2::ggplot(plot_data, ggplot2::aes(x=order_pos, y=a_ij)) +
    ggplot2::geom_line(color="#6C3483", linewidth=0.6) +
    ggplot2::labs(y="a_ij") +
    common_x_scale +
    ggplot2::scale_y_continuous(expand=c(0, 0)) +
    theme_strip +
    ggplot2::theme(panel.border=ggplot2::element_rect(color="black", fill=NA))
  
  # P6: Pair Weight - Line plot with original pair_weight values
  p6 <- ggplot2::ggplot(plot_data, ggplot2::aes(x=order_pos, y=pair_weight)) +
    ggplot2::geom_line(color="#BA4A00", linewidth=0.6) +
    ggplot2::labs(x="Ordered Features", y="Pair W.") +
    common_x_scale +
    ggplot2::scale_y_continuous(expand=c(0, 0)) +
    theme_strip +
    ggplot2::theme(
      plot.margin = ggplot2::margin(1, 5, 5, 1),
      axis.title.x = ggplot2::element_text(size=10),
      axis.text.x = ggplot2::element_text(size=8),
      axis.ticks.x = ggplot2::element_line(),
      panel.border = ggplot2::element_rect(color="black", fill=NA)
    )
  
  
  # --- 6. Assembly ---
  # P1 (Main) -> P2 (Strip) -> P3 (Med) -> P4 (Small) -> P5 (Small) -> P6 (Small)
  combined_plot <- (p1 / p2 / p3 / p4 / p5 / p6) + 
    patchwork::plot_layout(heights = c(3, 0.4, 1.2, 0.6, 0.6, 0.6))
  
  return(combined_plot)
}
