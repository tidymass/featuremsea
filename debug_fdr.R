# 调试 FDR 分布的代码
# 在运行 perform_fmsea_analysis 后使用此代码

debug_fdr_distribution <- function(result_object) {
  # 获取所有计算过的结果
  res_list <- result_object@res_list

  # 提取所有通路的 p_value 和 FDR
  all_results <- purrr::imap_dfr(
    res_list,
    ~ tibble::tibble(
      pathway_id = .y,
      ES = .x$ES %||% NA_real_,
      NES = .x$NES %||% NA_real_,
      p_value = .x$p_value %||% NA_real_
    )
  )

  # 计算 FDR
  all_results$FDR <- NA_real_
  idx <- is.finite(all_results$p_value) & !is.na(all_results$p_value)
  if (any(idx)) {
    all_results$FDR[idx] <- p.adjust(all_results$p_value[idx], method = "BH")
  }

  # 移除 NA 值
  valid_results <- all_results[!is.na(all_results$FDR), ]

  cat("=== FDR 分布统计 ===\n")
  cat("总通路数:", nrow(all_results), "\n")
  cat("有效 FDR 计算数:", nrow(valid_results), "\n")

  if (nrow(valid_results) > 0) {
    cat("\nFDR 分布:\n")
    print(summary(valid_results$FDR))

    cat("\n不同阈值下的显著通路数:\n")
    thresholds <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
    for (thr in thresholds) {
      count <- sum(valid_results$FDR < thr, na.rm = TRUE)
      cat(sprintf("FDR < %.2f: %d 个通路\n", thr, count))
    }

    cat("\n最显著的10个通路:\n")
    top10 <- valid_results[order(valid_results$FDR)[1:min(10, nrow(valid_results))], ]
    print(top10[, c("pathway_id", "p_value", "FDR")])
  }

  return(valid_results)
}

# 使用方法:
# result <- perform_fmsea_analysis(...)
# fdr_stats <- debug_fdr_distribution(result)