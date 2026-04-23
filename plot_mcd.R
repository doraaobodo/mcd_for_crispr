

get_group_col <- function(tbl) {
  if ("group_id" %in% colnames(tbl)) return("group_id")
  stop("No group_id column found.")
}

get_label_genes <- function(tbl, n = 15) {
  tbl <- tbl[order(tbl$p_value, -abs(tbl$mean_effect)), , drop = FALSE]
  tbl$gene[seq_len(min(n, nrow(tbl)))]
}

plot_group_volcano <- function(group_id,
                               outlier_summary_tbl,
                               alpha = 0.05,
                               label_top_n = 15) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  plot_df <- outlier_summary_tbl[outlier_summary_tbl$group_id == group_id, , drop = FALSE]
  if (nrow(plot_df) == 0) stop("No rows found for group_id: ", group_id)
  
  plot_df$neglog10p <- -log10(pmax(plot_df$p_value, .Machine$double.xmin))
  
  top_genes <- get_label_genes(plot_df, n = label_top_n)
  top_tbl <- plot_df[plot_df$gene %in% top_genes, , drop = FALSE]
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = mean_effect, y = neglog10p)) +
    ggplot2::geom_point(ggplot2::aes(color = outlier_flag), alpha = 0.7, size = 1.3) +
    ggplot2::geom_hline(yintercept = -log10(alpha), linetype = 2) +
    ggplot2::labs(
      title = paste("Volcano plot:", group_id),
      x = "Mean contrast effect",
      y = "-log10(p_value)",
      color = paste0("p_value < ", alpha)
    ) +
    ggplot2::theme_bw()
  
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = top_tbl,
      ggplot2::aes(label = gene),
      size = 3,
      max.overlaps = Inf
    )
  }
  
  p
}

write_mcd_report_rmd <- function(file = "mcd_report.Rmd",
                                 title = "MCD Analysis Report") {
  txt <- c(
    '---',
    paste0('title: "', title, '"'),
    'output:',
    '  html_document:',
    '    toc: true',
    '    toc_float: true',
    '    theme: cosmo',
    '---',
    '',
    '```{r setup, include=FALSE}',
    'knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)',
    'library(ggplot2)',
    '```',
    '',
    '# Inputs',
    '',
    '```{r load-data}',
    '# Update these paths as needed',
    'final_input <- readRDS("final_input.rds")',
    'res <- readRDS("mcd_results.rds")',
    '',
    '# source your plotting/helper functions',
    'source("mcd_report_plots.R")',
    '```',
    '',
    '# Run summary',
    '',
    '```{r run-summary}',
    'summary_df <- data.frame(',
    '  field = c("input_type", "selected_metric", "group_by", "contrast_mode", "control_treatment"),',
    '  value = c(',
    '    final_input$input_type,',
    '    final_input$selected_metric,',
    '    paste(res$group_by, collapse = ";"),',
    '    res$contrast_mode,',
    '    res$control_treatment',
    '  ),',
    '  stringsAsFactors = FALSE',
    ')',
    'knitr::kable(summary_df)',
    '```',
    '',
    '# Metadata overview',
    '',
    '```{r metadata-overview}',
    'knitr::kable(head(final_input$metadata, 20))',
    '```',
    '',
    '# Outlier counts by group',
    '',
    '```{r outlier-counts, fig.width=8, fig.height=4}',
    'plot_outlier_counts_by_group(res, alpha = 0.05)',
    '```',
    '',
    '# Heatmap across groups',
    '',
    '```{r cross-group-heatmap, fig.width=8, fig.height=10}',
    'plot_outlier_heatmap_across_groups(',
    '  res = res,',
    '  contrast_matrices = res$contrast_matrices,',
    '  alpha = 0.05,',
    '  n = 75',
    ')',
    '```',
    '',
    '# Per-group results',
    '',
    '```{r per-group-loop, results="asis"}',
    'for (g in names(res$mcd_results_by_group)) {',
    '  cat("\\n\\n## ", g, "\\n\\n", sep = "")',
    '  cat("Number of contrast columns: ", ncol(res$contrast_matrices[[g]]), "\\n\\n", sep = "")',
    '  ',
    '  print(plot_group_volcano(',
    '    group_id = g,',
    '    res = res,',
    '    contrast_matrices = res$contrast_matrices,',
    '    alpha = 0.05',
    '  ))',
    '  ',
    '  print(plot_group_mcd_scores(g, res, alpha = 0.05))',
    '  ',
    '  plot_group_outlier_heatmap(',
    '    group_id = g,',
    '    res = res,',
    '    contrast_matrices = res$contrast_matrices,',
    '    alpha = 0.05,',
    '    n = 50',
    '  )',
    '  ',
    '  if (ncol(res$contrast_matrices[[g]]) >= 2) {',
    '    print(plot_group_pca_genes(',
    '      group_id = g,',
    '      res = res,',
    '      contrast_matrices = res$contrast_matrices,',
    '      alpha = 0.05',
    '    ))',
    '    ',
    '    if (requireNamespace("uwot", quietly = TRUE)) {',
    '      print(plot_group_umap_genes(',
    '        group_id = g,',
    '        res = res,',
    '        contrast_matrices = res$contrast_matrices,',
    '        alpha = 0.05',
    '      ))',
    '    }',
    '  }',
    '}',
    '```'
  )
  
  writeLines(txt, con = file)
  invisible(file)
}