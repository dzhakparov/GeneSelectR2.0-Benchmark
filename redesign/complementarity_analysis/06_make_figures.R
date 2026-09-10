#!/usr/bin/env Rscript

# Create manuscript figures from saved CSV tables produced by scripts 01-04.
# No statistics are recalculated from rankings or expression matrices here.

source(file.path("redesign", "complementarity_analysis", "00_utils.R"))
require_packages(c("ggplot2", "patchwork"))

read_result <- function(filename, required_columns) {
  path <- file.path(complementarity_results_dir, filename)
  if (!file.exists(path)) {
    stop("Missing figure input: ", path, call. = FALSE)
  }
  table <- utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE
  )
  missing <- setdiff(required_columns, names(table))
  if (length(missing) > 0L) {
    stop(filename, " is missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  table
}

source_dir <- file.path(complementarity_results_dir, "figure_source_data")
dir.create(complementarity_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

method_colours <- c(
  "GeneSelectR" = "#168C82",
  "DGE" = "#D95F4A",
  "Random forest" = "#3478B8",
  "Boruta" = "#6A5BB5",
  "mRMR" = "#C79216",
  "LASSO" = "#66717A",
  "Elastic net" = "#9A6A4F"
)
group_colours <- c(
  "SHARED" = "#334E68",
  "DGE_ONLY" = "#D95F4A",
  "GS_ONLY" = "#168C82"
)

paper_theme <- ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = "#E4E9ED", linewidth = 0.3),
    strip.background = ggplot2::element_rect(fill = "#EEF2F5", colour = NA),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold", colour = "#17365D"),
    legend.position = "bottom"
  )

save_figure <- function(plot, stem, width, height) {
  pdf_path <- file.path(
    complementarity_figure_dir, paste0(stem, ".pdf")
  )
  png_path <- file.path(
    complementarity_figure_dir, paste0(stem, ".png")
  )
  ggplot2::ggsave(
    pdf_path, plot, width = width, height = height, units = "in",
    device = grDevices::pdf
  )
  ggplot2::ggsave(
    png_path, plot, width = width, height = height, units = "in", dpi = 300
  )
  created <- c(pdf_path, png_path)
  if (any(!file.exists(created)) || any(file.info(created)$size == 0L)) {
    stop("Figure output was not created: ", stem, call. = FALSE)
  }
}

# Figure 2: prediction and stability ------------------------------------------
stability <- read_result(
  "gene_set_stability_summary.csv",
  c(
    "dataset", "dataset_group", "method", "k",
    "nogueira_stability", "nogueira_exact_topk_union",
    "nogueira_fixed_common_genes_actual_selection",
    "nogueira_fixed_common_genes_reranked_topk",
    "mean_jaccard", "mean_auc", "mean_auc_minus_random"
  )
)
stability$k_label <- paste0("k = ", stability$k)
atomic_write_csv(
  stability,
  file.path(source_dir, "Figure_2A_dataset_prediction_stability.csv")
)

figure_2a <- ggplot2::ggplot(
  stability,
  ggplot2::aes(
    x = nogueira_stability, y = mean_auc,
    colour = method, shape = dataset_group
  )
) +
  ggplot2::geom_point(size = 2.0, alpha = 0.8) +
  ggplot2::facet_wrap(~k_label, nrow = 1) +
  ggplot2::scale_colour_manual(values = method_colours) +
  ggplot2::labs(
    title = "A  Prediction accuracy and cross-split gene-set stability",
    x = "Nogueira stability among genes available in every division",
    y = "Mean test AUC", colour = "Method", shape = "Dataset group"
  ) + paper_theme

mean_stability <- stats::aggregate(
  stability[, c(
    "nogueira_stability", "mean_auc", "mean_auc_minus_random"
  )],
  by = stability[, c("dataset_group", "method", "k")],
  FUN = mean
)
mean_stability$k_label <- paste0("k = ", mean_stability$k)
atomic_write_csv(
  mean_stability,
  file.path(source_dir, "Figure_2B_method_means.csv")
)

figure_2b <- ggplot2::ggplot(
  mean_stability,
  ggplot2::aes(
    x = nogueira_stability, y = mean_auc_minus_random,
    colour = method, group = method
  )
) +
  ggplot2::geom_hline(yintercept = 0, colour = "#788995", linewidth = 0.35) +
  ggplot2::geom_path(alpha = 0.55) +
  ggplot2::geom_point(
    ggplot2::aes(shape = k_label), size = 2.3, stroke = 0.7
  ) +
  ggplot2::facet_wrap(~dataset_group, scales = "free_y") +
  ggplot2::scale_colour_manual(values = method_colours) +
  ggplot2::labs(
    title = "B  Dataset-group means for each method",
    x = "Mean Nogueira stability",
    y = "Mean AUC minus same-size random genes",
    colour = "Method", shape = "Gene-set size"
  ) + paper_theme

figure_2 <- figure_2a / figure_2b +
  patchwork::plot_layout(heights = c(1.05, 1))
save_figure(figure_2, "Figure_2_prediction_and_stability", 11, 8.5)

stability_curve_source <- stability[, c(
  "dataset", "dataset_group", "method", "k",
  "nogueira_stability", "nogueira_exact_topk_union",
  "nogueira_fixed_common_genes_actual_selection",
  "nogueira_fixed_common_genes_reranked_topk",
  "mean_jaccard"
)]
atomic_write_csv(
  stability_curve_source,
  file.path(source_dir, "Supplementary_stability_by_gene_set_size.csv")
)
stability_curve <- ggplot2::ggplot(
  stability_curve_source,
  ggplot2::aes(
    k, nogueira_stability, group = dataset, colour = dataset_group
  )
) +
  ggplot2::geom_line(alpha = 0.45) +
  ggplot2::geom_point(size = 1.2) +
  ggplot2::facet_wrap(~method, ncol = 3) +
  ggplot2::scale_x_continuous(breaks = primary_gene_set_sizes) +
  ggplot2::labs(
    title = "Cross-split stability across gene-set sizes",
    x = "Number of genes", y = "Nogueira stability",
    colour = "Dataset group"
  ) + paper_theme
save_figure(
  stability_curve, "Supplementary_stability_by_gene_set_size", 10, 8
)

stability_heatmap_source <- stability[
  stability$k == 20L,
  c(
    "dataset", "dataset_group", "method", "k",
    "nogueira_stability"
  )
]
atomic_write_csv(
  stability_heatmap_source,
  file.path(source_dir, "Supplementary_stability_heatmap_k20.csv")
)
stability_heatmap <- ggplot2::ggplot(
  stability_heatmap_source,
  ggplot2::aes(dataset, method, fill = nogueira_stability)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
  ggplot2::facet_grid(~dataset_group, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_viridis_c(limits = c(0, 1)) +
  ggplot2::labs(
    title = "Cross-split stability for 20-gene sets",
    x = NULL, y = NULL, fill = "Nogueira"
  ) + paper_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(
    angle = 45, hjust = 1, vjust = 1
  ))
save_figure(stability_heatmap, "Supplementary_stability_heatmap_k20", 10, 5)

# Figure 3: DGE and GeneSelectR complementarity -------------------------------
gene_scores <- read_result(
  "dge_geneselectr_gene_scores.csv",
  c(
    "dataset", "dataset_group", "repeat_idx", "fold_idx", "gene",
    "geneselectr_rank_percentile", "dge_rank_percentile",
    "geneselectr_internal_recurrence", "geneselectr_raw_utility"
  )
)
figure_3a_source <- gene_scores[, c(
  "dataset", "dataset_group", "repeat_idx", "fold_idx", "gene",
  "geneselectr_rank_percentile", "dge_rank_percentile"
)]
atomic_write_csv(
  figure_3a_source,
  file.path(source_dir, "Figure_3A_rank_percentiles.csv")
)

figure_3a <- ggplot2::ggplot(
  figure_3a_source,
  ggplot2::aes(dge_rank_percentile, geneselectr_rank_percentile)
) +
  ggplot2::geom_bin_2d(bins = 24) +
  ggplot2::facet_wrap(~dataset, ncol = 4) +
  ggplot2::scale_fill_viridis_c(trans = "log10") +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "A  DGE and GeneSelectR rank percentiles",
    x = "DGE rank percentile (higher = stronger)",
    y = "GeneSelectR rank percentile (higher = stronger)",
    fill = "Genes"
  ) + paper_theme

overlap <- read_result(
  "dge_geneselectr_overlap_summary.csv",
  c("dataset", "dataset_group", "k", "mean_jaccard")
)
overlap$k_label <- paste0("k = ", overlap$k)
atomic_write_csv(
  overlap,
  file.path(source_dir, "Figure_3B_topk_overlap.csv")
)

figure_3b <- ggplot2::ggplot(
  overlap,
  ggplot2::aes(k_label, dataset, fill = mean_jaccard)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
  ggplot2::facet_grid(dataset_group ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_fill_viridis_c(limits = c(0, 1)) +
  ggplot2::labs(
    title = "B  Mean overlap of the top-ranked genes",
    x = "Gene-set size", y = NULL, fill = "Jaccard"
  ) + paper_theme

recurrence <- read_result(
  "gene_group_recurrence_long.csv",
  c("dataset", "dataset_group", "k", "gene", "group", "n_splits")
)
recurrence <- recurrence[recurrence$k == 20L, , drop = FALSE]
recurrence <- recurrence[order(
  recurrence$dataset, recurrence$group, -recurrence$n_splits,
  recurrence$gene
), ]
recurrence$within_group_order <- ave(
  seq_len(nrow(recurrence)),
  interaction(recurrence$dataset, recurrence$group, drop = TRUE),
  FUN = seq_along
)
figure_3c_source <- recurrence[recurrence$within_group_order <= 2L, ]
atomic_write_csv(
  figure_3c_source,
  file.path(source_dir, "Figure_3C_recurrent_gene_examples.csv")
)

figure_3c_source$label <- paste0(
  figure_3c_source$dataset, ": ", figure_3c_source$gene
)
figure_3c <- ggplot2::ggplot(
  figure_3c_source,
  ggplot2::aes(n_splits, reorder(label, n_splits), colour = group)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = n_splits, yend = reorder(label, n_splits)),
    colour = "#CCD5DB", linewidth = 0.5
  ) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~group, scales = "free_y") +
  ggplot2::scale_colour_manual(values = group_colours) +
  ggplot2::scale_x_continuous(limits = c(0, 15), breaks = c(0, 5, 10, 15)) +
  ggplot2::labs(
    title = "C  Most recurrent genes under a fixed selection rule (k = 20)",
    x = "Train/test divisions containing the gene", y = NULL,
    colour = "Ranking group"
  ) + paper_theme

figure_3 <- figure_3a / (figure_3b | figure_3c) +
  patchwork::plot_layout(heights = c(1.15, 1))
save_figure(figure_3, "Figure_3_DGE_GeneSelectR_complementarity", 12, 10)

component_source <- rbind(
  data.frame(
    figure_3a_source[, c(
      "dataset", "dataset_group", "repeat_idx", "fold_idx", "gene",
      "dge_rank_percentile"
    )],
    component = "GeneSelectR internal recurrence",
    component_value = gene_scores$geneselectr_internal_recurrence
  ),
  data.frame(
    figure_3a_source[, c(
      "dataset", "dataset_group", "repeat_idx", "fold_idx", "gene",
      "dge_rank_percentile"
    )],
    component = "GeneSelectR predictive utility",
    component_value = gene_scores$geneselectr_raw_utility
  )
)
atomic_write_csv(
  component_source,
  file.path(source_dir, "Supplementary_DGE_vs_GeneSelectR_components.csv")
)
component_plot <- ggplot2::ggplot(
  component_source,
  ggplot2::aes(dge_rank_percentile, component_value)
) +
  ggplot2::geom_bin_2d(bins = 24) +
  ggplot2::facet_grid(component ~ dataset, scales = "free_y") +
  ggplot2::scale_fill_viridis_c(trans = "log10") +
  ggplot2::labs(
    title = "DGE rank compared with GeneSelectR score components",
    x = "DGE rank percentile (higher = stronger)", y = "Component value",
    fill = "Genes"
  ) + paper_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(
    angle = 45, hjust = 1, vjust = 1
  ))
save_figure(
  component_plot, "Supplementary_DGE_vs_GeneSelectR_components", 16, 6
)

# Figure 4: biological description of ranking groups -------------------------
biology <- read_result(
  "gene_group_biology_summary.csv",
  c(
    "dataset", "dataset_group", "k", "group", "metric",
    "median_group_size", "median_observed_to_random_ratio", "n_missing"
  )
)
biology$ratio_for_plot <- ifelse(
  is.finite(biology$median_observed_to_random_ratio),
  pmax(biology$median_observed_to_random_ratio, 2^-8), NA_real_
)
biology$log2_ratio <- log2(biology$ratio_for_plot)
biology$k_label <- paste0("k = ", biology$k)
biology$metric_short <- sub(
  "Open Targets disease-association sum", "Open Targets", biology$metric,
  fixed = TRUE
)
biology$metric_short <- sub(
  "Hallmark shared-membership edges", "Hallmark", biology$metric_short,
  fixed = TRUE
)
biology$metric_short <- sub(
  "STRING protein-association edges", "STRING", biology$metric_short,
  fixed = TRUE
)
biology$metric_short <- sub(
  "GO semantic mean", "GO similarity", biology$metric_short,
  fixed = TRUE
)
atomic_write_csv(
  biology,
  file.path(source_dir, "Figure_4_group_biology.csv")
)

figure_4 <- ggplot2::ggplot(
  biology,
  ggplot2::aes(dataset, log2_ratio, colour = group, shape = dataset_group)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "#788995", linewidth = 0.4) +
  ggplot2::geom_point(
    position = ggplot2::position_dodge(width = 0.55), size = 1.9,
    alpha = 0.85
  ) +
  ggplot2::facet_grid(metric_short ~ k_label, scales = "free_y") +
  ggplot2::scale_colour_manual(values = group_colours) +
  ggplot2::labs(
    title = "Biological properties of shared and method-specific genes",
    subtitle = "Each value is compared with same-size random genes from the same training-specific candidate set",
    x = NULL, y = "log2(observed / random-set mean)",
    colour = "Ranking group", shape = "Dataset group"
  ) +
  paper_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(
    angle = 45, hjust = 1, vjust = 1
  ))
save_figure(figure_4, "Figure_4_gene_group_biology", 12, 10)

message("Wrote figures and source-data tables to ",
        complementarity_results_dir)
