#!/usr/bin/env Rscript

# Render paper-scope figures from the CSV tables produced by the
# complementarity analysis scripts. Base graphics are used so figure rendering
# does not depend on an additional plotting library. The input directory is
# supplied through GENESELECTR_COMPLEMENTARITY_RESULTS_DIR.

results_dir <- Sys.getenv(
  "GENESELECTR_COMPLEMENTARITY_RESULTS_DIR",
  unset = file.path("redesign", "complementarity_analysis", "results")
)
figure_dir <- file.path(results_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

paper_datasets <- c(
  "GSE101794", "GSE107994", "GSE13355", "GSE65682", "GSE69683",
  "imvigor210", "sosall"
)
method_order <- c(
  "GeneSelectR", "DGE", "Random forest", "Boruta", "mRMR", "LASSO",
  "Elastic net"
)
method_colours <- c(
  "GeneSelectR" = "#168C82", "DGE" = "#D95F4A",
  "Random forest" = "#3478B8", "Boruta" = "#6A5BB5",
  "mRMR" = "#C79216", "LASSO" = "#66717A", "Elastic net" = "#9A6A4F"
)
group_colours <- c(
  "SHARED" = "#334E68", "DGE_ONLY" = "#D95F4A", "GS_ONLY" = "#168C82"
)
dataset_colours <- setNames(
  grDevices::hcl.colors(length(paper_datasets), palette = "Dark 3"),
  paper_datasets
)
heatmap_colours <- grDevices::colorRampPalette(
  c("#440154", "#31688E", "#35B779", "#FDE725")
)(100)

read_table <- function(filename) {
  path <- file.path(results_dir, filename)
  if (!file.exists(path)) {
    stop("Missing figure input: ", path, call. = FALSE)
  }
  table <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if ("dataset" %in% names(table)) {
    extra <- setdiff(unique(table$dataset), paper_datasets)
    if (length(extra) > 0L) {
      stop("Figure input contains datasets outside paper scope: ",
           paste(extra, collapse = ", "), call. = FALSE)
    }
  }
  table
}

save_base_figure <- function(stem, width, height, draw) {
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  grDevices::pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  draw()
  grDevices::dev.off()
  pdftoppm <- Sys.getenv("PDFTOPPM_BIN", unset = Sys.which("pdftoppm"))
  if (!nzchar(pdftoppm) || !file.exists(pdftoppm)) {
    stop("pdftoppm is required to render the PNG preview", call. = FALSE)
  }
  png_prefix <- file.path(figure_dir, paste0(stem, ".png_tmp"))
  status <- system2(
    pdftoppm,
    c("-png", "-r", "300", shQuote(pdf_path), shQuote(png_prefix)),
    stdout = FALSE, stderr = FALSE
  )
  png_generated <- paste0(png_prefix, "-1.png")
  if (!identical(status, 0L) || !file.exists(png_generated) ||
      !file.rename(png_generated, png_path)) {
    stop("PNG preview was not created: ", stem, call. = FALSE)
  }
  created <- c(pdf_path, png_path)
  if (any(!file.exists(created)) || any(file.info(created)$size == 0L)) {
    stop("Figure output was not created: ", stem, call. = FALSE)
  }
}

assert_paper_rows <- function(table, filename) {
  if (!"dataset" %in% names(table)) return(invisible(table))
  missing <- setdiff(paper_datasets, unique(table$dataset))
  if (length(missing) > 0L) {
    stop(filename, " is missing paper dataset(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(table)
}

pretty_dataset <- function(x) {
  x <- sub("^GSE", "GSE", x)
  x
}

draw_bin2d <- function(x, y, title, xlab = "", ylab = "", y_limits = NULL,
                       show_x = TRUE, show_y = TRUE) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) == 0L) {
    plot.new()
    title(main = title, cex.main = 0.85)
    return(invisible(NULL))
  }
  xbreaks <- seq(0, 1, length.out = 25)
  if (is.null(y_limits)) y_limits <- range(y, finite = TRUE)
  if (diff(y_limits) == 0) y_limits <- y_limits + c(-0.5, 0.5)
  ybreaks <- seq(y_limits[[1L]], y_limits[[2L]], length.out = 25)
  xbin <- cut(x, breaks = xbreaks, include.lowest = TRUE, labels = FALSE)
  ybin <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
  counts <- table(
    factor(xbin, levels = seq_len(24L)),
    factor(ybin, levels = seq_len(24L))
  )
  z <- log1p(matrix(as.numeric(counts), nrow = 24L, ncol = 24L))
  image(
    seq(0, 1, length.out = 24),
    seq(y_limits[[1L]], y_limits[[2L]], length.out = 24),
    z,
    col = heatmap_colours, useRaster = TRUE, axes = FALSE,
    xlim = c(0, 1), ylim = y_limits, xlab = "", ylab = ""
  )
  if (show_x) axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                         cex.axis = 0.65)
  if (show_y) {
    y_ticks <- c(y_limits[[1L]], mean(y_limits), y_limits[[2L]])
    axis(2, at = y_ticks, labels = formatC(y_ticks, format = "fg", digits = 3),
         cex.axis = 0.65, las = 1)
  }
  box()
  title(main = title, xlab = xlab, ylab = ylab, cex.main = 0.78,
        cex.lab = 0.72)
  invisible(NULL)
}

# Figure 2 and stability supplements -----------------------------------------
stability <- read_table("gene_set_stability_summary.csv")
assert_paper_rows(stability, "gene_set_stability_summary.csv")
stability$method <- factor(stability$method, levels = method_order)
stability$dataset <- factor(stability$dataset, levels = paper_datasets)

save_base_figure("Figure_2_prediction_and_stability", 11, 8.5, function() {
  layout(rbind(c(1, 2, 3), c(4, 4, 4)),
         heights = c(1.2, 1), widths = c(1, 1, 1))
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mar = c(3, 3.5, 2.2, 0.8), oma = c(3.2, 0, 0, 0))
  for (k_value in c(10L, 20L, 50L)) {
    subset <- stability[stability$k == k_value, , drop = FALSE]
    plot(NA, xlim = c(0, 1), ylim = range(subset$mean_auc, finite = TRUE),
         xlab = "Nogueira stability", ylab = "Mean test AUC",
         main = paste0("A  k = ", k_value), cex.main = 0.9,
         cex.axis = 0.7, cex.lab = 0.78)
    for (method in method_order) {
      rows <- subset[subset$method == method, , drop = FALSE]
      points(rows$nogueira_stability, rows$mean_auc,
             pch = 16, col = method_colours[[method]], cex = 0.7)
    }
    grid(col = "#E4E9ED", lty = 1)
    box()
  }
  means <- aggregate(
    subset(stability,
           select = c(nogueira_stability, mean_auc_minus_random)),
    by = list(dataset_group = stability$dataset_group,
              method = stability$method, k = stability$k), FUN = mean
  )
  group <- as.character(unique(means$dataset_group)[1L])
  subset <- means[means$dataset_group == group, , drop = FALSE]
  xlim <- range(subset$nogueira_stability, finite = TRUE)
  ylim <- range(subset$mean_auc_minus_random, finite = TRUE)
  plot(NA, xlim = xlim + c(-0.03, 0.03), ylim = ylim + c(-0.02, 0.02),
       xlab = "Mean Nogueira stability", ylab = "Mean AUC minus random",
       main = paste0("B  ", group, " dataset-group means"),
       cex.main = 0.9, cex.axis = 0.7, cex.lab = 0.78)
  abline(h = 0, col = "#788995", lwd = 1)
  for (method in method_order) {
    rows <- subset[subset$method == method, , drop = FALSE]
    rows <- rows[order(rows$k), , drop = FALSE]
    lines(rows$nogueira_stability, rows$mean_auc_minus_random,
          col = method_colours[[method]], lwd = 1)
    points(rows$nogueira_stability, rows$mean_auc_minus_random,
           col = method_colours[[method]], pch = 16, cex = 0.8)
  }
  grid(col = "#E4E9ED", lty = 1)
  box()
  par(xpd = NA)
  legend("bottom", inset = c(0, -0.18), legend = method_order,
         col = unname(method_colours[method_order]), pch = 16,
         horiz = TRUE, bty = "n", cex = 0.72)
})

save_base_figure("Supplementary_stability_by_gene_set_size", 10, 8, function() {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(3, 3), mar = c(3, 3.5, 2, 0.6), oma = c(3, 0, 1, 0))
  for (method in method_order) {
    subset <- stability[stability$method == method, , drop = FALSE]
    ylim <- range(subset$nogueira_stability, finite = TRUE)
    plot(NA, xlim = c(10, 50), ylim = ylim + c(-0.03, 0.03),
         xlab = "Gene-set size", ylab = "Nogueira stability",
         main = method, xaxt = "n", cex.main = 0.85,
         cex.axis = 0.72, cex.lab = 0.78)
    axis(1, c(10, 20, 50), c("10", "20", "50"), cex.axis = 0.72)
    for (dataset in paper_datasets) {
      rows <- subset[subset$dataset == dataset, , drop = FALSE]
      rows <- rows[order(rows$k), , drop = FALSE]
      lines(rows$k, rows$nogueira_stability,
            col = dataset_colours[[dataset]], lwd = 1)
      points(rows$k, rows$nogueira_stability,
             col = dataset_colours[[dataset]], pch = 16, cex = 0.65)
    }
    grid(col = "#E4E9ED", lty = 1)
    box()
  }
  plot.new()
  par(xpd = NA)
  legend("center", legend = paper_datasets,
         col = unname(dataset_colours[paper_datasets]), lwd = 1, pch = 16,
         ncol = 2, bty = "n", cex = 0.8)
  mtext("Each line represents one manuscript dataset", side = 1,
        outer = TRUE, line = 1, cex = 0.8)
})

save_base_figure("Supplementary_stability_heatmap_k20", 10, 5, function() {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  subset <- stability[stability$k == 20L, , drop = FALSE]
  z <- matrix(NA_real_, nrow = length(method_order), ncol = length(paper_datasets),
              dimnames = list(method_order, paper_datasets))
  for (i in seq_len(nrow(subset))) {
    z[as.character(subset$method[i]), as.character(subset$dataset[i])] <-
      subset$nogueira_stability[i]
  }
  par(mar = c(7, 11, 3, 2))
  image(seq_along(paper_datasets), seq_along(method_order), t(z),
        col = heatmap_colours, zlim = c(0, 1), axes = FALSE,
        xlab = "", ylab = "Method")
  axis(1, seq_along(paper_datasets), paper_datasets, las = 2, cex.axis = 0.78)
  axis(2, seq_along(method_order), method_order, las = 2, cex.axis = 0.78)
  for (x in seq_along(paper_datasets)) {
    for (y in seq_along(method_order)) {
      if (is.finite(z[y, x])) text(x, y, sprintf("%.2f", z[y, x]), cex = 0.65)
    }
  }
  box()
  mtext("Dataset", side = 1, line = 5, cex = 0.9)
  title(main = "Gene-set stability at k = 20", cex.main = 0.95)
})

# Figure 3 and DGE-component supplement --------------------------------------
gene_scores <- read_table("dge_geneselectr_gene_scores.csv")
assert_paper_rows(gene_scores, "dge_geneselectr_gene_scores.csv")
gene_scores$dataset <- factor(gene_scores$dataset, levels = paper_datasets)

overlap <- read_table("dge_geneselectr_overlap_summary.csv")
assert_paper_rows(overlap, "dge_geneselectr_overlap_summary.csv")

recurrence <- read_table("gene_group_recurrence_long.csv")
assert_paper_rows(recurrence, "gene_group_recurrence_long.csv")
recurrence <- recurrence[recurrence$k == 20L, , drop = FALSE]
recurrence <- recurrence[order(recurrence$dataset, recurrence$group,
                               -recurrence$n_splits, recurrence$gene), ]
recurrence$within_group_order <- ave(
  seq_len(nrow(recurrence)), interaction(recurrence$dataset, recurrence$group),
  FUN = seq_along
)
recurrence <- recurrence[recurrence$within_group_order <= 2L, , drop = FALSE]

save_base_figure("Figure_3_DGE_GeneSelectR_complementarity", 12, 10, function() {
  mat <- matrix(seq_len(12L), nrow = 3L, byrow = TRUE)
  layout(mat, heights = c(1.3, 1.3, 1.1))
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mar = c(2.5, 3.3, 1.8, 0.5))
  for (dataset_index in seq_along(paper_datasets)) {
    dataset <- paper_datasets[[dataset_index]]
    par(mar = c(2.5, if (dataset_index %in% c(1L, 5L)) 5.8 else 1.2,
                1.8, 0.5))
    rows <- gene_scores[gene_scores$dataset == dataset, , drop = FALSE]
    draw_bin2d(rows$dge_rank_percentile, rows$geneselectr_rank_percentile,
               dataset,
               xlab = if (dataset_index >= 5L) "DGE percentile" else "",
               ylab = if (dataset_index %in% c(1L, 5L))
                 "GeneSelectR percentile" else "",
               y_limits = c(0, 1), show_x = dataset_index >= 5L,
               show_y = dataset_index %in% c(1L, 5L))
  }
  plot.new()
  title(main = "A  Rank percentiles", cex.main = 0.9)

  z <- matrix(NA_real_, nrow = length(paper_datasets), ncol = 3L,
              dimnames = list(paper_datasets, c("10", "20", "50")))
  for (i in seq_len(nrow(overlap))) {
    z[as.character(overlap$dataset[i]), as.character(overlap$k[i])] <-
      overlap$mean_jaccard[i]
  }
  par(mar = c(4.5, 9, 2.2, 2))
  image(seq_len(3L), seq_len(length(paper_datasets)), t(z),
        col = heatmap_colours, zlim = c(0, 1), axes = FALSE,
        xlab = "Gene-set size", ylab = "Dataset")
  axis(1, seq_len(3L), c("10", "20", "50"), cex.axis = 0.78)
  axis(2, seq_along(paper_datasets), paper_datasets, las = 2, cex.axis = 0.7)
  for (x in seq_len(3L)) for (y in seq_along(paper_datasets)) {
    if (is.finite(z[y, x])) text(x, y, sprintf("%.2f", z[y, x]), cex = 0.65)
  }
  box()
  title(main = "B  Mean top-k overlap", cex.main = 0.9)

  rows <- recurrence
  labels <- paste(rows$dataset, rows$gene, sep = ": ")
  y <- rev(seq_len(nrow(rows)))
  plot(NA, xlim = c(0, 15), ylim = c(0.5, length(rows) + 0.5),
       xaxt = "n", yaxt = "n", xlab = "Train/test divisions containing gene",
       ylab = "", main = "C  Most recurrent genes at k = 20",
       cex.main = 0.9, cex.axis = 0.72, cex.lab = 0.78)
  axis(1, c(0, 5, 10, 15), cex.axis = 0.72)
  axis(2, y, labels, las = 2, cex.axis = 0.38)
  for (i in seq_len(nrow(rows))) {
    colour <- group_colours[[as.character(rows$group[i])]]
    segments(0, y[i], rows$n_splits[i], y[i], col = "#CCD5DB")
    points(rows$n_splits[i], y[i], pch = 16, col = colour, cex = 0.75)
  }
  grid(col = "#E4E9ED", lty = 1)
  box()
  legend("bottomright", legend = names(group_colours),
         col = unname(group_colours), pch = 16, bty = "n", cex = 0.68)
})

save_base_figure("Supplementary_DGE_vs_GeneSelectR_components", 16, 6,
                 function() {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  layout(matrix(seq_len(14L), nrow = 2L, byrow = TRUE))
  par(mar = c(3, 3.2, 1.8, 0.3), oma = c(2, 0, 1, 0))
  components <- list(
    "Internal recurrence" = "geneselectr_internal_recurrence",
    "Predictive utility" = "geneselectr_raw_utility"
  )
  panel <- 0L
  for (component_name in names(components)) {
    value_column <- components[[component_name]]
    for (dataset_index in seq_along(paper_datasets)) {
      dataset <- paper_datasets[[dataset_index]]
      panel <- panel + 1L
      par(mar = c(3, if (dataset_index == 1L) 5.8 else 1.2, 1.8, 0.3))
      rows <- gene_scores[gene_scores$dataset == dataset, , drop = FALSE]
      draw_bin2d(rows$dge_rank_percentile, rows[[value_column]], dataset,
                 xlab = if (component_name == "Predictive utility")
                   "DGE percentile" else "",
                 ylab = if (dataset_index == 1L) component_name else "",
                 show_x = component_name == "Predictive utility",
                 show_y = dataset_index == 1L)
    }
  }
  mtext("DGE rank compared with GeneSelectR score components", side = 3,
        outer = TRUE, line = 0, cex = 0.95, font = 2)
})

# Figure 4: biological description of ranking groups -------------------------
biology <- read_table("gene_group_biology_summary.csv")
assert_paper_rows(biology, "gene_group_biology_summary.csv")
biology$ratio_for_plot <- ifelse(
  is.finite(biology$median_observed_to_random_ratio),
  pmax(biology$median_observed_to_random_ratio, 2^-8), NA_real_
)
biology$log2_ratio <- log2(biology$ratio_for_plot)
biology$metric_short <- biology$metric
biology$metric_short <- sub(
  "Open Targets disease-association sum", "Open Targets", biology$metric_short,
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
  "GO semantic mean", "GO similarity", biology$metric_short, fixed = TRUE
)
metric_order <- unique(biology$metric_short)

save_base_figure("Figure_4_gene_group_biology", 12, 10, function() {
  layout(matrix(seq_len(length(metric_order) * 2L),
                nrow = length(metric_order), byrow = TRUE),
         heights = rep(1, length(metric_order)))
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mar = c(6, 3.5, 2, 0.5), oma = c(2.2, 0, 0, 0))
  for (metric in metric_order) {
    for (k_value in c(10L, 20L)) {
      rows <- biology[biology$metric_short == metric & biology$k == k_value,
                      , drop = FALSE]
      finite_values <- rows$log2_ratio[is.finite(rows$log2_ratio)]
      ylim <- if (length(finite_values) > 0L) {
        range(finite_values)
      } else {
        c(-1, 1)
      }
      if (diff(ylim) == 0) ylim <- ylim + c(-1, 1)
      ylim <- ylim + c(-0.08, 0.08) * diff(ylim)
      plot(NA, xlim = c(0.5, length(paper_datasets) + 0.5), ylim = ylim,
           xaxt = "n", xlab = "", ylab = "log2(observed / random)",
           main = paste0(metric, ", k = ", k_value),
           cex.main = 0.82, cex.axis = 0.68, cex.lab = 0.72)
      axis(1, seq_along(paper_datasets), paper_datasets, las = 2, cex.axis = 0.68)
      abline(h = 0, col = "#788995", lwd = 1)
      for (group in names(group_colours)) {
        group_rows <- rows[rows$group == group, , drop = FALSE]
        positions <- match(group_rows$dataset, paper_datasets) +
          c(SHARED = 0, DGE_ONLY = -0.18, GS_ONLY = 0.18)[[group]]
        points(positions, group_rows$log2_ratio,
               pch = 16, col = group_colours[[group]], cex = 0.75)
      }
      grid(col = "#E4E9ED", lty = 1)
      box()
    }
  }
  par(xpd = NA)
  legend("bottom", inset = c(0, -0.14), legend = names(group_colours),
         col = unname(group_colours), pch = 16, horiz = TRUE,
         bty = "n", cex = 0.75)
  mtext("Biological properties of shared and method-specific genes",
        side = 3, outer = TRUE, line = -1, cex = 0.95, font = 2)
})

message("Wrote paper-scope figures to ", figure_dir)
