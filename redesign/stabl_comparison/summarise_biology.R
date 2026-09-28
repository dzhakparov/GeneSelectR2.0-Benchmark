#!/usr/bin/env Rscript

repo <- normalizePath(getwd(), mustWork = TRUE)
out <- Sys.getenv("STABL_BIOLOGY_OUT",
                  file.path(repo, "redesign/stabl_comparison/biology"))
datasets <- c("GSE101794", "GSE107994", "GSE13355", "GSE65682",
              "GSE69683", "imvigor210", "sosall")
paths <- file.path(out, paste0(datasets, "_stabl_biology.csv"))
if (!all(file.exists(paths))) stop("All seven dataset files are required")
detail <- do.call(rbind, lapply(paths, read.csv, check.names = FALSE))
metrics <- c("GO_enrichment", "hallmark_enrichment", "string_enrichment",
             "ot_0.05_enrichment", "ot_0.10_enrichment")
keys <- with(detail, paste(dataset, repeat_idx, fold_idx, k))
stopifnot(nrow(detail) == 630L, !anyDuplicated(keys),
          all(is.finite(as.matrix(detail[metrics]))),
          all(detail$pool_size == 2000L), all(detail$n_null == 1000L))
by_size <- aggregate(detail[metrics], detail[c("dataset", "k")], median)
stopifnot(nrow(by_size) == 42L)
primary <- aggregate(by_size[by_size$k %in% c(10L, 20L, 50L), metrics],
                     by_size[by_size$k %in% c(10L, 20L, 50L), "dataset",
                             drop = FALSE], mean)
stopifnot(nrow(primary) == 7L)
write.csv(detail, file.path(out, "stabl_biology_by_split.csv"), row.names = FALSE)
write.csv(by_size, file.path(out, "stabl_biology_by_size.csv"), row.names = FALSE)
write.csv(primary, file.path(out, "stabl_biology_by_dataset.csv"), row.names = FALSE)
cat("Assembled 630 split-size rows, 42 size medians, and seven dataset summaries.\n")
