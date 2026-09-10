#!/usr/bin/env Rscript

source(file.path("analysis", "config.R"))

required_packages <- c(
    "glmnet", "ranger", "xgboost", "withr", "pROC", "edgeR",
    "SummarizedExperiment", "AnnotationDbi", "GO.db", "org.Hs.eg.db",
    "msigdbr", "igraph"
)
package_status <- vapply(required_packages, requireNamespace, logical(1),
                         quietly = TRUE)

required_paths <- c(
    package_source = file.path("package", "GeneSelectR", "R"),
    sosall_expression = file.path("data", "normalized_logcpm.csv"),
    sosall_metadata = file.path("data", "metadata.csv"),
    imvigor210 = file.path("data", "IMvigor210.all.rds"),
    validation_registry = file.path("benchmarks", "validation_datasets.R"),
    gse65682_expression = file.path(
        "data", "GSE65682", "expression_prepared.csv"
    ),
    gse65682_metadata = file.path(
        "data", "GSE65682", "metadata_prepared.csv"
    ),
    gse69683_expression = file.path(
        "data", "GSE69683", "expression_prepared.csv"
    ),
    gse69683_metadata = file.path(
        "data", "GSE69683", "metadata_prepared.csv"
    ),
    gse13355_expression = file.path(
        "data", "GSE13355", "expression_prepared.csv"
    ),
    gse13355_metadata = file.path(
        "data", "GSE13355", "metadata_prepared.csv"
    ),
    gse107994_expression = file.path(
        "data", "GSE107994", "counts_prepared.csv"
    ),
    gse107994_metadata = file.path(
        "data", "GSE107994", "metadata_prepared.csv"
    ),
    gse101794_expression = file.path("data", "GSE101794", "expression.tsv"),
    gse101794_metadata = file.path("data", "GSE101794", "metadata.tsv")
)
path_status <- file.exists(required_paths)

cat("Software packages\n")
print(data.frame(package = names(package_status), available = package_status),
      row.names = FALSE)
cat("\nRequired local inputs\n")
print(data.frame(path = unname(required_paths), available = path_status),
      row.names = FALSE)

if (!all(package_status)) {
    cat("\nMissing packages:\n")
    cat(paste(names(package_status)[!package_status], collapse = ", "), "\n")
}
if (!all(path_status)) {
    cat("\nMissing local inputs:\n")
    cat(paste(unname(required_paths[!path_status]), collapse = "\n"), "\n")
}

if (!all(package_status) || !all(path_status)) {
  quit(save = "no", status = 1L)
}

invisible(TRUE)
