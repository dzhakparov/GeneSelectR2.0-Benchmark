# Index historical benchmark outputs after independent conclusions are frozen.
#
# This script records every file and summarizes each run directory. It does not
# select favorable runs. Comparability decisions are made from recorded
# configuration and integrity fields in the resulting index.

suppressPackageStartupMessages(library(dplyr))

project_root <- normalizePath(".", mustWork = TRUE)
historical_roots <- file.path(
  project_root,
  c("results_sosall", "results_imvigor210", "results_validation")
)
historical_roots <- historical_roots[dir.exists(historical_roots)]

output_root <- file.path(
  project_root,
  "independent_benchmark_runs",
  "historical_comparison"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

all_files <- unlist(lapply(historical_roots, function(root) {
  list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE,
             no.. = TRUE)
}), use.names = FALSE)
all_files <- all_files[file.info(all_files)$isdir %in% FALSE]

file_info <- file.info(all_files)
file_index <- data.frame(
  path = substring(all_files, nchar(project_root) + 2L),
  bytes = file_info$size,
  modified = format(file_info$mtime, tz = "Asia/Almaty", usetz = TRUE),
  md5 = unname(tools::md5sum(all_files)),
  stringsAsFactors = FALSE
)
write.csv(file_index, file.path(output_root, "historical_file_index.csv"),
          row.names = FALSE)

top_level_runs <- c(
  list.dirs(file.path(project_root, "results_sosall"),
            recursive = FALSE, full.names = TRUE),
  list.dirs(file.path(project_root, "results_imvigor210"),
            recursive = FALSE, full.names = TRUE)
)
validation_datasets <- list.dirs(
  file.path(project_root, "results_validation"),
  recursive = FALSE,
  full.names = TRUE
)
validation_runs <- unlist(lapply(validation_datasets, function(root) {
  list.dirs(root, recursive = FALSE, full.names = TRUE)
}), use.names = FALSE)
run_dirs <- sort(unique(c(top_level_runs, validation_runs)))

scalar_or_na <- function(x, name, default = NA) {
  value <- x[[name]]
  if (is.null(value) || length(value) == 0L) default else value[[1L]]
}

identify_dataset <- function(run_dir) {
  relative <- substring(run_dir, nchar(project_root) + 2L)
  parts <- strsplit(relative, "/", fixed = TRUE)[[1L]]
  if (parts[1L] == "results_sosall") return("SOS-ALL")
  if (parts[1L] == "results_imvigor210") return("IMvigor210")
  if (parts[1L] == "results_validation" && length(parts) >= 2L) {
    return(parts[2L])
  }
  "unknown"
}

summarise_run <- function(run_dir) {
  data_dir <- file.path(run_dir, "data")
  nested_file <- file.path(data_dir, "nested_results.csv")
  config_file <- file.path(data_dir, "config.rds")
  relative <- substring(run_dir, nchar(project_root) + 2L)

  nested <- if (file.exists(nested_file)) {
    tryCatch(read.csv(nested_file, stringsAsFactors = FALSE),
             error = function(e) NULL)
  } else {
    NULL
  }
  config <- if (file.exists(config_file)) {
    tryCatch(suppressWarnings(readRDS(config_file)), error = function(e) NULL)
  } else {
    NULL
  }

  methods <- if (!is.null(nested) && "Method" %in% names(nested)) {
    length(unique(nested$Method))
  } else {
    NA_integer_
  }
  evaluators <- if (!is.null(nested) && "Evaluator" %in% names(nested)) {
    length(unique(nested$Evaluator))
  } else if (!is.null(nested)) {
    1L
  } else {
    NA_integer_
  }
  split_columns <- c("Repeat", "Fold")
  splits <- if (!is.null(nested) && all(split_columns %in% names(nested))) {
    length(unique(paste(nested$Repeat, nested$Fold, sep = ":")))
  } else {
    NA_integer_
  }
  panel_sizes <- if (!is.null(nested) && "k" %in% names(nested)) {
    sort(unique(nested$k))
  } else {
    numeric(0)
  }
  expected_rows <- if (all(is.finite(c(methods, evaluators, splits))) &&
                       length(panel_sizes) > 0L) {
    methods * evaluators * splits * length(panel_sizes)
  } else {
    NA_integer_
  }

  key_columns <- intersect(
    c("Method", "Evaluator", "Repeat", "Fold", "k"),
    names(nested)
  )
  duplicate_keys <- if (!is.null(nested) && length(key_columns) >= 4L) {
    sum(duplicated(nested[key_columns]))
  } else {
    NA_integer_
  }

  log_files <- list.files(file.path(run_dir, "logs"), full.names = TRUE,
                          recursive = TRUE)
  log_text <- paste(unlist(lapply(log_files, function(path) {
    tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  }), use.names = FALSE), collapse = "\n")

  run_files <- all_files[startsWith(all_files, paste0(run_dir, "/"))]
  data.frame(
    dataset = identify_dataset(run_dir),
    run_dir = relative,
    files = length(run_files),
    bytes = sum(file.info(run_files)$size, na.rm = TRUE),
    nested_exists = file.exists(nested_file),
    nested_readable = !is.null(nested),
    nested_md5 = if (file.exists(nested_file)) {
      unname(tools::md5sum(nested_file))
    } else {
      NA_character_
    },
    rows = if (is.null(nested)) NA_integer_ else nrow(nested),
    expected_rows = expected_rows,
    complete_row_grid = !is.na(expected_rows) && !is.null(nested) &&
      nrow(nested) == expected_rows,
    methods = methods,
    evaluators = evaluators,
    splits = splits,
    panel_sizes = paste(panel_sizes, collapse = ","),
    duplicate_keys = duplicate_keys,
    missing_auc = if (!is.null(nested) && "AUC" %in% names(nested)) {
      sum(!is.finite(nested$AUC))
    } else {
      NA_integer_
    },
    config_exists = file.exists(config_file),
    top_variable_genes = scalar_or_na(config, "top_variable_genes"),
    variance_filter_scope = scalar_or_na(
      config, "variance_filter_scope", "unrecorded"
    ),
    gate_method = scalar_or_na(config, "gs_gate_method", "unrecorded"),
    calibration_mode = scalar_or_na(
      config, "gs_calibration_mode", "unrecorded"
    ),
    subsample_scheme = scalar_or_na(
      config, "subsample_scheme",
      scalar_or_na(config, "gs_subsample_scheme", "unrecorded")
    ),
    n_samples = scalar_or_na(config, "n_samples"),
    n_genes = scalar_or_na(config, "n_genes"),
    string_exists = file.exists(file.path(data_dir, "string_coherence.csv")),
    stability_exists = file.exists(
      file.path(data_dir, "nogueira_native_selected.csv")
    ),
    completion_marker = grepl("Run complete:|Done\\.", log_text),
    logged_failures = lengths(regmatches(
      log_text,
      gregexpr("!! .* FAILED", log_text, perl = TRUE)
    )),
    stringsAsFactors = FALSE
  )
}

run_index <- bind_rows(lapply(run_dirs, summarise_run)) %>%
  arrange(dataset, run_dir)
write.csv(run_index, file.path(output_root, "historical_run_index.csv"),
          row.names = FALSE)

print(run_index %>%
        select(dataset, run_dir, rows, expected_rows, complete_row_grid,
               methods, evaluators, splits, top_variable_genes,
               variance_filter_scope, gate_method, calibration_mode,
               string_exists, stability_exists, completion_marker,
               logged_failures),
      row.names = FALSE)
cat("\nHistorical index written to:", output_root, "\n")
