# ==============================================================================
#  Provenance controls for redesign benchmark checkpoints.
#
#  A checkpoint is valid only for the source files, input files, and run
#  configuration recorded in its manifest. Corrected runs use a separate
#  result root by default so the earlier benchmark files remain available for
#  comparison and cannot be resumed by the corrected implementation.
# ==============================================================================

REDESIGN_CACHE_VERSION <- "2026-08-26-corrected-v1"

redesign_results_root <- function() {
  Sys.getenv("GENESELECTR_REDESIGN_RESULTS_ROOT",
             file.path("redesign", "results_corrected"))
}

redesign_file_hashes <- function(paths) {
  paths <- sort(unique(paths[file.exists(paths)]))
  if (length(paths) == 0) return(character(0))
  hashes <- unname(tools::md5sum(paths))
  names(hashes) <- paths
  hashes
}

redesign_manifest <- function(config, source_files, input_files) {
  packages <- c("Boruta", "easierData", "edgeR", "ExperimentHub", "glmnet",
                "igraph", "mRMRe", "ranger", "SummarizedExperiment",
                "withr", "xgboost")
  versions <- vapply(packages, function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
      as.character(utils::packageVersion(package))
    } else {
      NA_character_
    }
  }, character(1))
  list(
    schema_version = 1L,
    cache_version = REDESIGN_CACHE_VERSION,
    configuration = config,
    runtime = list(R = as.character(getRversion()), packages = versions),
    source_hashes = redesign_file_hashes(source_files),
    input_hashes = redesign_file_hashes(input_files)
  )
}

redesign_object_hash <- function(object) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, version = 3, compress = FALSE)
  unname(tools::md5sum(path))
}

prepare_redesign_run <- function(out_dir, config, source_files,
                                 input_files = character(0)) {
  expected <- redesign_manifest(config, source_files, input_files)
  manifest_path <- file.path(out_dir, "run_manifest.rds")

  if (file.exists(manifest_path)) {
    observed <- readRDS(manifest_path)
    if (!identical(observed, expected)) {
      stop(paste0(
        "Checkpoint provenance does not match the current code, inputs, or ",
        "configuration: ", out_dir, ". Use a new ",
        "GENESELECTR_REDESIGN_RESULTS_ROOT."), call. = FALSE)
    }
    return(invisible(observed))
  }

  existing <- if (dir.exists(out_dir)) list.files(out_dir, all.files = TRUE,
                                                  no.. = TRUE) else character(0)
  if (length(existing) > 0) {
    stop(paste0(
      "Checkpoint directory has no provenance manifest: ", out_dir,
      ". Use a new GENESELECTR_REDESIGN_RESULTS_ROOT."), call. = FALSE)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(expected, manifest_path)
  writeLines(REDESIGN_CACHE_VERSION,
             file.path(out_dir, "run_manifest.version"))
  invisible(expected)
}

require_redesign_run <- function(out_dir) {
  manifest_path <- file.path(out_dir, "run_manifest.rds")
  if (!file.exists(manifest_path)) {
    stop("Missing corrected-run provenance manifest: ", out_dir,
         call. = FALSE)
  }
  manifest <- readRDS(manifest_path)
  if (!identical(manifest$cache_version, REDESIGN_CACHE_VERSION)) {
    stop("Run provenance version is incompatible: ", out_dir,
         call. = FALSE)
  }
  invisible(manifest)
}

require_redesign_extension <- function(out_dir, extension, source_files,
                                       config = list()) {
  require_redesign_run(out_dir)
  manifest_path <- file.path(
    out_dir, sprintf("extension_manifest_%s.rds", extension))
  version_path <- file.path(
    out_dir, sprintf("extension_manifest_%s.version", extension))
  if (!file.exists(manifest_path) || !file.exists(version_path)) {
    stop("Missing corrected extension provenance: ", extension, " in ",
         out_dir, call. = FALSE)
  }
  if (!identical(readLines(version_path, warn = FALSE),
                 REDESIGN_CACHE_VERSION)) {
    stop("Extension provenance version is incompatible: ", extension,
         call. = FALSE)
  }
  expected <- redesign_manifest(
    c(list(extension = extension), config), source_files,
    file.path(out_dir, "run_manifest.rds"))
  observed <- readRDS(manifest_path)
  if (!identical(observed, expected)) {
    stop("Extension provenance does not match current code: ", extension,
         " in ", out_dir, call. = FALSE)
  }
  invisible(observed)
}

prepare_redesign_extension <- function(
    out_dir, extension, source_files, config = list(),
    output_files = file.path(out_dir, sprintf("eval_%s.csv", extension))) {
  require_redesign_run(out_dir)
  expected <- redesign_manifest(
    c(list(extension = extension), config), source_files,
    file.path(out_dir, "run_manifest.rds"))
  manifest_path <- file.path(
    out_dir, sprintf("extension_manifest_%s.rds", extension))
  version_path <- file.path(
    out_dir, sprintf("extension_manifest_%s.version", extension))

  if (file.exists(manifest_path)) {
    observed <- readRDS(manifest_path)
    if (!identical(observed, expected)) {
      stop(paste0(
        "Extension checkpoint provenance does not match the current code: ",
        extension, " in ", out_dir, ". Use a new ",
        "GENESELECTR_REDESIGN_RESULTS_ROOT."), call. = FALSE)
    }
    return(invisible(observed))
  }

  if (any(file.exists(output_files))) {
    stop(paste0(
      "Extension output has no provenance manifest: ", extension, " in ",
      out_dir, ". Use a new GENESELECTR_REDESIGN_RESULTS_ROOT."),
      call. = FALSE)
  }

  saveRDS(expected, manifest_path)
  writeLines(REDESIGN_CACHE_VERSION, version_path)
  invisible(expected)
}

redesign_extension_sources <- function(driver) {
  sort(unique(c(
    driver,
    list.files(file.path("redesign", "R"), pattern = "[.]R$",
               full.names = TRUE),
    list.files(file.path("package", "GeneSelectR", "R"), pattern = "[.]R$",
               full.names = TRUE)
  )))
}

redesign_worker_count <- function(default = 2L) {
  value <- suppressWarnings(as.integer(Sys.getenv("GENESELECTR_N_CORES",
                                                   as.character(default))))
  if (is.na(value) || value < 1L) {
    stop("GENESELECTR_N_CORES must be a positive integer.", call. = FALSE)
  }
  value
}
