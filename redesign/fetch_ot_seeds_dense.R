#!/usr/bin/env Rscript

# Refetch Open Targets disease associations with a much larger seed budget
# (max_seeds = 3000, min_score = 0) for the seven older benchmark datasets.
#
# Why: the frozen files (ot_seeds_<id>_n100_s0.1.rds) are capped at exactly
# 100 targets by max_seeds = 100 in get_disease_seeds_opentargets(), not by
# the 0.1 score cutoff. With only 100 targets, the var2000 candidate pools
# contain 6-30 targets and most random and real panels have zero overlap,
# which makes the Open Targets biology axis near-degenerate.
#
# The GraphQL query in package/GeneSelectR/R/biology_network.R fetches a
# single page (index 0). The API caps a page at 3000 rows, so max_seeds=3000
# fits one page; this script still paginates (copied query, index loop) so a
# smaller cap would not silently truncate. The package file is NOT modified.
#
# Output: NEW files data/r_user_cache/R/GeneSelectR/ot_seeds_<id>_n3000_s0.rds
# in the same format as the frozen n100 files (columns ensembl_id, symbol,
# score; attrs efo_id, disease, resolved_name, fetch_date). The frozen n100
# files are left untouched.

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
})
options(warn = 1)

endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"
max_seeds <- 3000L
page_size <- 3000L
cache_dir <- file.path("data", "r_user_cache", "R", "GeneSelectR")

# Disease IDs are stored as executable configuration.
source(file.path("analysis", "config.R"))
config <- get_biology_config("older7")
diseases <- stats::setNames(config$dataset, config$ontology_ids)

# Copied from get_disease_seeds_opentargets() with a page-index variable added
# for pagination; otherwise identical fields.
assoc_query <- '
  query assoc($efoId: String!, $size: Int!, $index: Int!) {
    disease(efoId: $efoId) {
      id
      name
      associatedTargets(page: {index: $index, size: $size}) {
        count
        rows {
          target { id approvedSymbol }
          score
        }
      }
    }
  }'

post_graphql <- function(query, variables, what) {
  response <- tryCatch(
    httr::POST(endpoint, body = list(query = query, variables = variables),
               encode = "json", httr::timeout(120)),
    error = function(e) NULL
  )
  if (is.null(response)) stop(what, " request failed (no response)")
  status <- httr::status_code(response)
  if (status != 200) stop(what, " request returned HTTP ", status)
  parsed <- tryCatch(
    jsonlite::fromJSON(httr::content(response, as = "text",
                                     encoding = "UTF-8")),
    error = function(e) NULL
  )
  if (is.null(parsed)) stop(what, " response contained invalid JSON")
  if (!is.null(parsed$errors)) {
    messages <- tryCatch(paste(parsed$errors$message, collapse = "; "),
                         error = function(e) "unparseable error object")
    stop(what, " GraphQL error: ", messages)
  }
  parsed
}

fetch_page <- function(efo_id, index) {
  parsed <- post_graphql(
    assoc_query,
    list(efoId = efo_id, size = page_size, index = index),
    sprintf("Open Targets association fetch (%s, page %d)", efo_id, index)
  )
  assoc <- parsed$data$disease$associatedTargets
  list(name = parsed$data$disease$name,
       count = assoc$count,
       rows = assoc$rows)
}

for (efo_id in names(diseases)) {
  dataset <- diseases[[efo_id]]
  cache_file <- file.path(
    cache_dir,
    sprintf("ot_seeds_%s_n%d_s%g.rds",
            gsub("[^A-Za-z0-9]", "", efo_id), max_seeds, 0)
  )
  if (file.exists(cache_file)) {
    cat(sprintf("[skip] %s already fetched (%s)\n", efo_id,
                basename(cache_file)))
    next
  }

  cat(sprintf("[fetch] %s (dataset %s)\n", efo_id, dataset))
  pages <- list()
  index <- 0L
  disease_name <- NA_character_
  total_count <- NA_integer_
  repeat {
    page <- fetch_page(efo_id, index)
    disease_name <- page$name
    total_count <- page$count
    rows <- page$rows
    n_rows <- if (is.null(rows) || is.null(rows$target)) 0L else nrow(rows$target)
    cat(sprintf("  page %d: %d rows (total associations: %s)\n",
                index, n_rows, format(total_count, big.mark = ",")))
    if (n_rows == 0L) break
    pages[[length(pages) + 1L]] <- data.frame(
      ensembl_id = rows$target$id,
      symbol = rows$target$approvedSymbol,
      score = rows$score,
      stringsAsFactors = FALSE
    )
    if (n_rows < page_size ||
        sum(vapply(pages, nrow, integer(1))) >= max_seeds) break
    index <- index + 1L
    Sys.sleep(0.5)  # polite spacing between pages
  }

  seeds <- if (length(pages) == 0L) {
    data.frame(ensembl_id = character(0), symbol = character(0),
               score = numeric(0), stringsAsFactors = FALSE)
  } else {
    do.call(rbind, pages)
  }
  # Same post-processing as the package function: drop NA scores, keep
  # score >= min_score (0, so nothing is dropped), sort descending, cap.
  seeds <- seeds[!is.na(seeds$score) & seeds$score >= 0, , drop = FALSE]
  seeds <- seeds[order(seeds$score, decreasing = TRUE), , drop = FALSE]
  seeds <- head(seeds, max_seeds)
  rownames(seeds) <- NULL

  # Deeper pages contain duplicated approved symbols (several Ensembl IDs
  # share one symbol, e.g. readthrough loci). The scoring code keys on
  # symbol, so keep the highest-scoring Ensembl ID per symbol.
  if (anyDuplicated(seeds$symbol)) {
    n_before <- nrow(seeds)
    seeds <- seeds[!duplicated(seeds$symbol), , drop = FALSE]
    cat(sprintf("  collapsed %d duplicate-symbol rows (kept max score)\n",
                n_before - nrow(seeds)))
  }

  attr(seeds, "efo_id") <- efo_id
  attr(seeds, "disease") <- disease_name
  attr(seeds, "resolved_name") <- disease_name
  attr(seeds, "fetch_date") <- Sys.Date()
  attr(seeds, "dataset") <- dataset
  attr(seeds, "total_associations") <- total_count
  saveRDS(seeds, cache_file)
  cat(sprintf("  saved %d seeds -> %s (score range [%.4f, %.4f])\n",
              nrow(seeds), cache_file,
              if (nrow(seeds)) min(seeds$score) else NA_real_,
              if (nrow(seeds)) max(seeds$score) else NA_real_))
  Sys.sleep(1)  # polite spacing between diseases
}

cat("dense Open Targets fetch complete\n")
