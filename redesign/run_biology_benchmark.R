#!/usr/bin/env Rscript
# ==============================================================================
#  Biology assessment: STRING connectivity of top-50 panels (2026-08-24).
#
#  Per dataset, per arm, per split: map the top-50 panel to STRING v12
#  (score >= 400, frozen local cache), count observed edges in the
#  induced subgraph, compare against 1000 random panels of the same size
#  drawn from the SAME var2000 candidate pool. Metric = enrichment ratio
#  (observed / expected), mean over 15 splits.
#
#  Caveat from AGENTS.md stands: connectivity and prediction dissociate
#  (DGE on IMvigor210: 42x STRING, worst AUC). This is a reporting axis,
#  not a quality axis.
#
#  Every non-random arm in the completed evaluation tables is assessed. Each
#  extension persists the exact ranking it evaluated, which prevents biology
#  panels from being reconstructed with a different implementation or seed.
#
#  Usage: Rscript redesign/run_biology_benchmark.R <dataset> [budget_sec]
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
ds_arg <- if (length(args) >= 1) args[1] else stop("dataset required")
budget <- if (length(args) >= 2) as.numeric(args[2]) else 270

validation_ds <- c("GSE65682", "GSE69683", "GSE13355", "GSE107994",
                   "GSE101794")
stopifnot(ds_arg %in% c(validation_ds, "imvigor210", "sosall"))

t_start <- proc.time()[["elapsed"]]
over_budget <- function() (proc.time()[["elapsed"]] - t_start) > budget

suppressPackageStartupMessages(library(igraph))
source(file.path("redesign", "R", "run_provenance.R"))

is_validation <- ds_arg %in% validation_ds
results_root <- redesign_results_root()
split_dir <- if (is_validation)
  file.path(results_root, "validation_benchmark", ds_arg) else
  file.path(results_root, "grouped_benchmark", ds_arg)
out_dir <- if (is_validation) split_dir else
  file.path(results_root, "full_recipe", ds_arg)
require_redesign_run(split_dir)
if (!identical(out_dir, split_dir)) require_redesign_run(out_dir)
extension_drivers <- c(
  predfirst = "redesign/run_predfirst_benchmark.R",
  ensemble = "redesign/run_ensemble_eval.R",
  softprior = "redesign/run_softprior_benchmark.R",
  datadriven = "redesign/run_datadriven_benchmark.R",
  horseshoe = "redesign/run_horseshoe_benchmark.R",
  ensemble2 = "redesign/run_ensemble2_eval.R",
  wsweep = "redesign/run_wsweep.R",
  hsstab = "redesign/run_hsstab_benchmark.R",
  prune = "redesign/run_prune_benchmark.R",
  adapt = "redesign/run_adapt_benchmark.R"
)
for (extension in names(extension_drivers)) {
  require_redesign_extension(
    out_dir, extension,
    redesign_extension_sources(extension_drivers[[extension]]),
    config = list(dataset = ds_arg)
  )
}
prepare_redesign_extension(
  out_dir, "biology",
  redesign_extension_sources("redesign/run_biology_benchmark.R"),
  config = list(dataset = ds_arg),
  output_files = file.path(out_dir, "biology_string.csv")
)

nm_grouped   <- if (is_validation) "GS_full_grouped" else "full_grouped"
nm_ungrouped <- if (is_validation) "GS_full_ungrouped" else "full_ungrouped"

# ---- STRING background graph from the frozen cache (no network) -----------
info <- read.delim(gzfile("data/string_db_cache/9606.protein.info.v12.0.txt.gz"))
sym2id <- setNames(info$X.string_protein_id, info$preferred_name)
edges <- readRDS("data/string_db_cache/9606.protein.links.score400.v12.0.rds")

build_pool_graph <- function(pool_genes) {
  ids <- unique(na.omit(sym2id[pool_genes]))
  e <- edges[edges$protein1 %in% ids & edges$protein2 %in% ids, 1:2]
  g <- graph_from_data_frame(e, directed = FALSE)
  missing <- setdiff(ids, V(g)$name)
  if (length(missing) > 0) g <- add_vertices(g, length(missing),
                                             name = missing)
  g
}

panel_enrichment <- function(panel_genes, g, n_perm = 1000, seed = 7) {
  ids <- unique(na.omit(sym2id[panel_genes]))
  ids <- ids[ids %in% V(g)$name]
  if (length(ids) < 5) return(c(ratio = NA_real_, obs = length(ids) * 0,
                                exp = NA_real_, n_mapped = length(ids)))
  obs <- ecount(induced_subgraph(g, vids = ids))
  bg <- V(g)$name
  set.seed(seed)
  null_counts <- vapply(seq_len(n_perm), function(i)
    ecount(induced_subgraph(g, vids = sample(bg, length(ids)))),
    numeric(1))
  exp_ <- mean(null_counts)
  c(ratio = if (exp_ > 0) obs / exp_ else NA_real_, obs = obs, exp = exp_,
    n_mapped = length(ids))
}

# ---- rankings --------------------------------------------------------------
evaluation_paths <- sort(list.files(
  out_dir, pattern = "^eval.*[.]csv$", full.names = TRUE
))
evaluation_arms <- sort(unique(unlist(lapply(evaluation_paths, function(path) {
  tab <- read.csv(path, stringsAsFactors = FALSE)
  if (!"arm" %in% names(tab)) return(character(0))
  as.character(tab$arm)
}))))
ARMS_SAVED <- setdiff(evaluation_arms, "Random")
if (length(ARMS_SAVED) != 28L) {
  stop("Expected 28 non-random evaluated arms, found ",
       length(ARMS_SAVED), ": ", paste(ARMS_SAVED, collapse = ", "))
}

split_files <- sort(list.files(split_dir, pattern = "^split_r.*\\.rds$",
                               full.names = TRUE))
stopifnot(length(split_files) == 15)

out_path <- file.path(out_dir, "biology_string.csv")
done <- if (file.exists(out_path)) read.csv(out_path) else data.frame()
done_key <- if (nrow(done) > 0)
  paste(done$repeat_idx, done$fold_idx, done$arm) else character(0)

for (sf in split_files) {
  m <- regmatches(basename(sf),
                  regexec("split_r(\\d+)_f(\\d+)\\.rds", basename(sf)))[[1]]
  rep_idx <- as.integer(m[2]); fold_idx <- as.integer(m[3])

  need <- any(!(paste(rep_idx, fold_idx, ARMS_SAVED) %in% done_key))
  if (!need) next
  if (over_budget()) { cat("[budget] stop\n"); quit(save = "no") }

  sp <- readRDS(sf)
  pool <- sp$pools$var2000

  #  The var2000 pool differs across splits; the graph (and therefore the
  #  null) must be built from THIS split's pool. Costs ~2 s per split.
  g <- build_pool_graph(pool)

  rankings <- list()
  for (arm in ARMS_SAVED) {
    p <- file.path(out_dir, sprintf("ranking_r%d_f%d_%s.csv", rep_idx,
                                    fold_idx, arm))
    rk <- read.csv(p)$gene
    if (length(rk) == 0 && arm == nm_grouped) {
      rk <- read.csv(file.path(out_dir, sprintf("ranking_r%d_f%d_%s.csv",
                       rep_idx, fold_idx, nm_ungrouped)))$gene
    }
    rankings[[arm]] <- rk
  }

  rows <- list()
  for (arm in ARMS_SAVED) {
    key <- paste(rep_idx, fold_idx, arm)
    if (key %in% done_key) next
    panel <- head(rankings[[arm]], 50)
    enr <- panel_enrichment(panel, g)
    rows[[length(rows) + 1]] <- data.frame(
      repeat_idx = rep_idx, fold_idx = fold_idx, arm = arm,
      ratio = enr["ratio"], obs_edges = enr["obs"],
      exp_edges = enr["exp"], n_panel = length(panel),
      n_mapped = enr["n_mapped"], pool_size = length(pool))
    done_key <- c(done_key, key)
  }
  write.table(do.call(rbind, rows), out_path,
              append = file.exists(out_path), sep = ",",
              row.names = FALSE, col.names = !file.exists(out_path))
  cat(sprintf("[%s r%d f%d] biology done\n", ds_arg, rep_idx, fold_idx))
}

cat(sprintf("\nDone dataset=%s (%.0f s)\n", ds_arg,
            proc.time()[["elapsed"]] - t_start))
