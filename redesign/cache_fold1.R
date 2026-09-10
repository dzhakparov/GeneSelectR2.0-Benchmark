#!/usr/bin/env Rscript
#  One-off: build the repeat-1/fold-1 train matrix and cache it as RDS so the
#  T-Rex diagnostic runs don't each pay the TMM preprocessing cost.
suppressPackageStartupMessages(library(edgeR))
source(file.path("redesign", "R", "imvigor210_data.R"))

dat <- load_imvigor210()
folds <- make_stratified_folds(dat$outcome, k_folds = 5, seed = imv_random_seed)
test_idx <- folds[[1]]
train_idx <- setdiff(seq_along(dat$outcome), test_idx)
pp <- preprocess_split(dat$raw_counts, train_idx, test_idx, top_genes = 2000)
y <- as.integer(dat$outcome[train_idx]) - 1L

out <- file.path("redesign", "results", "trex_sanity", "fold1_train.rds")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(X = pp$train, y = y, train_idx = train_idx,
             test_idx = test_idx), out)
cat(sprintf("saved %s | n=%d p=%d responders=%d\n",
            out, nrow(pp$train), ncol(pp$train), sum(y)))
