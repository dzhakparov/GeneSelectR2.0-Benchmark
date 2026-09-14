# Supplementary methods

## Predictive ranking

For each outer split, the training expression matrix was standardized from
training quantities. Candidate genes were the 2,000 most variable genes in the
training division. SOS-ALL supplied the development analysis. GSE101794,
GSE107994, GSE13355, GSE65682, GSE69683, and IMvigor210 supplied validation
analyses.

The outer design was three repeats of five stratified folds. The base seed was
42, and split construction used `42 + 1000 * repeat_index`. The reported
GeneSelectR fit used 50 five-fold resamples (`B = 50`, `subsample_scheme =
"kfold"`, `subsample_k_folds = 5`) with elastic-net mixing parameter
`alpha = 0.5` or `1.0` selected by internal out-of-bag AUC where the runner
defines the alpha grid. The utility calculation used instance-level SHAP
contribution. Recurrence and contribution were combined with an equal-weight
geometric mean.

The primary panel sizes were `k = 10, 20, 50, 100, 200, 500`. Test-fold
prediction used the deterministic three-model evaluator retained in
`redesign/R/evaluator.R`. Matched random panels were sampled from the same
training-specific 2,000-gene pool. The evaluator seed was
`420000 + 1000 * repeat_index + 100 * fold_index + panel_index`, where
`panel_index` is the position of `k` in the primary panel-size vector. The
validation runner and saved-ranking reevaluation use this seed. The full-recipe
runner invokes the evaluator default seed of 42 and uses the same default for
its matched-Random model fits. The primary matched-Random panel seed was
`99 + 1000 * repeat_index + fold_index`, with three random panels per split and
panel size.

## Outcome-permutation references

Permutation adjustment used 20 shuffled-outcome references. Each reference used
20 null resamples, giving the documented `20 permutations x 20 null fits`
setting. The expression matrix and its correlation structure were retained
while the outcome labels were permuted. The same elastic-net and five-fold
resampling settings were used for observed and null calculations.
The permutation loop used base seed 42 with deterministic offsets by
permutation index.

The adjusted recurrence and contribution values are ratios relative to the
permutation references. They are ranking quantities. They have no probability,
p-value, or significance-test interpretation. A ratio above one indicates a
higher observed score relative to the corresponding reference under the stated
calculation.

The larger-reference sensitivity uses 100 permutations for GSE107994 and
GSE13355. The first 20 permutations are checked against the saved 20-permutation
checkpoint before the 100-permutation result is accepted.

## Ablation and sensitivity analyses

The component ablation reconstructs rankings from recurrence, SHAP
contribution, mutual information, and their geometric combination while using
the saved outer splits and evaluator seeds. Calibration on/off compares the
permutation-adjusted and unadjusted ranking quantities on identical splits.
Calibration diagnostics assess null centering and leave-one-permutation-out
ratios. Random-reference sensitivity repeats matched random panels with 30
draws; the first three draws reproduce the primary baseline.

The gene-set-size analyses retain the six primary panel sizes. The adaptive
weight analysis selects a boost weight from training-only five-fold inner AUC at
`k = 20` over `w = 0, 0.5, 1, 2, 4`. The compact and greedy analyses are
IMvigor210 method-development experiments with their own settings in the
configuration registry.

## Biological assessment

Biological assessment is performed after predictive ranking. It does not enter
the reported predictive score. The four-axis older-seven analysis assesses:

- GO semantic similarity using the saved semantic references;
- Hallmark shared-membership edges;
- Open Targets disease-association scores;
- STRING v12 protein-association edges at the frozen score-400 threshold.

Each observed panel is compared with 1,000 same-size random panels from the
same training-specific candidate pool. The random-panel seed was
`940000 + 10000 * repeat_index + 100 * fold_index + panel_index`. STRING null
panels use this seed plus `1000 * n_string_mapped`. The Open Targets dense assessment uses
the seven manuscript datasets and the frozen dense seed files. STRING results
remain a separate secondary axis. GO, Hallmark, Open Targets, and STRING are
reported as separate measurements.

## DGE/GeneSelectR complementarity

The complementarity analyses use saved outer training divisions. DGE and
GeneSelectR candidate genes must match exactly before comparison. The DGE score
is recomputed as `-log10(p) * abs(t)` from the saved training division, and the
complete rank order is checked against the saved ranking. Rank percentiles are
used for cross-dataset comparison.

Stability is calculated from the exact saved top-k gene sets for `k = 10, 20,
50`. Nogueira stability, pairwise Jaccard, Sørensen-Dice, overlap coefficient,
gene recurrence, saved AUC, and AUC minus matched Random are retained as
separate fields. The optional DGE prefilter is exploratory and excluded from
the paper supplement.

For gene-group biology, the raw Open Targets, GO, and Hallmark references use
`900000 + 100000 * dataset_index + 10000 * repeat_index + 1000 * fold_index +
k_index`. The STRING reference uses
`5900000 + 100000 * dataset_index + 10000 * repeat_index + 1000 * fold_index +
k_index` and matches the number of STRING-mapped genes.
