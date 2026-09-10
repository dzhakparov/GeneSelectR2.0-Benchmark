# GeneSelectR and differential-expression complementarity analyses

## Purpose

These scripts test whether differential-expression analysis (DGE) and
GeneSelectR provide related but different gene rankings. They use the saved
rankings and the saved train/test divisions from the existing benchmark. They
do not change the GeneSelectR score or refit any feature-selection method.

SOS-ALL, the six validation datasets and the four additional datasets retain
separate dataset-role labels in every primary output. The additional datasets
are described as follow-up evaluations because earlier approximate analyses
were performed before the final saved runs.

Run every command from the repository root.

For an isolated small check, set `GENESELECTR_COMPLEMENTARITY_RESULTS_DIR` to
a temporary directory. The default output remains `results/` under this
analysis directory.

## Files

### `00_utils.R`

Shared dataset locations, file readers, checks, gene-set similarity functions,
the repository's existing Nogueira calculation, predictive-result loading and
biological-resource loading.

This file maps the saved benchmark names to the following names used in the new
outputs:

- `GeneSelectR`
- `DGE`
- `Random forest`
- `Boruta`
- `mRMR`
- `LASSO`
- `Elastic net`

### `01_gene_set_stability.R` — light

Uses all saved top-k rankings for `k = 10, 20, 50`. For every dataset and
method it calculates:

- Nogueira stability;
- mean and median pairwise Jaccard similarity;
- mean Sørensen-Dice similarity;
- mean overlap coefficient;
- per-gene selection counts across the 15 train/test divisions;
- mean saved AUC and mean saved AUC minus same-size random genes.

Outputs:

- `results/gene_set_stability_summary.csv`
- `results/gene_cross_split_recurrence.csv`
- `results/top_gene_sets_by_split.csv`
- `results/candidate_gene_set_diagnostics.csv`

The training-specific candidate genes vary between divisions. The primary
Nogueira value therefore uses only genes that were available in all 15
divisions and retains the actual selected genes from each division. A second
fixed-gene calculation reranks within the common genes and selects exactly k
genes. The union-based value is retained as a labelled diagnostic because it
counts unavailable genes as unselected. These values are stored as
`nogueira_stability`, `nogueira_fixed_common_genes_reranked_topk`, and
`nogueira_exact_topk_union`. Jaccard and Dice use the exact saved top-k gene
sets and are the direct measures of overlap between the reported gene sets.

### `02_dge_geneselectr_concordance.R` — moderate

Compares DGE and GeneSelectR within each saved outer training division. The
candidate genes must match exactly before any comparison is calculated.

Outputs:

- `results/dge_geneselectr_gene_scores.csv`
- `results/dge_geneselectr_split_correlations.csv`
- `results/dge_geneselectr_correlation_summary.csv`
- `results/dge_geneselectr_topk_overlap.csv`
- `results/dge_geneselectr_overlap_summary.csv`

The saved DGE CSV files contain the ordered genes only. The script recomputes
the exact benchmark calculation, `-log10(p) * abs(t)`, from each saved training
division. The complete recomputed order must match the saved ranking before
the result is accepted. The output includes the ranking statistic, t statistic,
p-value, Benjamini-Hochberg adjusted p-value, and the existing adjusted-p-value
selection indicator. The existing DGE implementation did not calculate an
effect size, so `dge_log_fold_change` remains missing. Rank percentiles remain
the primary cross-dataset comparison because the numerical DGE scale depends
on the dataset.

### `03_gene_group_characterization.R` — light

Uses the output from script 02. At `k = 20` and `k = 50`, genes are classified
within each train/test division as:

- `SHARED`: top-k under both methods;
- `DGE_ONLY`: DGE top-k only;
- `GS_ONLY`: GeneSelectR top-k only.

Outputs:

- `results/gene_group_membership.csv`
- `results/gene_group_summary.csv`
- `results/gene_group_recurrence_long.csv`

These are ranking groups. The names do not indicate statistical significance,
causality or clinical validation.

### `04_gene_group_biology.R` — moderate

Uses the groups from script 03. Each group is compared with 1,000 same-size
random gene sets from the same training-specific 2,000 candidate genes. It
calculates:

- Open Targets disease-association score sum;
- mean saved GO semantic-similarity score;
- number of gene pairs sharing a Hallmark gene set;
- number of recorded STRING v12 protein associations.

STRING random sets are also matched on the number of genes mapped to STRING,
as in the existing benchmark. Empty groups and groups with fewer than two
mapped genes receive a missing value and an explicit reason.

Outputs:

- `results/gene_group_biology_by_split.csv`
- `results/gene_group_biology_summary.csv`

The script uses the association files recorded in `analysis/config.R`. These
are the common frozen Open Targets inputs available for all 11 datasets. They
contain up to 100
associations at the existing 0.1 retrieval threshold. The denser 3,000-record
files are available for the seven manuscript benchmark datasets only. Mixing them into
the primary 11-dataset comparison would change database coverage by dataset
group, so the script does not substitute them automatically.

### `05_exploratory_dge_prefilter.R` — heavy

For each saved outer training division, this script takes the top 100, 250 or
500 genes from the saved DGE ranking, orders those genes using the saved
GeneSelectR ranking, selects `k = 10, 20, 50`, and fits the unchanged predictive
evaluator on the resulting gene set. The DGE and GeneSelectR feature-selection
methods are not refitted.

Each train/test division is saved as a separate checkpoint. Existing saved DGE,
GeneSelectR and random-gene-set AUC values are joined for comparison.

Outputs after assembly:

- `results/exploratory_dge_prefilter_results.csv`
- `results/exploratory_dge_prefilter_gene_sets.csv`
- `results/exploratory_dge_prefilter_stability.csv`

The stability table contains the exact selected gene sets and a sensitivity
calculation restricted to genes available in every train/test division.

All output rows are labelled `exploratory_dge_prefilter`. The three values of
M and all three values of k must be reported. They must not be selected using
the test AUC.

### `06_make_figures.R` — light after scripts 01–04

Reads only the CSV files created above. It writes PDF and 300-dpi PNG figures,
plus the exact CSV used by each figure.

Outputs include:

- `results/figures/Figure_2_prediction_and_stability.*`
- `results/figures/Figure_3_DGE_GeneSelectR_complementarity.*`
- `results/figures/Figure_4_gene_group_biology.*`
- stability and component-comparison supplementary figures;
- `results/figure_source_data/*.csv`.

## Existing benchmark inputs

The scripts use:

- development rankings and train/test divisions under
  `redesign/results_corrected/validation_benchmark/`,
  `redesign/results_corrected/full_recipe/`, and
  `redesign/results_corrected/grouped_benchmark/`;
- additional-dataset rankings and train/test divisions under
  `redesign/results_frozen_external_exact_2026-08-31/validation_benchmark/`;
- `eval_deterministic.csv` for development prediction results;
- `eval_results.csv` and `competitor_eval_results.csv` for additional-dataset
  prediction results;
- `package/GeneSelectR/R/utils.R::compute_nogueira_stability()`;
- saved `older7_biology_semantic_reference.rds` and
  `biology_semantic_reference_v2.rds` files;
- frozen Open Targets files listed in the two biology configuration tables;
- MSigDB Hallmark genes through the existing `hallmark_sets_all()` function;
- frozen STRING v12 mapping and score-400 association files;
- `redesign/R/evaluator.R` and `redesign/R/bio_prior.R` for the optional
  predictive evaluation.

Every expected dataset has 15 saved train/test divisions and complete rankings
for all seven methods. Every saved ranking contains exactly 2,000 unique genes
and is checked against the training-specific candidate genes.

## Commands

Primary saved-result analyses:

```bash
Rscript redesign/complementarity_analysis/01_gene_set_stability.R
Rscript redesign/complementarity_analysis/02_dge_geneselectr_concordance.R
Rscript redesign/complementarity_analysis/03_gene_group_characterization.R
Rscript redesign/complementarity_analysis/04_gene_group_biology.R
Rscript redesign/complementarity_analysis/06_make_figures.R
```

Scripts 01 and 02 can be limited to named datasets for a small check, for
example:

```bash
Rscript redesign/complementarity_analysis/01_gene_set_stability.R GSE101794
Rscript redesign/complementarity_analysis/02_dge_geneselectr_concordance.R GSE101794
```

The output files from a limited run contain only the requested datasets and
will replace prior files in this new analysis directory. Run the full command
after any limited check.

Optional exploratory predictive evaluation:

```bash
Rscript redesign/complementarity_analysis/05_exploratory_dge_prefilter.R run GSE101794
Rscript redesign/complementarity_analysis/05_exploratory_dge_prefilter.R run GSE107994
Rscript redesign/complementarity_analysis/05_exploratory_dge_prefilter.R run GSE13355
Rscript redesign/complementarity_analysis/05_exploratory_dge_prefilter.R run GSE65682
Rscript redesign/complementarity_analysis/05_exploratory_dge_prefilter.R run GSE69683
Rscript redesign/complementarity_analysis/05_exploratory_dge_prefilter.R run imvigor210
Rscript redesign/complementarity_analysis/05_exploratory_dge_prefilter.R run sosall
Rscript redesign/complementarity_analysis/05_exploratory_dge_prefilter.R run GSE16879
Rscript redesign/complementarity_analysis/05_exploratory_dge_prefilter.R run GSE91061
Rscript redesign/complementarity_analysis/05_exploratory_dge_prefilter.R run GSE92415
Rscript redesign/complementarity_analysis/05_exploratory_dge_prefilter.R run GSE206285
Rscript redesign/complementarity_analysis/05_exploratory_dge_prefilter.R assemble
```

The optional step fits 1,485 prediction ensembles across 165 saved train/test
divisions. Run one dataset at a time. The split-level checkpoints allow a
stopped command to resume without repeating completed divisions.

## Methodological decisions required before interpretation

1. Candidate genes vary between training divisions. Use `nogueira_stability`,
   which has a fixed common gene universe, as the primary Nogueira value.
   Report exact-set Jaccard and Dice beside it. Treat the union-based Nogueira
   value as a diagnostic only.
2. AUC and stability describe different properties. There is no conventional
   p-value across the 15 repeated cross-validation results because the same
   participants occur in more than one test division.
3. DGE statistics were not saved and are recomputed from the saved training
   data using the exact benchmark formula. The existing benchmark did not
   calculate an effect size, so no fold change is reported.
4. Biological measurements partly reuse databases that were already examined
   during method development. They characterize the ranking groups. They do
   not establish that one feature-selection method is generally more
   biologically relevant.
5. The exploratory DGE-to-GeneSelectR workflow uses fixed M and k values and must
   report every result. A new default method would require a separate decision
   made before testing on new datasets.
