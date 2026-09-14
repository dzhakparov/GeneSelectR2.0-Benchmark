# Synthesized analysis report

**Scope:** GeneSelectR paper supplementary analyses  
**Report date:** 2026-09-14  
**Primary evidence:** deterministic benchmark and extension reports in the
external result archive; committed figures and reviewer-facing documentation
in this worktree.

## Summary

Seven datasets were evaluated: SOS-ALL (`sosall`) as the development cohort,
and GSE101794, GSE107994, GSE13355, GSE65682, GSE69683, and IMvigor210 as
validation cohorts. The deterministic comparison included 27 unique methods,
15 outer splits per dataset, six panel sizes, and 2,520 evaluation rows per
dataset.

GeneSelectR-derived configurations ranked first in five of seven datasets.
mRMR ranked first in two datasets. The best cross-dataset result among the two
base GeneSelectR configurations was obtained by `GS_full_ungrouped`, with mean
dataset rank 10.29 and mean AUC minus matched Random of +0.063 over
`k = 10, 20, 50`. Classical selectors remained competitive, with higher
primary deltas in sepsis and SOS-ALL.

The results support GeneSelectR as a reproducible panel-ranking method with
positive matched-Random performance across the evaluated cohorts. The data
support descriptive comparative conclusions. Repeat-cluster intervals for the
balanced comparison of `GS_full_ungrouped` with RF importance, DGE, Boruta,
mRMR, LASSO, and ElasticNet all included zero under the conservative
repeat-level analysis.

## Cohorts and scope

| Role | Dataset | Biological context | Scope |
|---|---|---|---|
| Development | SOS-ALL (`sosall`) | Atopic eczema | Primary |
| Validation | GSE101794 | Inflammatory bowel disease / Crohn's disease | Primary |
| Validation | GSE107994 | Tuberculosis | Primary |
| Validation | GSE13355 | Psoriasis | Primary |
| Validation | GSE65682 | Sepsis | Primary |
| Validation | GSE69683 | Asthma | Primary |
| Validation | IMvigor210 (`imvigor210`) | Urinary bladder carcinoma | Primary |
| Follow-up | GSE16879, GSE91061, GSE92415, GSE206285 | Additional cohorts | Exploratory; excluded from paper supplement |

Targeted-assay analyses, p009 analyses, and unrelated multi-omics analyses
are outside this report. The four follow-up datasets remain indexed in
[`exploratory_not_in_paper/README.md`](exploratory_not_in_paper/README.md).
The restored figure files preserve the earlier 11-dataset figure layout and
therefore include these follow-up cohorts where applicable. Primary numerical
conclusions remain restricted to the seven-dataset benchmark.

## Analysis design

Training expression values were standardized from training quantities. The
candidate pool contained the 2,000 most variable training genes. The outer
design used three repeats of five stratified folds, with base seed 42.
GeneSelectR used 50 five-fold subsample fits, instance-level SHAP contribution,
recurrence, permutation-adjusted utility, and an equal-weight geometric mean
for the current ranking. Elastic-net mixing parameters were selected from the
runner-defined alpha grid by internal out-of-bag AUC.

Panels contained 10, 20, 50, 100, 200, or 500 genes. Test-fold prediction
used the deterministic three-model evaluator retained in
`redesign/R/evaluator.R`. Matched Random panels were sampled from the same
split-specific 2,000-gene pool. Primary comparisons use mean AUC minus the
matched Random reference.

The permutation-adjustment reference used 20 outcome permutations and 20 null
subsample fits per permutation. The ratios are ranking quantities. They do not
have p-value or significance-test interpretations.

## Primary deterministic benchmark

### Dataset-level results

| Dataset | Highest-ranked method | Mean AUC minus matched Random |
|---|---|---:|
| GSE101794, Crohn's disease | mRMR | +0.045 |
| GSE107994, tuberculosis | `soft_gs` | +0.062 |
| GSE13355, psoriasis | `cbgs_prune` | +0.004 |
| GSE65682, sepsis | mRMR | +0.081 |
| GSE69683, asthma | `GS_full_grouped` | +0.110 |
| IMvigor210 | `cb_gs_w2` | +0.134 |
| SOS-ALL | `cbgs_prune` | +0.128 |

`GS_full_ungrouped` had mean primary delta +0.063 and mean dataset rank 10.29.
Relative to the strongest classical selector, its primary delta was higher by
0.027 in asthma, 0.003 in IMvigor210, and 0.003 in tuberculosis. Differences
were below 0.001 in Crohn's disease and psoriasis. Its primary delta was lower
by 0.018 in sepsis and 0.030 in SOS-ALL.

The deterministic evaluator correction changed primary deltas by a mean
absolute value of 0.0035 and dataset ranks by a mean absolute value of 1.42.
The largest method-level primary-delta change was 0.0157. The dataset-level
interpretation was retained after correction, with changes in exact method
order.

Cross-validation selected biological weight zero in 99 of 105 outer splits.
The mean method-level Spearman correlation between STRING enrichment and
predictive delta was -0.010. The corresponding disease-association
correlation was +0.034.

### Primary configuration

`GS_full_ungrouped` was selected as the primary base configuration from the
deterministic development comparison. Consensus panels of 10, 20, 50, and 100
genes were locked from the 15 GSE107994 outer-split rankings using mean
reciprocal rank. The lock was created after access to GSE19442. Evaluation of
the newly locked configuration requires a future independent cohort.

## Method-development comparisons

The 27-method comparison included the current GeneSelectR configurations,
classical selectors, prediction-first rankings, soft-prior and data-driven
module variants, rank ensembles, weight variants, horseshoe rankings, and
redundancy-pruned rankings.

The best method varied by dataset. The current evidence supports a
configuration-level comparison across cohorts. It does not support a single
universally dominant selector. The leading GeneSelectR configuration had
mean primary delta +0.063; the leading classical comparator had mean rank
11.29 in the cross-dataset summary.

Additional method-development analyses were retained for provenance:

| Analysis family | Short result | Interpretation |
|---|---|---|
| Prediction-first ranking | Included as `predfirst` and `predfirst_raw` methods | Descriptive method comparison |
| Soft-prior and data-driven modules | Included in `soft_*`, `dc_*`, and `cb*` methods | Dataset-dependent performance |
| Rank ensembles | `ens2`, `ens3`, `ens_rank`, and `ens_rank3` evaluated | Descriptive combinations |
| Weight sweep | `cb_gs_w2` was the leading method in IMvigor210 | Weight effects varied by dataset |
| Adaptive weight | Training-only inner AUC selection evaluated over the configured weight grid | No uniform cross-dataset advantage was established |
| Horseshoe and horseshoe stability | Included as `hs` and `hs_stab` methods | Exploratory Bayesian and stability comparisons |
| Redundancy pruning | Included as `cbgs_prune` and `dcpf_prune` methods | Pruned methods ranked first in psoriasis and SOS-ALL |
| Compact and greedy IMvigor210 panels | Fold-level method-development experiments | Exploratory; outside the primary configuration lock |
| T-Rex diagnostics | IMvigor210 fold-1 diagnostic probes | Exploratory; no primary conclusion |

The adaptive nested table contained 630 rows across six datasets. For
`GS_full_ungrouped`, the mean selected panel size was 20.7, 26.0, 32.0, 22.7,
44.0, and 96.7 genes for GSE101794, GSE107994, GSE13355, GSE69683, IMvigor210,
and SOS-ALL, respectively. The corresponding mean AUC-minus-Random values
were +0.051, +0.074, +0.002, +0.080, +0.126, and +0.065.

## Component ablation

The saved outer-training fits were used to reconstruct five rankings:
recurrence-only, SHAP-only, MI-only, SHAP × MI, and the current combined
recipe. Each ranking was evaluated at `k = 10, 20, 50, 100, 200, 500`.

For the primary short-result view at `k = 10, 20, 50`:

- The current recipe was the top variant in 4 of 21 dataset-by-panel-size
  cells.
- The current recipe was within 0.01 AUC of the best variant in 17 of 21
  cells.
- MI was strongest on selected SOS-ALL and IMvigor210 cells.
- SHAP-only performance was lower than the current recipe by 0.056 at
  IMvigor210 `k = 10` and 0.065 at SOS-ALL `k = 10`.
- Recurrence-only was competitive across all seven datasets.
- 79 of 84 dataset-by-panel-size-variant mean contrasts had absolute
  differences below 0.03.

GSE13355 had matched-Random AUC values of 0.994–0.998. Its ablation deltas
were within approximately ±0.005. At `k = 20` and `k = 50`, recurrence-only
and current rankings produced identical panels on all 15 splits.

The conservative repeat-cluster inference was applied to the dataset-balanced
variant contrasts. All 12 reported ablation intervals included zero under the
repeat-level t analysis.

## Calibration analyses

### Calibration on/off

Calibration-on and calibration-off rankings were reconstructed from identical
saved fits and evaluator seeds. The dataset-balanced primary endpoint over
`k = 10, 20, 50` gave the following raw-minus-calibrated differences:

| Comparison | Raw-minus-calibrated AUC | Short result |
|---|---:|---|
| Combined ranking | -0.0062 | Calibrated ranking had higher mean AUC |
| Recurrence-only ranking | -0.0044 | Difference was within 0.005 AUC |
| SHAP × MI ranking | -0.0067 | Calibrated ranking had higher mean AUC |

The current ranking reproduced all 630 saved dataset/split/panel AUC cells
with maximum absolute difference 5.55e-16.

### Null-reference diagnostics

The 20-permutation leave-one-permutation-out assessment showed:

- Utility ratios were well centered on five datasets, with median ratios
  0.95–1.00 and 90–99% of ratios within [0.5, 2].
- GSE107994 had median utility ratio 0.0625, standard deviation of log2 ratio
  2.08, and 43% of ratios within [0.5, 2] in the earlier reference.
- GSE13355 had median utility ratio 1.00, standard deviation of log2 ratio
  2.24, and 51% of ratios within [0.5, 2].
- Stability-ratio medians ranged from 0.39 to 1.00, with 29–68% of ratios
  within [0.5, 2].
- Null-model failure fraction was zero in all seven datasets.

A reproducibility error in the calibration-diagnostic null rebuild was fixed.
The null rebuild lacked the ambient RNG seed used by the saved benchmark fits.
Wrapping the rebuild in `withr::with_seed(42L, ...)` reproduced the saved
stability and utility ratios exactly. The extraction stage used saved fits and
was unaffected. All seven post-fix QA assemblies passed.

### 100-permutation follow-up

The larger-reference confirmation was completed for the two exceptional
datasets. The first 20 permutations matched the saved checkpoints exactly;
each run contained 2,000 null fits with zero failures.

| Dataset | Pillar | Median LOPO ratio | SD of log2 ratio | Fraction within [0.5, 2] |
|---|---|---:|---:|---:|
| GSE107994 | Stability | 0.643 | 1.002 | 0.590 |
| GSE107994 | Utility | 0.0625 | 2.084 | 0.161 |
| GSE13355 | Stability | 0.832 | 1.050 | 0.679 |
| GSE13355 | Utility | 0.0625 | 2.204 | 0.271 |

The larger reference reduced stability-ratio dispersion. Utility-ratio
instability persisted. The diagnostic source attributes this behavior to
sparse observed utility combined with a small epsilon relative to the null
mean for floor genes. A revised epsilon rule was not implemented in the
reported workflow.

## Random-reference sensitivity

Thirty matched Random panels were sampled per split, panel size, and dataset.
The first three draws reproduced the primary baseline with maximum absolute
difference 5.6e-16.

Mean Random AUC increased with panel size on every dataset. The range across
datasets was 0.542–0.994 at `k = 10`, 0.550–0.997 at `k = 20`, and
0.574–0.998 at `k = 50`. Within-split Monte Carlo standard errors for the
30-draw mean ranged from 0.0004 to 0.026. Across-split standard deviation
ranged from 0.004 to 0.056 and was the larger source of variation.

GSE13355 had the highest Random AUC because the expression array contains
strong shared signal. Matched AUC-minus-Random remains the appropriate
descriptive scale for the benchmark.

## Gene-set size and stability

The six fixed panel sizes were evaluated in the primary and extension analyses.
The mean GeneSelectR Nogueira stability across the seven datasets was:

| Panel size | Mean stability | Median stability | Range across dataset summaries |
|---:|---:|---:|---:|
| 10 | 0.439 | 0.431 | 0.280–0.665 |
| 20 | 0.454 | 0.411 | 0.304–0.652 |
| 50 | 0.458 | 0.505 | 0.331–0.645 |

Stability increased modestly with panel size in the balanced means. The
top-k recurrence and utility rankings remained correlated because both were
derived from the same subsample fits.

## Biological assessment

Biological assessment was performed after predictive ranking. Biological
measurements were retained as separate axes and excluded from the predictive
score. The older-seven analysis evaluated GO semantic similarity, Hallmark
shared-membership edges, Open Targets disease-association scores, and STRING
protein-association edges against 1,000 matched Random panels.

The gene-group summary contained 168 dataset/group/panel/metric cells. Across
the 42 dataset/group/panel summaries for each axis, the median observed-to-
Random ratios were 1.38 for GO semantic similarity, 1.42 for Hallmark edges,
0.00 for Open Targets disease-association sum, and 4.09 for STRING edges.
These values are axis-specific descriptive summaries. They do not define a
combined biological score.

The dense Open Targets follow-up refetched associations at a newer platform
snapshot, data version 26.06, with up to 3,000 associations. At the primary
score cutoff of 0.05, zero-overlap fractions at `k <= 50` ranged from 0.3% to
51.7% across datasets. At the continuity cutoff of 0.1, the range was 1.3% to
74.6%. SOS-ALL retained sparse coverage. The dense run is a secondary coverage
assessment because its platform snapshot differs from the frozen 2026-08
resources.

## DGE and GeneSelectR complementarity

DGE and GeneSelectR rankings were reconstructed from the same saved outer
training divisions. The DGE rank order was independently checked after
recomputation from `-log10(p) * abs(t)`.

Median Spearman correlation between DGE and GeneSelectR rank across the seven
datasets ranged from 0.116 to 0.620. The mean top-k Jaccard overlap between the
two rankings was 0.292 at `k = 10`, 0.290 at `k = 20`, and 0.258 at `k = 50`.
The two procedures therefore selected overlapping gene sets with substantial
rank differences.

GeneSelectR Nogueira stability had means of 0.439, 0.454, and 0.458 at
`k = 10, 20, 50`, respectively. The DGE comparison retained rank correlation,
top-k overlap, pairwise stability, recurrence, AUC, and AUC-minus-Random as
separate measurements. A DGE prefilter is classified as exploratory and is
excluded from the paper supplement.

## External holdout

GSE19442 was evaluated as a secondary locked external holdout. It contained
51 independent pretreatment participants: 20 with active tuberculosis and 31
with latent tuberculosis. The prespecified manifest contained 17 earlier
GeneSelectR variants and classical comparators at panel sizes 10, 20, 50, and
100.

The run produced 4,080 fold-level metrics and 41,616 sample predictions with
complete AUC output. A second complete run was byte-identical. Across panel
sizes 10, 20, and 50, ensemble mean AUC was 0.992–0.995 for every named
selector, while the locked Random mean was 0.846.

The cohort supports portability of the earlier locked panels. Its near-ceiling
performance provides limited separation between GeneSelectR and classical
selectors. The latest redesign configurations were absent from the manifest,
so this holdout does not evaluate the current `GS_full_ungrouped` lock.

## Reproducibility and quality assurance

The following checks were completed in the archived analyses:

| Check | Result |
|---|---|
| Saved-ranking AUC reconstruction | Maximum absolute difference 4.4e-16 |
| Component-ablation AUC reconstruction | Maximum absolute difference 5.6e-16; all QA passed |
| Random-reference first-three-draw reproduction | Maximum absolute difference 5.6e-16; all QA passed |
| Calibration ratio reconstruction after seed fix | Maximum difference 0; all seven QA assemblies passed |
| `glmnet` selector replay | 20,400 fits bit-identical at the saved coefficient matrices |
| `glmnet` path warnings | 37 selector fits of 20,400, 0.18%; no boundary selector fits |
| Deep-path evaluator refit | Three warning-affected rows; maximum AUC change 3.331e-16 |
| Repository-content check | Passed |

The clean worktree contains the scripts, committed rendered figures, figure
manifest, numerical summary, methods, provenance, and this synthesized report.
Expression matrices, metadata, annotation caches, fitted objects, and source
result tables remain under ignored paths in the external analysis archive. A
fresh full numerical rerun was therefore not performed from this checkout.

The package interface used by the maintained runners was verified at revision
`621be0c1`. The simplified package revision `db1f6528` has a different
`geneselectr2_fit` interface and is incompatible with the archived reported
workflow.

## Interpretation and reporting recommendations

1. Report the primary endpoint as mean AUC minus matched Random, with dataset,
   panel size, and method counts stated explicitly.
2. State that GeneSelectR-derived configurations ranked first in five of seven
   primary datasets and that classical selectors remained competitive.
3. Use `GS_full_ungrouped` as the current primary base configuration, with
   external evaluation of its newly locked panels pending.
4. Describe component and biological analyses as descriptive sensitivity and
   interpretation analyses. Biological measurements are separate from the
   predictive score.
5. Include the GSE13355 ceiling behavior and the GSE107994 utility-calibration
   instability as dataset-specific limitations.
6. Preserve the 100-permutation follow-up as a completed diagnostic for the
   two exceptional datasets. Treat epsilon-rule revision as future method
   development.

## Linked supplementary artifacts

- [Numerical results summary](results_summary.md)
- [Supplementary methods](methods/README.md)
- [Configuration registry](configurations/registry.md)
- [Figure index and captions](figures/README.md)
- [Figure manifest and QA](figures/manifest.md)
- [Provenance and reproducibility checklist](provenance/README.md)
