# GeneSelectR configuration registry

The registry records maintained scripts and their status. Parameters are taken
from the script headers and configuration blocks in this archive. Output data
are ignored by Git.

## Reported workflow

| Analysis | Maintained entry point | Configuration summary | Status |
|---|---|---|---|
| Primary seven-dataset workflow | `redesign/run_corrected_benchmark.sh`; `redesign/run_grouped_benchmark.R`; `redesign/run_full_recipe.R`; `redesign/run_validation_benchmark.R` | 3 outer repeats x 5 folds; seed 42; candidate pool 2,000 variable genes; panel sizes 10, 20, 50, 100, 200, 500; GeneSelectR `B=50`, five-fold subsampling, instance SHAP, stability + utility, geometric combination; 20 outcome permutations x 20 null fits | Reported |
| Saved-ranking reevaluation | `redesign/reevaluate_saved_rankings.R`; `redesign/validate_deterministic_results.R` | Reuses saved rankings, splits, and deterministic evaluator seeds | Reported |
| Post-ranking biological assessment | `redesign/run_older7_biology.R`; `redesign/run_older7_biology_ot_dense.R` | GO, Hallmark, Open Targets, and STRING; 1,000 matched random panels; dense Open Targets files for the seven manuscript datasets | Reported supplementary assessment |

## Ablation and sensitivity analyses

| Analysis | Entry point | Exact settings recorded in source | Status |
|---|---|---|---|
| Component ablation | `redesign/run_older7_component_ablation.R` | Saved outer splits; k = 10, 20, 50, 100, 200, 500; deterministic evaluator seeds; recurrence, SHAP, MI, and combined ranking variants | Ablation |
| Permutation adjustment on/off | `redesign/run_older7_calibration_onoff.R` | Same split, fit, and evaluator settings; adjusted and unadjusted score components | Ablation |
| Permutation-reference diagnostics | `redesign/run_older7_calibration_diagnostics.R` | Extract all 15 outer fits; optional r1f1 null check; leave-one-permutation-out diagnostics | Sensitivity analysis |
| Larger permutation reference | `redesign/run_older7_null100_confirmation.R` | 100 permutations for GSE107994 and GSE13355; first 20 rows must match the saved 20-permutation checkpoint | Sensitivity analysis |
| Random-reference sensitivity | `redesign/run_older7_random_sensitivity.R` | 30 matched random panels; first three reproduce the primary baseline; six panel sizes | Sensitivity analysis |
| Gene-set-size/adaptive analysis | `redesign/run_older7_adaptive.R` | Training-based adaptive gene-set-size selection; same older-seven split design | Sensitivity analysis |
| `glmnet` convergence replay | `redesign/run_older7_glmnet_convergence_check.R` | Replays saved model conditions and checks convergence diagnostics | Sensitivity analysis |

## Method-development configurations

The following experiments document alternatives. Their results do not define the
reported method.

| Configuration family | Entry points | Configuration details |
|---|---|---|
| Hallmark-restricted and grouped | `redesign/run_grouped_benchmark.R` | `var2000_ungrouped`, `var2000_grouped`, `bio_ungrouped`, `bio_grouped`, and `varlarge_ungrouped`; Hallmark pool and module gate; `gate_B=1000`; percentile calibration; B=50; seed 42 |
| Full recipe alpha comparison | `redesign/run_full_recipe.R` | Evidence-ratio calibration; 20 permutations x 20 null fits; alpha grid 0.5 and 1.0; internal OOB AUC selection |
| Prediction-first | `redesign/run_predfirst_benchmark.R` | Score ranking plus redundancy filter; raw score-only arm; panel sizes 10 to 500; B=50, five folds, alpha grid 0.5 and 1.0, gamma 0.25 |
| Soft prior | `redesign/run_softprior_benchmark.R` | Hallmark module z-scores; B=1000 label permutations; multiplier `1 + 0.5 * max(0, z)`; unannotated multiplier 1 |
| Data-driven modules | `redesign/run_datadriven_benchmark.R` | Correlation clusters with dynamic tree cut; module size 10 to 200; B=1000 label permutations; Hallmark and data-driven combined arms |
| Rank ensembles | `redesign/run_ensemble_eval.R`; `redesign/run_ensemble2_eval.R` | Rank-average combinations; `kswitch` cutoff 20; post-result `ens2` and `ens3` descriptive combinations |
| Weight sweep | `redesign/run_wsweep.R` | Combined soft prior at `w=2.0`; no model or null refit |
| Adaptive weight | `redesign/run_adapt_benchmark.R` | Training-only inner AUC; `W_GRID = 0, 0.5, 1, 2, 4`; `K_SEL=20`; arm `cbgs_adapt` |
| Horseshoe | `redesign/run_horseshoe_benchmark.R` | Horseshoe probit; 4,000 iterations; burn-in 1,000; seed 42; ranking by posterior mean absolute beta times PIP |
| Horseshoe stability | `redesign/run_hsstab_benchmark.R` | 50 five-fold subsamples; 2,000 iterations; burn-in 500; seed 42; two workers by default |
| Redundancy pruning | `redesign/run_prune_benchmark.R` | Greedy absolute-correlation filter `|cor| < 0.7`; combined Hallmark/data-driven boosted rankings |
| Compact and greedy panels | `redesign/run_gs_slim_imvigor210.R`; `redesign/run_greedy_imvigor210.R` | IMvigor210; preprocessing and evaluator matched to incumbent; GS_slim modes smoke/pilot/full; greedy pool 200, `k_max=50`, inner k=5 |
| T-Rex diagnostics | `redesign/run_trex_sanity.R`; `redesign/run_trex_diagnose.R`; `redesign/run_trex_diagnose2.R`; `redesign/run_group_probe.R` | IMvigor210 fold 1; T-Rex K and tFDR diagnostics; group probe `trex+GVS` and Hallmark permutation checks |

## Historical analyses

| Entry points | Classification |
|---|---|
| `benchmarks/sosall_benchmark.R`; `benchmarks/sosall_benchmark_v2.0.R`; `benchmarks/imvigor210_benchmark.R`; `imvigor210_benchmark_v2.0.R` | Historical benchmark implementations |
| `benchmarks/synthetic_benchmark.R` | Historical synthetic experiments |
| `benchmarks/compare_historical_primary.R`; `benchmarks/index_historical_results.R` | Historical result inspection |
| `redesign/run_biology_disease_benchmark.R` | Historical disease-specific biological assessment |

Historical scripts remain available for provenance. They are excluded from the
reported configuration.

## Exploratory and excluded from the paper

Dataset-role notes and the excluded entry points are also collected in
[`../exploratory_not_in_paper/README.md`](../exploratory_not_in_paper/README.md).

The following analyses use the four follow-up datasets or test method concepts
outside the reported paper scope:

| Entry points | Dataset or scope | Classification |
|---|---|---|
| `redesign/run_validation_benchmark.R` with `GSE16879`, `GSE91061`, `GSE92415`, or `GSE206285` | Four follow-up datasets | Exploratory; excluded from paper |
| `redesign/run_frozen_external_parallel.sh`; `redesign/run_frozen_external_competitors.R`; `redesign/run_frozen_external_adaptive.R`; `redesign/run_frozen_external_biology.R`; `redesign/summarise_frozen_external*.R` | Four follow-up datasets | Exploratory; excluded from paper |
| `redesign/run_random_baseline_sensitivity.R` | Four follow-up datasets | Exploratory; excluded from paper |
| `redesign/complementarity_analysis/05_exploratory_dge_prefilter.R` | All saved datasets, DGE prefilter | Exploratory; excluded from paper |

No result from this section is included by `supplementary/scripts/run_all.sh`.
