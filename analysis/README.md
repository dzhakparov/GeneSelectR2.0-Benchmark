# Analysis navigation

## Scope

The code covers the complete recorded GeneSelectR analysis. SOS-ALL is the
development dataset. GSE101794, GSE107994, GSE13355, GSE65682, GSE69683 and
IMvigor210 are the six validation datasets used in the manuscript. GSE16879,
GSE91061, GSE92415 and GSE206285 are retained as exploratory follow-up
datasets and are kept separate from the manuscript analysis.

No result table, fitted object, cache, log, spreadsheet or raw expression file
is tracked on this branch. Every output is generated under ignored directories.
The `older7` prefix is retained in established script and output names. It
refers to the seven manuscript benchmark datasets: SOS-ALL for development and
six datasets for validation.

## Start here

1. Run `bash analysis/bootstrap_package.sh` to obtain the package source from
   the pinned package revision.
2. Run `Rscript analysis/check_inputs.R` to list missing inputs and software.
3. Run `bash analysis/check_repository_contents.sh` to verify that the branch
   tracks only code, documentation and optional images.
4. Run `bash analysis/run_all.sh tests` for the code checks.
5. Run the required analysis scope. `bash analysis/run_all.sh all` executes the
   full sequence and can require many hours.

The default is two workers. Set `GENESELECTR_N_CORES` to a safe value for the
machine before running the benchmarks. On macOS,
`OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES` must be present in `~/.Renviron`
before R is launched when more than one worker is used.

The package bootstrap uses the exact package revision recorded in
`analysis/config.R`. Set `GENESELECTR_PACKAGE_REF` only when intentionally
testing another package revision.

## Directory map

| Location | Contents |
|---|---|
| `analysis/config.R` | Dataset roles, shared settings and biological annotation identifiers |
| `analysis/run_all.sh` | Ordered entry point for the complete analysis |
| `analysis/check_inputs.R` | Input and package checks before a run |
| `redesign/R/` | Shared fitting, evaluation, competitor and biological functions |
| `redesign/run_corrected_benchmark.sh` | Primary seven-dataset benchmark |
| `redesign/run_older7_*` | Ablation, permutation, random-reference, size-selection and biological analyses |
| `redesign/complementarity_analysis/` | Stability and DGE/GeneSelectR complementarity analyses |
| `redesign/tests/` | Unit and integration checks for the analysis code |
| `benchmarks/` | Dataset preparation, validation and historical comparison scripts |

The available scopes are `tests`, `primary`, `ablations`, `biology`,
`complementarity`, `exploratory`, `review` and `all`. The `exploratory` scope
contains the four follow-up datasets and writes to its own result directory.

`fetch_geo_data.py`, `benchmarks/validation_prepare.R` and
`benchmarks/gse101794_build_expression.R` prepare the public GEO inputs. The
SOS-ALL expression and metadata files are study inputs. IMvigor210 is loaded
through `easierData` when available, with `data/IMvigor210.all.rds` as the
local fallback. All downloaded and prepared inputs remain under the ignored
`data/` directory.

## Configuration groups

### Reported configuration

The reported ranking uses repeated elastic-net fitting, 50 five-fold
resamples, within-training recurrence, excluded-sample predictive
contribution, 20 outcome permutations with 20 null fits, and a geometric mean
of the two adjusted components. Biological information is assessed after
ranking. Candidate genes are restricted to the 2,000 most variable genes
within each training division.

Primary scripts:

- `redesign/run_grouped_benchmark.R`
- `redesign/run_full_recipe.R`
- `redesign/run_validation_benchmark.R`
- `redesign/reevaluate_saved_rankings.R`
- `redesign/run_corrected_benchmark.sh`

### Feature-selection configurations

- Hallmark-restricted and grouped variants: `run_grouped_benchmark.R`
- Prediction-first variants: `run_predfirst_benchmark.R`
- Soft biological prior: `run_softprior_benchmark.R`
- Data-driven Hallmark modules: `run_datadriven_benchmark.R`
- Score ensembles: `run_ensemble_eval.R`, `run_ensemble2_eval.R`
- Weight sweep: `run_wsweep.R`
- Horseshoe and stability variants: `run_horseshoe_benchmark.R`,
  `run_hsstab_benchmark.R`
- Pruning and adaptive weights: `run_prune_benchmark.R`,
  `run_adapt_benchmark.R`
- Compact and greedy alternatives: `run_gs_slim_imvigor210.R`,
  `run_greedy_imvigor210.R`
- T-Rex and group-selection diagnostics: `run_trex_sanity.R`,
  `run_trex_diagnose.R`, `run_trex_diagnose2.R`, `run_group_probe.R`

These configurations are retained as method-development analyses. The reported
configuration is identified explicitly in the manuscript and package.

### Ablation and sensitivity analyses

- Component ablation: `run_older7_component_ablation.R`
- Permutation adjustment on/off: `run_older7_calibration_onoff.R`
- Permutation-null diagnostics: `run_older7_calibration_diagnostics.R`
- Larger permutation reference: `run_older7_null100_confirmation.R`
- Repeated random-gene reference: `run_older7_random_sensitivity.R`
- Training-based gene-set size selection: `run_older7_adaptive.R`
- `glmnet` convergence replay: `run_older7_glmnet_convergence_check.R`
- Combined summaries and uncertainty: `summarise_older7_*.R`

### Biological assessment

- GO, Hallmark and STRING summaries: `run_older7_biology.R`
- Dense Open Targets assessment: `run_older7_biology_ot_dense.R`
- Open Targets retrieval: `fetch_ot_seeds_dense.R`
- Disease-specific historical assessment: `run_biology_disease_benchmark.R`

The disease identifiers and GO terms are defined in `analysis/config.R`.
Biological information does not enter the reported predictive ranking.

### Complementarity and stability

Run the numbered scripts in `redesign/complementarity_analysis/` from 01 to 06.
They calculate gene-set stability, DGE/GeneSelectR agreement, gene-group
characteristics, biological annotations and the exploratory DGE prefilter.

## Generated directories

The scripts create `data/`, `cache/`, `redesign/results*/`,
`independent_benchmark_runs/`, `locked_panels/` and analysis-specific result
directories. These locations are ignored by Git. Each runner writes its own
configuration or provenance record into its output directory.
