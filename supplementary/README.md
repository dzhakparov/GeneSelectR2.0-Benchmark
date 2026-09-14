# GeneSelectR supplementary analyses

## Scope

This directory indexes the analyses that support the GeneSelectR paper. The
reported predictive analysis uses SOS-ALL for development and six validation
datasets:

| Role | Datasets |
|---|---|
| Development | SOS-ALL (`sosall`) |
| Validation | GSE101794, GSE107994, GSE13355, GSE65682, GSE69683, IMvigor210 (`imvigor210`) |
| Exploratory follow-up | GSE16879, GSE91061, GSE92415, GSE206285 |

The four exploratory datasets are indexed under `exploratory_not_in_paper/`
in the configuration registry. Their result files are excluded from the paper
supplement.

Targeted-assay, p009, and unrelated multi-omics analyses are outside this
scope.

## Navigation

| Path | Purpose |
|---|---|
| `methods/README.md` | Analysis methods, resampling, null references, evaluation, and biological assessment |
| `configurations/registry.md` | Configuration and classification registry for reported, sensitivity, development, historical, and exploratory work |
| `scripts/` | Repository-relative entry points that call the maintained runners under `redesign/` and `benchmarks/` |
| `figures/README.md` | Rendered figure files, captions, source tables, and verification state |
| `results_summary.md` | Numerical primary and component-ablation results with provenance |
| `provenance/README.md` | Package revision decision, script provenance, and reproducibility checklist |

The established implementation remains under `redesign/` and `benchmarks/`.
The scripts here provide a stable reviewer-facing entry point and preserve the
existing output locations.

## Recommended execution order

Run commands from the repository root.

1. `bash supplementary/scripts/check_inputs.sh`
2. `bash supplementary/scripts/run_tests.sh`
3. `bash supplementary/scripts/run_primary.sh`
4. `bash supplementary/scripts/run_method_development.sh`
5. `bash supplementary/scripts/run_biology.sh`
6. `bash supplementary/scripts/run_dge_complementarity.sh`
7. `bash supplementary/scripts/run_figures.sh`

The complete paper-scope sequence is:

```bash
bash supplementary/scripts/run_all.sh
```

The sequence uses two workers by default. Set
`GENESELECTR_N_CORES` before execution to change the worker limit. Set
`GENESELECTR_CALL_BUDGET` to change the per-run checkpoint budget.

The full sequence is computationally intensive. Primary seven-dataset fitting
and extension analyses can require many hours and substantial temporary disk
space. The default worker limit is two. Lightweight checks should be run before
any benchmark scope.

## Inputs

The required inputs are private or downloaded files under ignored paths. They
are listed by `bash supplementary/scripts/check_inputs.sh`, which delegates to
`analysis/check_inputs.R`.

Expected input groups include:

- SOS-ALL expression and metadata: `data/normalized_logcpm.csv` and
  `data/metadata.csv`.
- Prepared validation matrices and metadata under `data/GSE*/`.
- IMvigor210 data from `easierData` or `data/IMvigor210.all.rds`.
- Frozen annotation resources under `data/` for GO, Hallmark, Open Targets,
  and STRING assessments.
- The package source extracted into the ignored `package/` directory from the
  verified revision recorded in `provenance/README.md`.

Generated data, fitted objects, result tables, logs, caches, and reports remain
under ignored directories. They are not part of the repository commit.

## Analysis classes

The reported workflow is repeated five-fold elastic-net modelling, gene
recurrence, excluded-sample predictive contribution, outcome-permutation
references, equal-weight geometric-mean ranking, and downstream biological
assessment. Biological information is evaluated after the predictive ranking.

Sensitivity analyses quantify component removal, permutation-reference
variation, random-reference variation, gene-set-size choices, and convergence.
Method-development experiments record alternative priors, score combinations,
weight rules, stability procedures, pruning, compact panels, greedy selection,
and T-Rex diagnostics. Historical scripts preserve earlier implementations.
The registry identifies each class and states whether it contributes to the
paper supplement.

## Outputs

Primary and extension runners write to `redesign/results_corrected/`.
Complementarity analyses write to
`redesign/complementarity_analysis/results/` unless
`GENESELECTR_COMPLEMENTARITY_RESULTS_DIR` is set.

Figure generation writes PDFs, PNGs, and source CSV files under the ignored
complementarity result directory. The six paper-scope rendered figures are
committed under `figures/`. Their source tables remain external and are listed
in `figures/README.md`.

## Reproducibility limits

The repository contains scripts, rendered paper-scope figures, and a numerical
results summary. Required expression data, annotation caches, fitted objects,
and source result tables are unavailable in this clean worktree. The documented
workflow therefore provides a reproducible execution procedure after those
inputs are restored. It does not establish a fresh reproduction of numerical
paper results from this worktree alone.
