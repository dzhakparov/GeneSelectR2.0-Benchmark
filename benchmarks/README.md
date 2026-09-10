# Benchmark script index

## Current dataset preparation

- `validation_datasets.R` records cohort definitions, outcomes, exclusions,
  grouping variables, preprocessing scales and biological terms.
- `validation_prepare.R` creates analysis-ready validation inputs.
- `gse101794_build_expression.R` assembles and checks the per-sample
  GSE101794 TPM files.
- `count_preprocessing.R` contains the count-data preprocessing checks.
- `worker_budget.R` limits concurrent workers across benchmark processes.

## Validation and external-holdout analyses

- `validation_benchmark.R` runs validation-dataset feature selection and
  predictive evaluation.
- `lock_validation_panels.R` records gene sets before an external holdout is
  evaluated.
- `validation_holdout.R` evaluates previously recorded gene sets.
- `independent_preflight.R` checks cohort availability and study design.
- `independent_analyse_results.R` summarizes completed independent runs.

## Historical analyses

`sosall_benchmark.R`, `sosall_benchmark_v2.0.R`,
`imvigor210_benchmark.R`, the repository-root
`imvigor210_benchmark_v2.0.R`, and
`synthetic_benchmark.R` preserve earlier benchmark implementations. They are
retained to document method development and are not the reported workflow.
`compare_historical_primary.R` and `index_historical_results.R` inspect these
older runs when their ignored output directories are present.

Targeted-assay and p009 multi-omics analyses are excluded from this branch.
