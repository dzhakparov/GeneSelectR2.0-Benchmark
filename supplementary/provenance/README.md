# Supplementary provenance

## Repository and package revisions

The repository source is `codex/analysis-archive` at commit `13314a54`.
The requested package reference is `621be0c1`, recorded in `analysis/config.R`
and used by `analysis/bootstrap_package.sh`.

Revision `db1f6528` was checked against `621be0c1` on 2026-09-14. The
comparison used the formal arguments of `geneselectr2_fit` from each package
revision. `621be0c1` provides the interface used by the maintained runners:
`gate_method`, `subsample_scheme`, `subsample_k_folds`, `utility_method`,
`components`, `score_formula`, `calibration_mode`, `calibration_n_permutations`,
`calibration_null_B`, `use_cache`, and `n_cores`. `db1f6528` provides the
simplified interface with `B`, `k_folds`, `permutations`, `null_B`, and
`workers`; the runner arguments are absent. The cleaned package revision is
therefore incompatible with the archived reported workflow.

The comparison was a call-interface check. Numerical reproduction of the
reported calculation was unavailable because the clean worktree contains no
expression matrices, metadata, annotation caches, or saved result files.
`621be0c1` remains the verified workflow reference. The revision check can be
repeated with:

```bash
bash supplementary/provenance/check_package_revision.sh
```

## Script provenance

The reviewer-facing wrappers call these maintained sources:

- `analysis/config.R` and `analysis/run_all.sh` for shared dataset roles and
  settings;
- `redesign/run_corrected_benchmark.sh`, `run_grouped_benchmark.R`,
  `run_full_recipe.R`, and `run_validation_benchmark.R` for the primary
  workflow;
- `redesign/run_older7_*.R` for the archived ablations, diagnostics, and
  biological analyses;
- `redesign/complementarity_analysis/01` through `06` for stability,
  complementarity, and figures;
- `benchmarks/` scripts for data preparation and historical comparisons.

Each maintained runner writes its own configuration or provenance record to its
ignored output directory when inputs are available.

## Reproducibility checklist

- [ ] Checkout `codex/paper-supplementary-materials`.
- [ ] Restore the required expression matrices, metadata, and annotation
      resources under ignored `data/` paths.
- [ ] Confirm `git rev-parse 621be0c1` succeeds.
- [ ] Run `bash supplementary/scripts/check_inputs.sh`.
- [ ] Run `bash supplementary/scripts/run_tests.sh`.
- [ ] Set `GENESELECTR_N_CORES` to a safe value; the default is two.
- [ ] Run the primary workflow and record its output provenance files.
- [ ] Run ablation and biological analyses after primary checkpoints exist.
- [ ] Run complementarity scripts 01 to 04 before script 06.
- [ ] Inspect every figure source table and rendered output.
- [ ] Run `bash analysis/check_repository_contents.sh`.
- [ ] Confirm `git status --short` contains only intended code and Markdown
      files.
- [ ] Preserve the exact package revision, seeds, worker limit, and input
      checksums in the run archive.

## Missing inputs in the current worktree

The following required inputs were absent during repository preparation:

- package source extracted under `package/GeneSelectR/`;
- SOS-ALL expression and metadata;
- all prepared validation matrices and metadata;
- IMvigor210 data;
- GO, Hallmark, Open Targets, and STRING resources;
- saved benchmark splits, rankings, fits, and result tables.

The missing inputs prevent numerical validation and figure acceptance in this
clean worktree.

The requested source document `redesign/full_benchmark_plan.md` is absent from
the archive branch. The maintained plan is distributed across
`analysis/README.md`, `analysis/config.R`, `analysis/run_all.sh`,
`benchmarks/README.md`, and the script headers under `redesign/`.
