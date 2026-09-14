# Exploratory analyses excluded from the paper

This section records the four follow-up datasets and their entry points:

- GSE16879
- GSE91061
- GSE92415
- GSE206285

The frozen-external workflow is maintained under `redesign/` for method
development and follow-up evaluation. Its outputs are excluded from the paper
supplement and are not called by `supplementary/scripts/run_all.sh`.

Relevant entry points are:

- `redesign/run_validation_benchmark.R` for preparation and evaluation;
- `redesign/run_frozen_external_parallel.sh`;
- `redesign/run_frozen_external_competitors.R`;
- `redesign/run_frozen_external_adaptive.R`;
- `redesign/run_frozen_external_biology.R`;
- `redesign/run_random_baseline_sensitivity.R`;
- `redesign/complementarity_analysis/05_exploratory_dge_prefilter.R`.

No exploratory result file is committed.
