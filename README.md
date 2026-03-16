# GeneSelectR 2.0 — Benchmark

Reproducible benchmarking scripts for **GeneSelectR 2.0**, a stability-based gene selection framework for RNA-seq data. These scripts accompany the manuscript and evaluate GeneSelectR 2.0 against established feature-selection methods on both synthetic and real-world transcriptomic datasets.

## Repository structure

```
├── benchmarks/
│   ├── synthetic_benchmark.R       # Synthetic RNA-seq validation (Experiments 1–5)
│   ├── sosall_benchmark.R          # SOS-ALL atopic dermatitis cohort benchmark
│   └── imvigor210_benchmark.R      # IMvigor210 cancer immunotherapy benchmark
├── data/                           # Input data directory (not tracked; see below)
├── .gitignore
└── README.md
```

## Benchmarks overview

### 1. Synthetic validation (`synthetic_benchmark.R`)

Evaluates GeneSelectR's gene recovery accuracy on simulated RNA-seq data with known ground truth. Data is generated via negative binomial count simulation with log₂(CPM + 1) normalization, heterogeneous fold changes, and optional co-expression modules.

**Experiments:**

| # | Question | Varies | Fixed |
|---|----------|--------|-------|
| 1 | Signal strength | Mean \|log₂FC\| ∈ {0.5, 1.0, 1.5, 2.5} | 20K genes, 100 DE, n = 100 |
| 2 | Sample size scaling | n ∈ {30, 50, 100, 200} | 20K genes, 100 DE, \|log₂FC\| = 1.5 |
| 3 | Signal sparsity | DE genes ∈ {50, 100, 200, 500} | 20K genes, n = 100, \|log₂FC\| = 1.5 |
| 4 | LASSO vs elastic net | α × correlation regime | 20K genes, 100 DE, n = 100 |
| 5 | Co-expression modules | Module count × within-module r | 20K genes, 100 DE, n = 100 |

**Metrics:** Recall, precision, F1 at multiple top-k cutoffs; CV AUC; selection stability (π).

**Output:** `results_realistic/` containing per-experiment CSVs, timing logs, individual experiment figures, and a combined 3×2 supplementary figure (PDF + PNG).

### 2. SOS-ALL atopic dermatitis (`sosall_benchmark.R`)

Nested cross-validation benchmark on the SOS-ALL cohort comparing nine feature-selection methods with hyperparameter optimisation for every method. Includes confounder residualisation (tissue location) and a variance pre-filter.

**Methods compared:**

- GeneSelectR (targeted — supervised GO terms for immune response)
- GeneSelectR (data-driven — enrichment-based biological scoring)
- Differential expression (t-test, BH-adjusted)
- LASSO
- Elastic net (α tuned via inner CV)
- mRMR (minimum redundancy, maximum relevance)
- Boruta (wrapper around random forest)
- Random forest variable importance
- Random baseline

**Evaluation:** 5×5 repeated stratified outer CV; for each fold, genes are ranked by each method, then the top-k genes (k ∈ {10, 20, 50, 100, 200, 500}) are used to train a tuned glmnet classifier and evaluated on the held-out fold. Pairwise Wilcoxon tests (BH-corrected) at k = 50 and k = 200.

**Output:** `results_targeted_sosall/<date>/` with nested results, summary tables, parsimony curves, AUC box plots, pairwise statistical tests, and a run configuration RDS.

### 3. IMvigor210 cancer immunotherapy (`imvigor210_benchmark.R`)

Same nested CV framework applied to the Mariathasan et al. (2018) metastatic urothelial carcinoma dataset. Binary outcome: responder (CR/PR) vs non-responder (SD/PD) to anti-PD-L1 (atezolizumab). Targeted GO terms include immune response and TGFβ signalling pathways.

**Data access** (auto-detected by the script):

- **Option A — easierData** (Bioconductor, maintained, ~192 patients):
  ```r
  BiocManager::install("easierData")
  ```
- **Option B — GitHub RDS** (full cohort, ~348 patients):
  ```r
  download.file(
    "https://github.com/snijeshvp/IMvigor210/raw/main/IMvigor210.all.rds",
    destfile = "data/IMvigor210.all.rds"
  )
  ```

**Output:** `results_cancer/<date>/` with the same structure as the SOS-ALL benchmark plus a formatted performance table.

## Requirements

**R ≥ 4.1** with the following packages:

Core (all scripts):
`GeneSelectR`, `glmnet`, `ggplot2`, `dplyr`, `tidyr`

Synthetic benchmark additionally:
`parallel`, `patchwork`

Real-data benchmarks additionally:
`mRMRe`, `Boruta`, `randomForest`, `ranger`

IMvigor210 benchmark additionally:
`edgeR`, `ExperimentHub`, `easierData`, `SummarizedExperiment` (for Bioconductor data access)

Missing optional packages are installed automatically at runtime.

## Data setup

### SOS-ALL cohort

Place the following files in `data/`:

- `normalized_logcpm.csv` — log₂ CPM expression matrix (rows = genes, columns = samples)
- `metadata.csv` — sample metadata with a `treatment` column formatted as `location_diagnosis` (e.g., `skin_AD`)

### IMvigor210

No manual setup needed — the script auto-detects installed sources. For fastest setup, install `easierData` via Bioconductor.

### Synthetic

No input data required; all data is generated programmatically.

## Running

```bash
# Synthetic validation (~hours depending on hardware)
Rscript benchmarks/synthetic_benchmark.R

# SOS-ALL benchmark
Rscript benchmarks/sosall_benchmark.R

# IMvigor210 benchmark
Rscript benchmarks/imvigor210_benchmark.R
```

All scripts auto-detect available CPU cores (using `n - 1` cores). Output directories are created automatically.

## Reproducibility

All scripts set `set.seed(42)` at the top. Stratified CV folds are seeded per repeat to ensure identical splits across runs. Date-stamped output directories (real-data benchmarks) allow preserving results from multiple runs.

## Citation

If you use these benchmarking scripts, please cite the GeneSelectR 2.0 manuscript (reference to be added upon publication).

## License

These scripts are provided as supplementary material to the GeneSelectR 2.0 manuscript. See the [GeneSelectR package repository](https://github.com/dgeorgiou3/GeneSelectR) for package licensing.
