# Stabl comparison on the GeneSelectR 2.0 benchmark datasets

## Analysis

Stabl was evaluated on the seven datasets used in the manuscript: GSE101794, GSE107994, GSE13355, GSE65682, GSE69683, IMvigor210, and SOS-ALL. Each dataset used the 15 saved outer train–test divisions (three repetitions of five folds). Selection used the 2,000-gene variance-filtered candidate pool and standardized outer-training expression matrix saved for each division. Training means and standard deviations were applied to the corresponding outer-test samples. Psoriasis participants remained grouped during outer splitting and during Stabl subsampling.

The implementation is the authors' Stabl tag `v1.0.1-lw` (commit `2f048470eb1276c84310a3b7a6d53b6134f60205`). Its optional `knockpy` import was made conditional in the temporary installation because the run uses random-permutation artificial features; the selection algorithm was unchanged. The configuration follows the authors' cell-free RNA example: balanced L1 logistic regression, 150 subsamples per regularization value, 10 equally spaced `C` values from 0.01 to 1, 50% sample fraction without replacement, 50% random-permutation artificial features, and the data-derived Stabl threshold. For GSE13355, subsampling used the implementation's group-aware function.

To prepare the Python environment, check out the stated Stabl commit, apply `stabl_random_permutation_import.patch` from this directory, and install that checkout in a Python 3.11 environment with the versions below. The patch only makes the unused knockoff sampler import optional; random-permutation generation and selection code remain those of the tagged release.

For comparison at the manuscript's fixed sizes, genes were ranked by their maximum Stabl selection frequency over the 10 `C` values. Equal scores were resolved using the candidate-pool order, which was determined by outer-training expression variance. The top 10, 20, 50, 100, 200, and 500 genes were evaluated. These fixed-size panels are a benchmark adaptation of Stabl. The native Stabl threshold and selected-set size are recorded separately in the per-division metadata.

The original `glmnet` + `xgboost` + `ranger` soft-voting evaluator and the original evaluation seeds were applied to each Stabl panel. The baseline GeneSelectR 2.0 and classical-method AUCs were read from saved `eval_deterministic.csv` files; those methods were not fitted again. The original evaluator was checked by reproducing an existing GeneSelectR 2.0 AUC of 0.9745098 for GSE101794, outer division 1/1, and 10 genes.

AUC summaries average over 15 outer-test divisions within dataset and gene-set size. The primary validation summary averages sizes 10, 20, and 50 within each of the six independent datasets, then weights the six datasets equally. Sizes 100, 200, and 500 are reported separately. Nogueira stability uses the selected genes projected onto the intersection of the 15 candidate pools in each dataset, matching the manuscript's stated procedure.

## Biological enrichment

The saved Stabl gene sets were subsequently evaluated for GO-BP semantic similarity, Hallmark sharing, STRING interactions, and Open Targets disease associations at score thresholds 0.05 and 0.10. The analysis covers all seven datasets, 15 outer divisions, and all six panel sizes. Each set was compared with the same 1,000-set matched reference rules used for the original benchmark. The primary dataset summary is the mean of the size-specific medians at 10, 20, and 50 genes. All six size-specific results and the 630 split-size rows are retained in `biology/`.

Run `Rscript redesign/stabl_comparison/run_biology.R` and then `Rscript redesign/stabl_comparison/summarise_biology.R` from the repository root. `biology/reference_inputs/` contains compact frozen copies of the candidate pools and annotation inputs used for scoring. The full method, checks and output files are described in `biology/README.md`. Figure 4E of the manuscript includes the 20-gene asthma panels from Stabl.

## Running the comparison

`run_comparison.R` extracts one saved outer division at a time, runs `rank_stabl.py`, and evaluates all six panel sizes with the common evaluator. `summarise.R` combines the Stabl output with the saved benchmark results and computes the AUC and stability tables. Both scripts read the saved outer splits and candidate pools produced by the primary benchmark (`supplementary/scripts/run_primary.sh`).

From the repository root, set `GENESELECTR_STABL_PYTHON` to a Python 3.11 environment containing the tagged Stabl installation, and set `GENESELECTR_R_LIB` if `xgboost` and `ranger` are in a nonstandard R library. Then run `Rscript redesign/stabl_comparison/run_comparison.R` followed by `Rscript redesign/stabl_comparison/summarise.R`, or `bash supplementary/scripts/run_stabl_comparison.sh`. The scripts write rankings, per-division metadata, AUCs and summaries to this directory.

MVFS-SHAP was not evaluated because no executable implementation linked to the publication was identified.

Sources: [Stabl implementation](https://github.com/gregbellan/Stabl), [Stabl publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC10002850/), [MVFS-SHAP publication record](https://pubmed.ncbi.nlm.nih.gov/41289809/).

## Software versions used in this run

The original benchmark manifest records R 4.5.2. This extension used R 4.6.1 with the same recorded evaluator-package versions: `glmnet` 5.0, `xgboost` 3.2.1.1, `ranger` 0.18.0, and `withr` 3.0.3. Python 3.11; NumPy 1.26.4; SciPy 1.11.4; scikit-learn 1.3.2; pandas 2.1.4. The Stabl package identifies itself as 1.0.0 at tag `v1.0.1-lw`.

Two `glmnet` warnings about incomplete convergence at late lambda values were emitted during the final SOS-ALL evaluator pass. All panel predictions and AUCs were finite.
