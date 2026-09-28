# Stabl biological enrichment on the seven manuscript datasets

The analysis scored the saved Stabl rankings in each of 15 outer-training divisions per dataset at 10, 20, 50, 100, 200, and 500 genes. It produced 630 split-size measurements. Gene selection was completed before this analysis. No predictive model or classical feature selector was rerun.

## Measurements

The five reported ratios are GO Biological Process semantic similarity, Hallmark shared-membership edges, STRING protein interaction edges, and Open Targets disease-association sums at score thresholds 0.05 and 0.10. Each observed value was divided by the mean from 1,000 matched random sets. The random sets and annotation snapshots follow `redesign/run_older7_biology.R` and `redesign/run_older7_biology_ot_dense.R`. STRING random sets were additionally matched on the number of mapped proteins.

GO, Hallmark, and Open Targets null means were reused from the original benchmark for the same dataset, outer division, set size, and source. Those references do not depend on the selection method. For STRING, 354 Stabl cells had a saved null mean for the same mapped-protein count. The remaining 276 cells used 1,000 new draws with the original seed formula and the saved candidate-pool order. Ratios are descriptive enrichment measures. No multiple-testing or independent-validation claim is made from them.

Dataset-level results use the median over 15 outer divisions separately at 10, 20, and 50 genes, followed by the mean of those three medians. `stabl_biology_by_size.csv` also reports all six sizes. The manuscript table pairs these results with the previously saved GeneSelectR 2.0 results; the asthma case study compares the median ratios at 20 genes.

## Files

- `*_stabl_biology.csv`: one row per outer division and gene-set size for one dataset.
- `stabl_biology_by_split.csv`: all 630 rows.
- `stabl_biology_by_size.csv`: median over 15 divisions for each of 42 dataset-size cells.
- `stabl_biology_by_dataset.csv`: primary 10/20/50 summary for each dataset.
- `reference_inputs/`: compact copies of the exact candidate-pool gene lists, biological references, dense Open Targets scores, Hallmark sets, and STRING v12 data used for scoring. `SHA256SUMS` records file hashes.
- `ranking_SHA256SUMS`: hashes of the 105 input Stabl ranking files in the parent directory.

The reference data are frozen copies of the benchmark candidate pools and annotation cache, stored in `reference_inputs`. The original benchmark source revision was `621be0c1` (GeneSelectR 2.0 version 0.99.2). This extension ran with R 4.6.1, igraph 2.3.3, withr 3.0.3, and msigdbr 26.1.1.

## Reproduce and verify

From the benchmark repository root, run `Rscript redesign/stabl_comparison/run_biology.R`, then `Rscript redesign/stabl_comparison/summarise_biology.R`. Both scripts write to this directory by default. Set `STABL_BIOLOGY_OUT` to an alternate directory that contains a copy of `reference_inputs` to write elsewhere.

The scoring code checks each Stabl ranking against its 2,000-gene outer-division candidate pool and validates all 630 output keys. Recalculation of GO, Hallmark, and STRING observed values for GeneSelectR 2.0 at 10, 20, and 50 genes in the first outer division of each dataset matched the saved benchmark values to less than `5e-15`. Repeating the asthma Stabl analysis from the archived references produced a byte-identical output CSV. The saved GeneSelectR 2.0 primary biology summaries matched independent aggregation of the frozen split-level results to less than `3e-14`.
