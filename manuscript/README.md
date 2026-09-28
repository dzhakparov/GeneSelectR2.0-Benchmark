# Manuscript source tables and figures

This directory contains the aggregated results behind the manuscript tables and
figures and the scripts that render the figures. The figures can be regenerated
without the raw expression data:

```bash
python3 make_figures.py          # Figures 2–4 and the enrichment-count figure (NumPy, pandas, Matplotlib)
python3 make_workflow_figure.py  # Figure 1 (NumPy, Matplotlib, cairosvg)
```

`make_figures.py` checks the expected numbers of datasets, methods and gene-set
sizes before drawing, and rewrites the derived summaries
`asthma_k20_biology_summary.csv`, `biology_enrichment_counts.csv` and
`biology_dataset_summary.csv`.

## Tables and figures

| Manuscript item | Files in `source_data/` | Contents |
|---|---|---|
| Table 2 | `prediction_table2_verified.csv` | SOS-ALL, six-validation and seven-dataset mean AUC per method, and the difference from the 30-panel matched random reference |
| | `random30_by_dataset_size.csv` | Matched random-panel reference AUC per dataset and gene-set size (30 panels per outer division) |
| Figure 2 | `by_dataset_primary.csv`, `by_validation_primary.csv` | Mean outer-test AUC over 10, 20 and 50 genes, per dataset and across the six validation datasets |
| | `by_dataset_size.csv`, `by_validation_size.csv` | The same at all six gene-set sizes (10–500) |
| Figure 3 | `stability_classical.csv` | Nogueira stability per dataset, method and gene-set size for the classical methods and GeneSelectR 2.0 |
| | `stability_by_dataset_primary.csv` | Nogueira stability for GeneSelectR 2.0 and Stabl, averaged over 10, 20 and 50 genes |
| Section 4.3 | `biology_comparison_by_dataset.csv` | Dataset-level GO-BP, Hallmark, STRING and Open Targets (0.05, 0.10) enrichment ratios for all eight methods |
| | `biology_by_dataset.csv`, `ot_by_dataset.csv`, `ot_dense_by_split.csv` | GO-BP, Hallmark and STRING ratios for the seven non-Stabl methods; Open Targets ratios per dataset and per outer division |
| | `stabl_biology_by_split.csv`, `stabl_biology_by_size.csv`, `stabl_biology_by_dataset.csv` | Stabl biological ratios per outer division, per size and per dataset |
| | `biology_enrichment_counts.csv` | Number of datasets with an enrichment ratio above one, per method and source |
| Figure 4 | `asthma_eight.csv` | Selection counts, recurrence, utility and Open Targets scores for the eight most recurrent GeneSelectR 2.0 genes in GSE69683 |
| | `asthma_biology.csv` | GO-BP, Hallmark and STRING ratios per outer division for the asthma gene sets |
| | `asthma_k20_biology_comparison.csv` | Median 20-gene enrichment ratios for the eight methods (panel E) |
| Figure 1 | — | Drawn by `make_workflow_figure.py` |

Dataset-level biological ratios are the mean of the median ratios over the 15
outer divisions at 10, 20 and 50 genes. Enrichment ratios are observed values
divided by the mean of 1,000 matched random gene sets.

## Notes on the tables

- Method labels: `GS_full_ungrouped` is GeneSelectR 2.0, `RF_importance` is
  random-forest importance and `ElasticNet` is elastic net.
- `stability_classical.csv` also contains four follow-up datasets (GSE16879,
  GSE91061, GSE92415, GSE206285) that are not part of the manuscript. Its
  `dataset_group` column uses an earlier labelling in which all seven
  manuscript datasets are marked `development`; the manuscript roles are given
  in Table 1.
- The `open_targets_*` columns in `asthma_biology.csv` come from an earlier
  Open Targets query and are not used. The Open Targets ratios in Figure 4E
  are taken from `ot_dense_by_split.csv`.
- The training-based gene-set size selection (Section 4.4) and the component
  analysis (Section 4.5) are produced by
  `supplementary/scripts/run_gene_set_size_sensitivity.sh`,
  `run_component_ablation.sh` and `run_permutation_adjustment.sh`.
  Per-dataset component results are listed in
  [`supplementary/results_summary.md`](../supplementary/results_summary.md).
