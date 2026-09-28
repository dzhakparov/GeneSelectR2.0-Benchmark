# Numerical results summary

## Provenance and scope

This summary records numerical results for the seven manuscript datasets,
transcribed from the verified deterministic benchmark and extension reports.

> **Relation to the manuscript.** The manuscript compares the locked GeneSelectR 2.0 configuration
> (`GS_full_ungrouped`) with seven comparators: Stabl, mRMR, random-forest importance, DGE, Boruta,
> LASSO and elastic net. This document also covers method-development variants (for example `soft_gs`,
> `cbgs_prune`, `cb_gs_w2`, `GS_full_grouped`) that were evaluated during development and are not part
> of the manuscript comparison. Their status is listed in the
> [configuration registry](configurations/registry.md). Source tables for the manuscript results are in
> [`manuscript/source_data/`](../manuscript/source_data/).


The primary benchmark contains seven datasets, 27 unique methods, 15 outer
splits, six panel sizes, and 2,520 evaluation rows per dataset. Independent
reconstruction of one saved-ranking result per dataset reproduced the written
AUC values with a maximum absolute difference of `4.4e-16`.

## Primary performance

The primary metric is mean AUC minus matched Random across panel sizes 10, 20,
and 50. Values are descriptive summaries over repeated cross-validation
splits.

| Dataset | Highest-ranked method | Primary delta |
|---|---|---:|
| Crohn's disease, GSE101794 | mRMR | +0.045 |
| Tuberculosis, GSE107994 | soft_gs | +0.062 |
| Psoriasis, GSE13355 | cbgs_prune | +0.004 |
| Sepsis, GSE65682 | mRMR | +0.081 |
| Asthma, GSE69683 | GS_full_grouped | +0.110 |
| IMvigor210 | cb_gs_w2 | +0.134 |
| SOS-ALL | cbgs_prune | +0.128 |

GeneSelectR-derived configurations ranked first in five of seven datasets.
mRMR ranked first in two datasets. `GS_full_ungrouped` had the best
cross-dataset result among the two base GeneSelectR configurations, with mean
dataset rank 10.29 and mean primary delta +0.063. The deterministic correction
changed individual primary deltas by a mean absolute value of 0.0035 and
dataset ranks by a mean absolute value of 1.42.

## Component ablation

The table reports mean AUC minus the 30-draw matched Random mean for the fixed
current GeneSelectR recipe and four component variants. Panel sizes are 10, 20,
and 50. No variant was consistently superior across datasets.

| Dataset | k | Recurrence | SHAP | MI | SHAP x MI | Current |
|---|---:|---:|---:|---:|---:|---:|
| GSE101794 | 10 | 0.058 | 0.060 | 0.055 | 0.055 | 0.055 |
| GSE101794 | 20 | 0.038 | 0.035 | 0.032 | 0.038 | 0.038 |
| GSE101794 | 50 | 0.020 | 0.020 | 0.014 | 0.020 | 0.022 |
| GSE107994 | 10 | 0.087 | 0.087 | 0.086 | 0.087 | 0.084 |
| GSE107994 | 20 | 0.045 | 0.039 | 0.047 | 0.041 | 0.046 |
| GSE107994 | 50 | 0.018 | 0.018 | 0.015 | 0.017 | 0.015 |
| GSE13355 | 10 | -0.000 | 0.004 | 0.005 | 0.005 | -0.001 |
| GSE13355 | 20 | 0.001 | 0.002 | 0.002 | 0.002 | 0.001 |
| GSE13355 | 50 | 0.001 | 0.001 | 0.002 | 0.002 | 0.001 |
| GSE65682 | 10 | 0.060 | 0.079 | 0.049 | 0.073 | 0.071 |
| GSE65682 | 20 | 0.051 | 0.057 | 0.041 | 0.075 | 0.076 |
| GSE65682 | 50 | 0.043 | 0.031 | 0.030 | 0.044 | 0.045 |
| GSE69683 | 10 | 0.096 | 0.096 | 0.082 | 0.099 | 0.098 |
| GSE69683 | 20 | 0.079 | 0.071 | 0.067 | 0.070 | 0.076 |
| GSE69683 | 50 | 0.032 | 0.027 | 0.038 | 0.034 | 0.035 |
| imvigor210 | 10 | 0.110 | 0.070 | 0.122 | 0.136 | 0.126 |
| imvigor210 | 20 | 0.101 | 0.083 | 0.115 | 0.124 | 0.124 |
| imvigor210 | 50 | 0.086 | 0.069 | 0.110 | 0.072 | 0.094 |
| sosall | 10 | 0.103 | 0.032 | 0.124 | 0.097 | 0.097 |
| sosall | 20 | 0.072 | 0.024 | 0.093 | 0.070 | 0.066 |
| sosall | 50 | 0.073 | 0.043 | 0.075 | 0.043 | 0.051 |

The fixed current recipe was within 0.01 of the best component in 17 of 21
dataset-by-panel-size cells. The strongest component varied by dataset.

## Validation checks

The extension analyses were reconstructed from saved outer-training fits.

| Analysis | Rows or checks per dataset | Result |
|---|---:|---|
| Component ablation | 450 evaluation; 37,500 panel; 30 summary; 90 QA | All QA passed; maximum AUC reconstruction difference `5.6e-16` |
| Calibration diagnostics | 30,000 gene-score; 15 dependence; 4 null-validation; 15 QA | All QA passed; saved ratios reproduced with maximum difference 0 |
| Random-baseline sensitivity | 2,700 draw; 90 by-split; 6 summary; 90 QA | All QA passed; first three draws reproduced with maximum difference `5.6e-16` |

All findings are descriptive. Repeated cross-validation splits share
observations, so significance claims are not made in this summary.
